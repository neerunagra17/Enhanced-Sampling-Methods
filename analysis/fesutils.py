"""fesutils.py - shared helpers for all analysis scripts.

Conventions used everywhere in this project
  * energies in kJ/mol (PLUMED's default unit), angles in rad
  * phi and psi are periodic on [-pi, pi)
  * every FES is written on the SAME regular grid of NBINS x NBINS bin centres
    in a 3-column text file  "phi  psi  F"  (blank line after each phi block,
    so gnuplot 'splot ... with pm3d' works); unsampled bins are written as nan
"""
import re
import numpy as np

KB = 0.0083144626          # Boltzmann constant, kJ/(mol K)
NBINS = 72                 # default grid: 72 x 72 bins = 5 degrees per bin
TWO_PI = 2.0 * np.pi


# --------------------------------------------------------------------------
# reading PLUMED COLVAR files
# --------------------------------------------------------------------------
def load_colvar(path):
    """Read a PLUMED COLVAR file into a dict {field_name: np.array}.

    Uses the '#! FIELDS' header, so columns are found by name, not position.
    If the file was appended after a restart (a second FIELDS line), the
    second block is appended to the first.
    """
    fields, rows = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#! FIELDS"):
                new = line.split()[2:]
                if fields is not None and new != fields:
                    raise ValueError(f"{path}: FIELDS changed inside the file")
                fields = new
                continue
            if line.startswith("#") or not line.strip():
                continue
            rows.append(line.split())
    if fields is None:
        raise ValueError(f"{path}: no '#! FIELDS' header found")
    ncol = len(fields)
    data = np.array([r[:ncol] for r in rows if len(r) >= ncol], dtype=float)
    cv = {name: data[:, i] for i, name in enumerate(fields)}
    if "time" in cv:
        cv["time"] = continuous_time(cv["time"])
    return cv


def continuous_time(t):
    """Make time monotonic across restarts.

    When a long run is done in segments, PLUMED may restart its clock at 0 in
    every segment while appending to the same COLVAR/HILLS file. Every time the
    clock jumps backwards, the rest of the series is shifted so that it
    continues one stride after the previous frame.
    """
    t = np.array(t, dtype=float)
    if len(t) < 2:
        return t
    dt = np.median(np.diff(t[: min(len(t), 1000)]))
    for k in np.where(np.diff(t) <= 0)[0]:
        t[k + 1:] += t[k] + dt - t[k + 1]
    return t


def select_time(cv, tmin=None, tmax=None, time_key="time"):
    """Keep frames with tmin <= time < tmax (times in ps). Returns a new dict."""
    t = cv[time_key]
    mask = np.ones_like(t, dtype=bool)
    if tmin is not None:
        mask &= t >= tmin
    if tmax is not None:
        mask &= t < tmax
    return {k: v[mask] for k, v in cv.items()}


# --------------------------------------------------------------------------
# grids, histograms and free energies
# --------------------------------------------------------------------------
def grid(nbins=NBINS):
    """Bin edges and bin centres on [-pi, pi)."""
    edges = np.linspace(-np.pi, np.pi, nbins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    return edges, centres


def wrap(x):
    """Map angles to [-pi, pi)."""
    return (np.asarray(x) + np.pi) % TWO_PI - np.pi


def hist2d(phi, psi, nbins=NBINS, weights=None):
    edges, _ = grid(nbins)
    H, _, _ = np.histogram2d(wrap(phi), wrap(psi), bins=[edges, edges],
                             weights=weights)
    return H


def fes_from_prob(P, kT):
    """F = -kT ln P, shifted so that min(F) = 0; empty bins -> nan."""
    F = np.full(P.shape, np.nan)
    ok = P > 0
    F[ok] = -kT * np.log(P[ok])
    F -= np.nanmin(F)
    return F


def log_weights_to_weights(logw):
    """exp(logw) without overflow (normalised so the largest weight is 1)."""
    logw = np.asarray(logw, dtype=float)
    return np.exp(logw - logw.max())


# --------------------------------------------------------------------------
# writing / reading FES files in the common format
# --------------------------------------------------------------------------
def write_fes(path, F, header=""):
    nb = F.shape[0]
    _, c = grid(nb)
    with open(path, "w") as fh:
        fh.write(f"# phi(rad) psi(rad) F(kJ/mol)  grid {nb}x{nb}  {header}\n")
        for i in range(nb):
            for j in range(nb):
                fh.write(f"{c[i]:10.5f} {c[j]:10.5f} {F[i, j]:12.5f}\n")
            fh.write("\n")


def read_fes(path, nbins=NBINS, energy_scale=1.0, xcol=0, ycol=1, fcol=2):
    """Read ANY 2D FES text file (ours, sum_hills, WHAM output, ...) and
    return it on the common nbins x nbins grid (min shifted to 0).

    The file may use a different grid: values are interpolated periodically.
    energy_scale converts units (use 4.184 for files in kcal/mol).
    Lines starting with '#' or '@' are skipped; inf/nan/huge values become nan.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in "#@":
                continue
            p = s.split()
            try:
                rows.append((float(p[xcol]), float(p[ycol]), float(p[fcol])))
            except (ValueError, IndexError):
                continue
    a = np.array(rows)
    x, y, f = wrap(a[:, 0]), wrap(a[:, 1]), a[:, 2] * energy_scale
    f[~np.isfinite(f)] = np.nan
    f[np.abs(f) > 1e5] = np.nan

    xs, ys = np.unique(np.round(x, 6)), np.unique(np.round(y, 6))
    Fin = np.full((len(xs), len(ys)), np.nan)
    ix = np.searchsorted(xs, np.round(x, 6))
    iy = np.searchsorted(ys, np.round(y, 6))
    Fin[ix, iy] = f

    _, c = grid(nbins)
    if len(xs) == nbins and len(ys) == nbins and np.allclose(xs, c, atol=1e-4) \
            and np.allclose(ys, c, atol=1e-4):
        Fout = Fin
    else:
        Fout = _interp_periodic(xs, ys, Fin, c)
    return Fout - np.nanmin(Fout)


def _interp_periodic(xs, ys, Fin, c):
    """Bilinear periodic interpolation of Fin(xs, ys) onto centres c x c.
    A target point is nan if any of its 4 neighbours is nan."""
    from scipy.interpolate import RegularGridInterpolator
    # pad one point on each side using periodicity
    xp = np.concatenate([[xs[-1] - TWO_PI], xs, [xs[0] + TWO_PI]])
    yp = np.concatenate([[ys[-1] - TWO_PI], ys, [ys[0] + TWO_PI]])
    Fp = np.pad(Fin, 1, mode="wrap")
    interp = RegularGridInterpolator((xp, yp), Fp, bounds_error=False,
                                     fill_value=np.nan)
    X, Y = np.meshgrid(c, c, indexing="ij")
    return interp(np.stack([X.ravel(), Y.ravel()], axis=1)).reshape(X.shape)


# --------------------------------------------------------------------------
# misc
# --------------------------------------------------------------------------
def parse_list(s):
    """'10,20,50' -> [10.0, 20.0, 50.0]"""
    return [float(v) for v in s.split(",") if v.strip()]


def teff_from_name(path):
    """Extract t_eff from a file name like 'fes_opes_120ns.dat' -> 120.0"""
    m = re.search(r"_([0-9]+(?:\.[0-9]+)?)ns", path)
    return float(m.group(1)) if m else np.nan


def fmt_ns(x):
    """120.0 -> '120', 12.5 -> '12.5' (for file names)"""
    return f"{x:g}"
