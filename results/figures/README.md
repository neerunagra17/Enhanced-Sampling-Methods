# Result figures

Alanine dipeptide in explicit water (ff14SB, TIP3P, 300 K). Interim results for the four completed methods: umbrella sampling (US), TASS, well-tempered metadynamics (WT-MetaD) and TAMD. Each method has 500 ns of production. WT-MetaD is shown with two estimators from the same trajectory:
- **bias:** sum of the deposited Gaussians, $F = -\sum(\text{stored hills})$;
- **reweighted:** reweighting with $V - c(t)$.

OPES results are pending.

**Conventions:**
- Grid and units: common 72 × 72 grid (5° bins); surfaces shifted to a minimum of zero and given in kJ/mol; $L^2$ errors in kcal/mol.
- Unvisited bins are white.
- $t_{\text{eff}}$ is the total production time of a method. For US and TASS it is $32 \times$ the time per window.

## `fes2d_500ns_all_methods.png`
**Free-energy surfaces $F(\phi,\psi)$ at $t_{\text{eff}} = 500$ ns.** Crosses mark the minima of the β, P<sub>II</sub>, α<sub>R</sub> and α<sub>L</sub> basins. All estimates show the same topology, with the global minimum in the P<sub>II</sub> basin.
- Umbrella sampling, restrained along $\phi$ only, leaves part of the region near the $\phi \approx 0$ barrier at negative $\psi$ unvisited.
- TAMD is noisier, and its surface is scaled up relative to the other methods.

## `fes1d_profiles_phi_psi_500ns.png`
**One-dimensional free-energy profiles along (a) $\phi$ and (b) $\psi$ at 500 ns.** Each profile is obtained by Boltzmann integration of the 2D surface over the other angle at 300 K.
- US, TASS and both WT-MetaD estimates nearly coincide.
- TAMD gives systematically higher barriers.

## `fes_evolution_20_100_500ns.png`
**Development of the surfaces with $t_{\text{eff}}$** (columns: 20, 100 and 500 ns; rows: methods and estimators). Each US and TASS window contributes the interval from 0.5 ns to $t_{\text{eff}}/32$; the first 0.5 ns of every window is discarded.
- The bias-derived WT-MetaD estimate is defined everywhere from the start, because it is a sum of Gaussians of width 0.35 rad.
- The window methods and the histogram-based estimates fill in progressively.

## `l2_internal_convergence_and_coverage.png`
**Internal convergence.**
- **(a)** $L^2$ error of each estimate relative to the same method's 500 ns surface, over the points with $F < 14.58$ kcal/mol, mean-aligned. The dotted line is $k_BT$.
- **(b)** Percentage of those points not yet visited; they are excluded from (a).

Internal convergence measures how fast an estimate stops changing, not its accuracy.

## `l2_pairwise_200ns_500ns.png`
**Pairwise $L^2$ errors between estimates at (a) 200 ns and (b) 500 ns** (kcal/mol, mean-aligned). Each row is the reference surface that defines the mask, so the matrix is not exactly symmetric.
- At 500 ns, US, TASS and WT-MetaD agree with one another to 0.28–0.35 kcal/mol.
- The two WT-MetaD estimators agree to 0.11 kcal/mol; they come from the same trajectory.
- TAMD differs from all other estimates by 0.67–1.11 kcal/mol.

These are single runs with one parameter set per method, without statistical uncertainties. Full discussion: UGP report, Chapter 6.
