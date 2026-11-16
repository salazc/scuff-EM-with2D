Analyzed CreateC2DWorkspace.cc:L151-320

`scuff-cas2DLCN` is a 2D boundary-element-method (BEM) Casimir solver for geometries that are translationally invariant along the $z$-axis (cross-sections defined in `.tdgeo` files and managed by `libTDRT`). It computes Casimir **energy** and **x/y forces** (optionally separated into **TE** and **TM** polarizations) across a list of geometrical transformations.

---

### 1. Input & Configuration ([scuff-cas2DLCN.cc](file:///home/salazc/scuff-EM-with2D/applications/scuff-cas2DLCN/scuff-cas2DLCN.cc))
`main()` reads an XML configuration file (see [sample_config.xml](file:///home/salazc/scuff-EM-with2D/applications/scuff-cas2DLCN/sample_config.xml)) containing:
- **Geometry & Transforms**: `<Geometry>` (`.tdgeo` file), `<TransList>` (list of displacements/rotations with tags), `<LengthUnit>`, `<GroundPlane>`, and optional `<Rectangle>x1 y1 x2 y2</Rectangle>` (to restrict force/energy trace sums to internal vertices inside a box).
- **Calculation Mode** (`<CalculationMode type="...">`):
  - `FullIntegral`: 0 K integral over the $(\xi, q)$ quarter-plane.
  - `SinglePoint`: Single $(\xi, q)$ evaluation.
  - `QIntegral`: Integral over $q \in [0, \infty)$ at fixed $\xi$.
  - `XiIntegral`: Integral over $\xi \in [0, \infty)$ at fixed $q$.
  - `MatsubaraSum`: Finite-temperature sum at `<T>` (in Kelvin; `T=0` falls back to `FullIntegral`).
  - `XQList`: Evaluates $(\xi, q)$ pairs from a file.
  - `VisualizeOnly`: Generates `gpmesh/` and `ppmesh/` visualization files and exits.
- **Outputs, Cubature & Runtime**: `<Energy>`, `<XForce>`, `<YForce>`, `<TETM>`, `<WriteHDF5>`, `<IntCache>`, `<ProfileFile>`, `<AbsTol>`, `<RelTol>`, `<XQMin>`, `<NumThreads>`, and `<LogLevel>`.

---

### 2. Workspace Initialization ([CreateC2DWorkspace.cc](file:///home/salazc/scuff-EM-with2D/applications/scuff-cas2DLCN/CreateC2DWorkspace.cc))
`CreateC2DWorkspace()` builds the `C2DWorkspace` ([scuff-cas2DLCN.h](file:///home/salazc/scuff-EM-with2D/applications/scuff-cas2DLCN/scuff-cas2DLCN.h)):
1. **Transform parsing & visualization**: Iterates over `TransList`, applies each transform to `TDRTGeometry *G`, writes `gpmesh/<base>.<tag>` and `ppmesh/<base>.pp`, and records `TransLines` and `Tags`.
2. **Ground-plane image objects** ([ImageObjects.cc](file:///home/salazc/scuff-EM-with2D/applications/scuff-cas2DLCN/ImageObjects.cc)): If `GroundPlane` is enabled, `CreateImageObject()` creates a mirrored copy of each object across $y = 0$ (rotations in `TransList` are forbidden in this mode).
3. **Matrix allocation**:
   - `T[no]` ($N_{\text{bf},a} \times N_{\text{bf},a}$, Hermitian): Self-interaction block per unique object (reused across identical mates via `G->Mate[no]`).
   - `TI[no]` (if `GroundPlane`): Self-image interaction block.
   - `Uab[no][nop]` ($b > a$): Off-diagonal coupling blocks between object $a$ and object $b$.
   - `M` ($N \times N$, `LHM_NORMAL` complex): Full BEM matrix (stored non-compressed for faster LAPACK linear algebra).
   - `dM` ($N \times N_{\text{bf},0}$) and `dU0bdX[nop]`, `dU0bdY[nop]`, `dT0dY`: Derivative blocks with respect to object 0's position (allocated only when forces are requested).
4. **Static SSI tables**: Precomputes frequency-independent segment-segment integral tables (`TSSSIDataTables` per unique object and `USSSIDataTables` per object pair and transform) to speed up matrix assembly at every $(\xi, q)$.

---

### 3. BEM Matrix & Physics Kernel ([XQIntegrand.cc](file:///home/salazc/scuff-EM-with2D/applications/scuff-cas2DLCN/XQIntegrand.cc))
Because the geometry is invariant in $z$, fields are Fourier-transformed along $z$ with wavenumber $q$ and evaluated at imaginary frequency $\omega = i\xi$. At each $(\xi, q)$ point, `XQIntegrand()` performs:

1. **Cache check**: Calls `CacheRead()` (when `WriteHDF5 == 0`) to look up $(\xi, q)$ in `<GeoBase>.byXQ` and skip recomputation if all transforms are present.
2. **External medium & $T$ blocks**: Queries $\varepsilon_{\text{out}}(i\xi), \mu_{\text{out}}(i\xi)$ via `G->MP->GetEpsMu(cdouble(0, Xi), ...)` and assembles `T[no]` for each unique object via `AssembleT()`.
3. **$M_\infty$ reference for energy**: If energy is requested, each unique `T[no]` is LU-factorized with LAPACK `zgetrf_` and the magnitudes of its diagonal entries are stored in `DRMInf` (representing the isolated-objects reference matrix $M_\infty$).
4. **Loop over transformations**:
   - Applies `TransLines[nt]` to `G` (and displaces `ImageObjects` to $(x, -y)$ if `GroundPlane` is active).
   - Assembles `TI[no]` (and `dT0dY` for object 0) via `AssembleTI()`, and `Uab[no][nop]` (plus `dU0bdX[nop]`, `dU0bdY[nop]` for object 0) via `AssembleU()` and `AddImageContributionToU()`—skipping blocks where neither object moved since the previous transform.
   - **Assembly & LU factorization** (`Factorize()`): Stamps `T[no]` on the diagonal, `Uab[no][nop]` and its transpose on the off-diagonals, and adds `TI[no]` to the diagonal blocks if `GroundPlane` is enabled, then calls `M->LUFactorize()`.
   - **Energy** (`GetLNDetMInvMInf()`): Computes
     $$\mathcal{E}(\xi, q) = -\frac{1}{2\pi^2} \ln\det(M^{-1} M_\infty) = -\frac{1}{2\pi^2} \sum_n \ln\left(\frac{\text{DRMInf}_n}{|M_{nn}|}\right)$$
   - **Force** (`GetTraceMInvdM()`): Stamps `dU0b` (and `dT0dY` for y-force with a ground plane) into `dM`, solves $M^{-1} \partial M$ via `M->LUSolve(dM)`, and sums the diagonal entries of object 0's block (or only `ContributingIVIndices` if `<Rectangle>` is set):
     $$\mathcal{F}_{X,Y}(\xi, q) = -\frac{2}{2\pi^2} \operatorname{Tr}_{\text{obj } 0}\!\left(M^{-1} \frac{\partial M}{\partial X,Y}\right)$$
     (The factor of $2$ accounts for the fact that object 0's block gives half the full trace.)
   - **TE/TM decomposition**: Basis functions are grouped per internal vertex—2 per vertex for `AllPEC` (even = TE, odd = TM) and 4 per vertex for dielectrics (`0, 3` = TE: $J_\parallel, M_z$; `1, 2` = TM: $J_z, M_\parallel$). Optional `SCUFF_ELECTRIC_ONLY=1` / `SCUFF_MAGNETIC_ONLY=1` env vars isolate electric or magnetic current contributions.
   - Writes the $(\xi, q)$ row to `<GeoBase>.byXQ` and undoes the transform.

---

### 4. Integration & Matsubara Summation ([CubatureMethods.cc](file:///home/salazc/scuff-EM-with2D/applications/scuff-cas2DLCN/CubatureMethods.cc))
- **`EvaluateQIntegral` / `EvaluateXiIntegral`**: Use `libSGJC`'s `pcubature` with $q = x/(1-x)$ and $\xi = \text{XQMin} + x/(1-x)$ (assuming a constant integrand on $[0, \text{XQMin}]$ for $\xi$).
- **`EvaluateXQIntegral`** (0 K `FullIntegral`):
  - For `AllPEC` geometries, exploits rotational symmetry in $(\xi, q)$ space (`RXQIntegrand`) to reduce the 2D integral to a 1D radial integral at $q=0$ with Jacobian $\frac{\pi}{2}\,\text{RXQ}/(1-x)^2$.
  - For general dielectric geometries, runs a 2D `pcubature_log` over $(\xi, q)$.
- **`EvaluateMatsubaraSum`**: Sums $2\pi k_B T \sum'_{n} \int dq\, I(\xi_n, q)$ with $\xi_n = 2\pi k_B T n$ ($k_B T = 4.36763\times 10^{-4}\,T$ in internal units, $n=0$ weighted by $0.5$ at $\xi = \text{XQMin}$), stopping once all quantities converge within tolerance for 3 consecutive terms.

---

### 5. Output Files
- `<GeoBase>.byXQ`: Raw integrand values at every evaluated $(\xi, q)$ point (also serves as the cache for `CacheRead()`).
- `<GeoBase>.byXi` / `<GeoBase>.byQ`: Intermediate 1D integral values at each $\xi$ or $q$.
- `<GeoBase>.out` (or `.out.N`): Final integrated/summed values and error estimates per transform tag.
- Optional `<ProfileFile>` (per-vertex TE/TM breakdown), `<GeoBase>.log`, and MATLAB/HDF5 matrix exports (`WriteHDF5`).