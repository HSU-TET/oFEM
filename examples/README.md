# oFEM examples

The examples are sorted first by the finite element space and then by the
kind of time dependence:

- **Static**: no time dependence (electrostatics, stationary currents,
  steady flow).
- **QuasiStatic**: time dependent, but without wave propagation: diffusion
  type problems (heat, diffusion, electro- and magneto-quasistatics, eddy
  currents) and flows that follow the driving force instantaneously.
- **FullWave**: the full Maxwell equations, i.e. wave propagation and
  resonators.

Run every script from its own folder (the mesh and geometry paths are
relative to it). Meshes (`*.msh`) and `export` folders are not tracked; most
scripts create a missing mesh with gmsh from the `.geo` file, otherwise run
gmsh on it yourself. Examples marked *Cubit* use meshes exported from
Cubit/Trelis (`.inp`, loaded with `load_from_inp`).

## H1 (`H1Element`)

| Example | Time | Description |
|---|---|---|
| [H1/Static/Capacitor2D](H1/Static/Capacitor2D) | static | 2D capacitor with two dielectrics; `tutorial.mlx` is a step by step introduction to oFEM based on it |
| [H1/Static/Capacitor3D](H1/Static/Capacitor3D) | static | 3D planar capacitor for several element orders |
| [H1/Static/CubitPlanarCapacitor2D](H1/Static/CubitPlanarCapacitor2D) | static | 2D planar capacitor, *Cubit* mesh |
| [H1/Static/CubitPlanarCapacitor3D](H1/Static/CubitPlanarCapacitor3D) | static | 3D planar capacitor, *Cubit* mesh |
| [H1/Static/CubitChargedSphere](H1/Static/CubitChargedSphere) | static | potential of a charged sphere, *Cubit* mesh |
| [H1/Static/PlanarCapacitor](H1/Static/PlanarCapacitor) | static | planar capacitor (`planarCapacitor.m`); `Main.m`: Poisson problem with Dirichlet and Neumann data and a volume force on the same geometry |
| [H1/Static/PlanarCapacitorCUDA](H1/Static/PlanarCapacitorCUDA) | static | planar capacitor solved with `pcg` on the GPU (CPU if no GPU is available) |
| [H1/QuasiStatic/Heat2D](H1/QuasiStatic/Heat2D) | quasi-static | transient heat conduction in a steel rod (`parabolic_solve`, order 3) |
| [H1/QuasiStatic/Diffusion2D](H1/QuasiStatic/Diffusion2D) | quasi-static | transient diffusion through paper and textile layers with Dirichlet and Neumann boundaries (`parabolic_solve`) |
| [H1/QuasiStatic/RobinBoundary](H1/QuasiStatic/RobinBoundary) | quasi-static | transient heat conduction with Robin boundary conditions in 2D and 3D |
| [H1/QuasiStatic/Taser](H1/QuasiStatic/Taser) | quasi-static | electro-quasistatic voltage pulse in tissue, `div((kappa + epsilon d/dt) grad phi) = 0`: `Taser2D.m` (forearm cross section), `TaserTestGeo.m` (3D leg with ring electrodes), `simpleTest.m` (two layers, exports E and J), `Simple2D.m` (purely capacitive), `Taser.m` (3D leg); `eqsTimeStepping.m` is the shared implicit Euler scheme |

## H(curl) (`HCurlElement`)

| Example | Time | Description |
|---|---|---|
| [HCurl/QuasiStatic/FerriteCoil](HCurl/QuasiStatic/FerriteCoil) | quasi-static | stationary current in a coil (H1), then the time-harmonic eddy current problem `S + i*omega*M` for the vector potential |
| [HCurl/QuasiStatic/Coil](HCurl/QuasiStatic/Coil) | quasi-static | the same workflow for an air coil — **incomplete**: the eddy current part uses `v_potential`/`dofsNe`, which are not set up in the script |
| [HCurl/QuasiStatic/AirCoil](HCurl/QuasiStatic/AirCoil) | quasi-static | geometry only (`aircoil.geo`, `.brep`) |
| [HCurl/FullWave/MaxwellEigenvalue2D](HCurl/FullWave/MaxwellEigenvalue2D) | full wave | eigenvalues of a 2D PEC resonator, convergence study (`runConvergence.m`) |
| [HCurl/FullWave/RectangularCavity](HCurl/FullWave/RectangularCavity) | full wave | eigenvalues of a 3D cavity resonator, live script `rect.mlx` |
| [HCurl/FullWave/RectangularCavityConvergence](HCurl/FullWave/RectangularCavityConvergence) | full wave | convergence of the 3D cavity eigenvalues in mesh size and order (`runConvergence.m`) |

## H(div) (`HDivElement`)

| Example | Time | Description |
|---|---|---|
| [HDiv/Static/CylinderFlow](HDiv/Static/CylinderFlow) | static | potential (Darcy) flow past a cylinder in a channel, mixed formulation with an exactly divergence-free velocity |
| [HDiv/QuasiStatic/CylinderTransport](HDiv/QuasiStatic/CylinderTransport) | quasi-static | pulsating flow past an off-center cylinder with a conservative tracer transport and UCD export for ParaView |

## Tools

| Example | Description |
|---|---|
| [Tools/PointLocation](Tools/PointLocation) | locating points in 2D and 3D meshes |

