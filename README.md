# FlowFPM++

Force-partitioning influence fields and vortex-induced forces on CPUs and
NVIDIA GPUs, in Fortran 2018 with MPI and OpenACC.

![Vorticity and vortex-induced force density around a stair-step cylinder](docs/images/demo_cylinder_vortices.png)

*Demo (`examples/`): a stair-step cylinder in a synthetic street of
Lamb-Oseen vortices. Top: vorticity. Middle and bottom: vortex-induced force
density for drag and lift, with the influence potentials phi_x and phi_y as
contour lines (dashed where negative). The flow is a prescribed kinematic
input, not a Navier-Stokes solution.*

## What it computes

The force-partitioning method (FPM) splits the pressure force on a body into
contributions from vorticity, added mass, viscous effects, and so on. It does
this by projecting the momentum equation onto one auxiliary *influence
potential* per force direction. For direction `i`, `phi_i` solves

```
laplacian(phi_i) = 0          in the fluid
n . grad(phi_i)  = n_i        on the body surface
n . grad(phi_i)  = 0          on the outer boundary
```

The vortex-induced force is then a volume integral over the fluid:

```
F_i(vortex) = -rho * integral_fluid (omega x u) . grad(phi_i) dV
```

FlowFPM++ solves for `phi_i` and evaluates the local densities
`-(omega x u) . grad(phi_i)` and `-a . grad(phi_i)` (for a supplied field `a`)
together with their integrals (per unit density). It does not include the
surface terms of the full FPM decomposition.

### Discretization

- **Grid.** Cell-centred, second order on stretched Cartesian grids.
- **Body: stair-step (blocked cells).** Cells whose centres lie inside the
  body are solid. They are excluded from the solve and hold `phi = 0`. Every
  face between a fluid and a solid cell carries the body condition. The
  condition `n . grad(phi_i) = n_i` does not change when `n` flips sign, so on
  a grid-aligned face it reduces to prescribing the face gradient along that
  face's axis: 1 on faces normal to axis `i`, 0 on the others. Near the body
  the geometry is only first-order accurate (see Validation).
- **Compatibility.** A pure Neumann problem needs the source and the boundary
  flux to balance. One uniform correction to the outer-boundary data enforces
  this (it is zero to roundoff for a closed body).
- **Gauge.** The fluid-volume mean of `phi` is removed.
- **Solver.** Weighted Jacobi,
  `phi <- 0.2 phi + 0.8 (sum(a_nb phi_nb) - f) / sum(a_nb)`. The damping
  removes the checkerboard two-cycle of plain Jacobi on this operator.
- **Force terms.** Cell gradients are the average of the two face gradients
  on each axis, using the prescribed values on outer and body faces.

## Design

- **One code path.** Every solver routine takes a `use_openacc` switch. The
  CPU reference is the same source with OpenACC disabled, so CPU/GPU
  comparisons use the same algorithm and memory layout.
- **Mask-free sweep.** Prescribed outer and body face gradients are folded
  once into a source array, and solid and outer neighbours are dropped from a
  precomputed diagonal. Solid cells stay at zero, so the Jacobi kernel reads
  no masks and has no branches.
- **Device-resident frame.** The operator arrays are built on the device, and
  the sweeps ping-pong between two buffers. Unknowns already on the device are
  not copied again, so a caller's data region can keep the flow fields and
  `phi` resident from the solve through the force terms (the benchmark does
  this). Per sweep, only the one-layer x-y halo faces cross the host-device
  boundary, and noncontiguous x faces are packed on the device into one
  transfer per face.
- **MPI.** 2-D (x-y) Cartesian decomposition with host-staged halo messages;
  z is not decomposed.
- **Double precision** throughout.
- **Instrumented.** NVTX ranges mark setup, sweeps, gauge, and the force-term
  kernel. They compile to no-ops in CPU builds.

## Validation

Every test generates its own inputs; none needs external data.

**Exact discrete solution** (`tests/test_stairstep_exact_3d.f90`). The grid is
stretched, 24x20x10, with four blocked bodies, one of them touching the outer
boundary. With a constant face gradient `g` on the body and outer faces, the
field `phi = g . x` solves the discrete problem exactly, and its gradient and
force densities have closed forms. The test checks this on 1, 4, and 6 ranks
on the CPU, and on 1 and 4 ranks on the GPU:

| Quantity | Max error |
|---|---|
| phi (converged to residual 1e-12) | 1.0e-13 |
| Gradient | 3.4e-13 |
| Vortex density | 2.9e-13 |
| Integrated terms | 3.3e-15 |
| Distributed vs one-rank field (600 fixed sweeps) | 3.6e-16 |

**Sphere convergence** (`tests/test_stairstep_sphere_3d.f90`). The x
influence potential of a sphere of radius `a` is the dipole
`phi = -a^3 x / (2 r^3)`. The test imposes its gradient on the outer faces of
the box `[-3a, 3a]^3` and compares the fluid field with it. The added-mass
coefficient `C_m = -(1/V) * integral(phi n_x dA)` has the exact value 1/2.

| N^3 | Sweeps | RMS error | Far-field (r > 2a) max error | C_m |
|---|---|---|---|---|
| 12 | 2,100 | 1.08e-2 | 2.24e-2 | 0.897 |
| 24 | 8,800 | 9.27e-3 | 1.87e-2 | 0.770 |
| 48 | 35,600 | 4.56e-3 | 9.21e-3 | 0.639 |
| 96 | 143,000 | 2.14e-3 | 4.38e-3 | 0.566 |

The observed order approaches 1 for all three quantities (1.09, 1.07, and
1.06 on the finest refinement). That is the expected accuracy of a
stair-step body: usable, but a sharp-interface method would be needed for
precise surface quantities.

**CPU/GPU parity.** On every benchmark grid the GPU field matches the CPU
field to a maximum difference of 2.8e-17 (relative L2 at most 1.4e-15).
Results are reproducible only to roundoff, because the order of GPU
reductions varies from run to run.

**Toolchains.** The CPU suite passes with GNU Fortran 13 and Open MPI 4.1,
and with Intel `ifx` and Intel MPI; CI runs the GNU build.

## Performance

Hardware: one NVIDIA A100-PCIE-40GB; Intel Xeon Gold 6248R (2 x 24 cores),
of which 1 or 16 ranks are bound to cores. Software: NVIDIA HPC SDK 25.5; one
binary built with `-O2 -acc -gpu=cc80` serves both the CPU and GPU runs.
Workload: `benchmark_stairstep_3d`, the x influence field of a stair-step
sphere on a stretched grid followed by the FPM terms. Times are medians of
three runs after a warm-up.

**Per-sweep cost** is what a converged solve spends its time on, since those
take 10^4 to 10^6 sweeps. It is measured as the difference between the
2,500-sweep and 250-sweep solve times, divided by 2,250.

| Grid | Cells | 1 core | 16 ranks | A100 | A100 vs 16 ranks | A100 vs 1 core |
|---|---|---|---|---|---|---|
| 64x64x16 | 65,536 | 88 us | 23 us | 20 us | 1.1x | 4.4x |
| 128x128x32 | 524,288 | 721 us | 84 us | 26 us | 3.2x | 28x |
| 256x256x64 | 4,194,304 | 13.4 ms | 926 us | 115 us | 8.0x | 117x |
| 512x512x64 | 16,777,216 | - | 5.09 ms | 426 us | 12.0x | - |

On the two largest grids the A100 sustains 36-39 billion cell updates per
second. That is about 1.2 TB/s of memory traffic at 32 bytes per cell, or
75-81% of the card's 1.555 TB/s peak, so the sweep is bandwidth-bound.

**A complete short frame** (compatibility, 250 sweeps, gauge, and force
terms, including every host-device transfer):

| Grid | 1 core | 16 ranks | A100 |
|---|---|---|---|
| 64x64x16 | 0.0267 s | 0.0063 s | 0.0070 s |
| 128x128x32 | 0.218 s | 0.024 s | 0.028 s |
| 256x256x64 | 3.66 s | 0.261 s | 0.186 s |
| 512x512x64 | 15.2 s | 1.36 s | 0.704 s |

For short solves the GPU is limited by pageable host-device transfers of the
input fields, and 16 CPU ranks win on the small grids, which fit in the CPU
caches. Raw data: [`benchmarks/stairstep_benchmark.csv`](benchmarks/stairstep_benchmark.csv).

### Optimization log

Each change was profiled first (NVHPC `NVCOMPILER_ACC_TIME`), then measured
against the previous binary on the same node, in two interleaved passes. The
fields agreed to roundoff.

| Profile finding | Change | Measured effect (one A100) |
|---|---|---|
| The Jacobi kernel already ran at about 1.1 TB/s; a separate copy-back kernel cost another 64 us per sweep at 256x256x64 | Ping-pong between two padded buffers | Per sweep 182 us to 118 us (256x256x64) and 641 us to 420 us (512x512x64) |
| Solver setup ran on the host, and its source, diagonal, and inverse-diagonal arrays were uploaded on every solve | Build them on the device from the solid flags and face data | Included in the frame times below |
| The force-term call uploaded 11 fields (about 77 ms) and downloaded 5 (about 43 ms) around a 0.7 ms kernel | Frame-level data region: fields and phi stay resident between the solve and the terms | Terms 0.113 s to 0.022 s (256x256x64) and 0.442 s to 0.089 s (512x512x64) |
| Combined | | 250-sweep frame 0.300 s to 0.180 s (256x256x64) and 1.158 s to 0.708 s (512x512x64) |

## Limitations

- **Stair-step geometry.** The body is represented to first order near its
  surface (see the sphere study); surface-sensitive quantities converge slowly.
- **Jacobi.** The solver needs O(N^2) sweeps to converge; the demo takes
  about two million sweeps. Multigrid or Krylov acceleration is the obvious
  next step and would change the performance picture.
- **Scaling not measured beyond one GPU.** MPI messages are host-staged (no
  GPU-aware MPI), and multi-GPU scaling has not been measured.
- **Volume terms only.** The code computes the vortex-induced and advection
  volume terms, not the surface terms of the full FPM decomposition.
- **Synthetic demo.** The demo flow is kinematic, not a flow solution.

## Background

FlowFPM++ is the public, rewritten core of the author's GPU port of an FPM
postprocessor for VICAR3D immersed-boundary simulations. That private code
base is built around the VICAR3D code and is not included here. On a production 3-D frame of about 87 million cells, its
OpenACC solve agreed with the CPU reference to 5.7e-17 and ran in 0.91 s on
one A100, against 12.2 s for a single CPU rank. The body treatment, solver
core, benchmark, and tests in this repository are new.

![Vorticity and x-direction vortex-induced force density around a flapping foil](docs/images/vicar3d_foil_frame.png)

*A production case: one frame (timestep 18400) of a 2-D flapping-foil
simulation computed with the VICAR3D code. Top: vorticity. Bottom: the
x-direction vortex-induced force density, evaluated from the saved fields by
the private postprocessor with the same FPM volume-term formula as here. Its
integral matches the stored VICAR3D result to 2e-15. The simulation data are
not distributed.*

## Build and run

CPU tests use GNU Fortran and Open MPI (Intel `ifx`/`mpiifort` also work):

```sh
make test FC=gfortran
make mpi-tests FC=gfortran MPIFC=mpifort MPIEXEC="mpirun --oversubscribe"
```

GPU build (NVIDIA HPC SDK with its bundled MPI):

```sh
make acc ACCGPU=cc80                     # cc80 = A100
mpirun -np 1 ./build_acc/test_stairstep_exact_3d openacc
mpirun -np 1 ./build_acc/test_stairstep_sphere_3d openacc extended
mpirun -np 1 ./build_acc/benchmark_stairstep_3d benchmark 256 256 64 250
mpirun -np 1 ./build_acc/demo_cylinder_vortices gpu 360 240 demo.bin
python examples/plot_demo.py demo.bin demo.png      # needs numpy + matplotlib
```

Benchmark usage is `benchmark_stairstep_3d [benchmark|cpu|gpu] [nx ny nz [sweeps]]`.
`cpu` mode runs on any rank count whose 2-D layout divides the grid. Example
Slurm scripts with placeholders are in [`cluster/`](cluster/).

## Repository layout

| Path | Contents |
|---|---|
| `src/` | Stair-step solver and force terms, halo exchange, serial 3-D kernels, NVTX wrappers |
| `tests/` | Self-contained tests and the benchmark |
| `examples/` | Cylinder/vortex-street demo and its plotting script |
| `benchmarks/` | Measured benchmark results |
| `docs/` | Provenance statement and figures |
| `cluster/` | Example Slurm scripts |
| `tools/` | Pre-publication scan |

## License and citation

BSD 3-Clause (see [`LICENSE`](LICENSE)); provenance is described in
[`docs/PROVENANCE.md`](docs/PROVENANCE.md). If you use this code, please cite
it (see [`CITATION.cff`](CITATION.cff)) and the works it builds on:

- H. S. Raut, J.-H. Seo, and R. Mittal, "Hydrodynamic performance and scaling
  laws for a modelled wave-induced flapping-foil propulsor," *Journal of
  Fluid Mechanics* 999, A1 (2024). https://doi.org/10.1017/jfm.2024.530
- R. Mittal, H. Dong, M. Bozkurttas, F. M. Najjar, A. Vargas, and
  A. von Loebbecke, "A versatile sharp interface immersed boundary method for
  incompressible flows with complex boundaries," *Journal of Computational
  Physics* 227, 4825-4852 (2008). https://doi.org/10.1016/j.jcp.2008.01.028
