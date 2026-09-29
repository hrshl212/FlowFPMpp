# Provenance

All source code in this repository was written by the author for this
project and is released under the BSD 3-Clause License (see `LICENSE`).

- The stair-step immersed-body solver (`src/fpm_stairstep_mpi_3d.f90`),
  tests, benchmark, demo, and tools are original. The blocked-cell (stair-step)
  treatment of the body is a standard textbook discretization.
- The serial 3-D Neumann Poisson and FPM-term kernels (`src/fpm_laplace_3d.f90`,
  `src/fpm_terms_3d.f90`), the halo exchange (`src/fpm_halo_mpi_3d.f90`), and
  the NVTX wrappers come from the author's earlier GPU port of an FPM
  postprocessor. They implement general numerics and the published FPM
  volume-term formulas.
- No code from the VICAR3D immersed-boundary solver is included, and neither
  are reimplementations of its ghost-cell/image-point interpolation, its
  immersed-boundary operators, or readers for its file formats. That work
  exists as a separate private code base.
- Tests, the benchmark, and the cylinder demo generate their own inputs.
- `examples/data/foil_frame_18400.bin` is one cropped frame (grid, u, v,
  omega_z, and solid flag) from the author's 2-D flapping-foil simulation
  computed with the VICAR3D code, published by the author as example data. It
  was converted from the simulation output by a private tool; no VICAR3D file
  reader is included. The figure `docs/images/vicar3d_foil_frame.png` shows the
  same frame as evaluated by the private postprocessor.
- The force-partitioning method follows the literature cited in `CITATION.cff`
  and the README.
