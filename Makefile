# FlowFPM++ build. CPU defaults use GNU Fortran and an MPI wrapper; override
# FC, MPIFC, FFLAGS, and MODFLAGS for other compilers, for example
#   make test FC=ifx MPIFC=mpiifort MODFLAGS='-module $(BUILD_DIR)'
# OpenACC targets use NVIDIA HPC SDK (nvfortran and its MPI wrapper).
ifeq ($(origin FC),default)
FC = gfortran
endif
MPIFC ?= mpifort
MPIEXEC ?= mpirun
FFLAGS ?= -O2 -std=f2018
BUILD_DIR ?= build
MODFLAGS ?= -J$(BUILD_DIR)
ACC_MPIFC ?= mpifort
ACCGPU ?= cc80
ACC_DIR ?= build_acc
ACCFLAGS ?= -O2 -acc -gpu=$(ACCGPU) -Minfo=accel

CORE := src/fpm_trace.f90 src/fpm_halo_mpi_3d.f90 src/fpm_stairstep_mpi_3d.f90 tests/fpm_test_grid.f90
CORE_ACC := src/fpm_trace_nvtx.f90 src/fpm_halo_mpi_3d.f90 src/fpm_stairstep_mpi_3d.f90 tests/fpm_test_grid.f90
MPI_PROGRAMS := test_stairstep_exact_3d test_stairstep_sphere_3d benchmark_stairstep_3d

.PHONY: all test mpi-tests acc foil clean

all: $(BUILD_DIR)/test_fpm_3d $(BUILD_DIR)/test_fpm_3d_nonuniform $(BUILD_DIR)/test_mpi_halo_3d \
	$(addprefix $(BUILD_DIR)/,$(MPI_PROGRAMS)) $(BUILD_DIR)/demo_cylinder_vortices $(BUILD_DIR)/foil_frame_fpm

$(BUILD_DIR)/test_fpm_3d: src/fpm_laplace_3d.f90 src/fpm_terms_3d.f90 tests/test_fpm_3d.f90
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $(MODFLAGS) -I$(BUILD_DIR) $^ -o $@

$(BUILD_DIR)/test_fpm_3d_nonuniform: src/fpm_laplace_3d.f90 src/fpm_terms_3d.f90 tests/test_fpm_3d_nonuniform.f90
	@mkdir -p $(BUILD_DIR)
	$(FC) $(FFLAGS) $(MODFLAGS) -I$(BUILD_DIR) $^ -o $@

$(BUILD_DIR)/test_mpi_halo_3d: src/fpm_halo_mpi_3d.f90 tests/test_mpi_halo_3d.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIFC) $(FFLAGS) $(MODFLAGS) -I$(BUILD_DIR) $^ -o $@

$(BUILD_DIR)/%: $(CORE) tests/%.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIFC) $(FFLAGS) $(MODFLAGS) -I$(BUILD_DIR) $^ -o $@

$(BUILD_DIR)/demo_cylinder_vortices: $(CORE) examples/demo_cylinder_vortices.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIFC) $(FFLAGS) $(MODFLAGS) -I$(BUILD_DIR) $^ -o $@

$(BUILD_DIR)/foil_frame_fpm: $(CORE) examples/foil_frame_fpm.f90
	@mkdir -p $(BUILD_DIR)
	$(MPIFC) $(FFLAGS) $(MODFLAGS) -I$(BUILD_DIR) $^ -o $@

# Real-frame example (examples/data/foil_frame_18400.bin): writes foil_frame_fpm.bin.
foil: $(BUILD_DIR)/foil_frame_fpm
	$(MPIEXEC) -np 1 $(BUILD_DIR)/foil_frame_fpm cpu

# Serial manufactured-solution tests.
test: $(BUILD_DIR)/test_fpm_3d $(BUILD_DIR)/test_fpm_3d_nonuniform
	$(BUILD_DIR)/test_fpm_3d
	$(BUILD_DIR)/test_fpm_3d_nonuniform

# MPI tests (host execution); about one minute on four cores.
mpi-tests: $(BUILD_DIR)/test_mpi_halo_3d $(BUILD_DIR)/test_stairstep_exact_3d $(BUILD_DIR)/test_stairstep_sphere_3d
	$(MPIEXEC) -np 4 $(BUILD_DIR)/test_mpi_halo_3d
	$(MPIEXEC) -np 1 $(BUILD_DIR)/test_stairstep_exact_3d
	$(MPIEXEC) -np 4 $(BUILD_DIR)/test_stairstep_exact_3d
	$(MPIEXEC) -np 4 $(BUILD_DIR)/test_stairstep_sphere_3d

# OpenACC builds; the same programs take an "openacc" or "gpu" argument.
acc: $(ACC_DIR)/test_mpi_halo_3d_acc $(addprefix $(ACC_DIR)/,$(MPI_PROGRAMS)) $(ACC_DIR)/demo_cylinder_vortices \
	$(ACC_DIR)/foil_frame_fpm

$(ACC_DIR)/test_mpi_halo_3d_acc: src/fpm_halo_mpi_3d.f90 tests/test_mpi_halo_3d_acc.f90
	@mkdir -p $(ACC_DIR)
	$(ACC_MPIFC) $(ACCFLAGS) -module $(ACC_DIR) -I$(ACC_DIR) $^ -o $@

$(ACC_DIR)/%: $(CORE_ACC) tests/%.f90
	@mkdir -p $(ACC_DIR)
	$(ACC_MPIFC) $(ACCFLAGS) -cudalib=nvtx -module $(ACC_DIR) -I$(ACC_DIR) $^ -o $@

$(ACC_DIR)/demo_cylinder_vortices: $(CORE_ACC) examples/demo_cylinder_vortices.f90
	@mkdir -p $(ACC_DIR)
	$(ACC_MPIFC) $(ACCFLAGS) -cudalib=nvtx -module $(ACC_DIR) -I$(ACC_DIR) $^ -o $@

$(ACC_DIR)/foil_frame_fpm: $(CORE_ACC) examples/foil_frame_fpm.f90
	@mkdir -p $(ACC_DIR)
	$(ACC_MPIFC) $(ACCFLAGS) -cudalib=nvtx -module $(ACC_DIR) -I$(ACC_DIR) $^ -o $@

clean:
	rm -rf $(BUILD_DIR) $(ACC_DIR)
