# Parallel NERDSS Development

Structure-Resolved Reaction-Diffusion Simulation Software by Johnson Lab, JHU
Parallelization with MPI

In addition to the serial version, NERDSS can also be parallelized using MPI (Message Passing Interface).

### Version Information

- Current Version: 1.0.0

### Installation

To build NERDSS with MPI support, you need:

1. A C++ compiler with MPI support:
    - For macOS, install OpenMPI with Homebrew: *brew install open-mpi*
    - For Ubuntu, install OpenMPI through apt: *sudo apt install openmpi-bin libopenmpi-dev*
2. GNU Scientific Library (GSL) version 1.0 or higher:
    - For macOS, install via Homebrew: *brew install gsl*
    - For Ubuntu, install via apt: *sudo apt install libgsl-dev*
3. To compile using make:
    - Navigate to the main directory
    - Run *make mpi* (with profiling support: *make mpi ENABLE_PROFILING=1*)
    - The executable will be placed in the ./bin directory

### Running Simulations

Example input files are located in the subdirectories within the sample_inputs folder. They can also be generated using the python GUI program, which is included in the ioNERDSS tool.

To start the simulation, use the command *mpirun -np $nprocs  ./nerdss_mpi -f parms.inp*.

To debug the code, use the command *mpirun -np 2 xterm -e gdb --ex 'b error' --ex r --args ./nerdss_mpi -f parms.inp -s 123*.

### Limitations

Note that the size of the largest complex in the simulation cannot be larger than half of the size of one rank, due to the parallelization scheme used in NERDSS. If this limitation is exceeded, the simulation will produce incorrect results.

### Analyzing Results

1. Use the ioNERDSS PyPi library for visualizing simulation results.
2. Install ioNERDSS with *pip install ioNERDSS*.
3. Refer to the [ioNERDSS repository](https://github.com/mjohn218/io_nerdss) for more details.

### NERDSS Parallel Developer Guide
 [NERDSS_Parallel_Developer_Guide](./NERDSS_Parallel_Developer_Guide.pdf)