## NERDSS

Structure-Resolved Reaction-Diffusion Simulation Software by Johnson Lab, JHU

### Version Information

- Current Version: 1.2.0
    - This version includes the latest stable features
    - For release notes and changes, see the [Changelog](./Changelog.md)

- Parallel NERDSS
    - A parallelized version of NERDSS is in the [mpi](https://github.com/mjohn218/NERDSS/tree/mpi) branch.

### Installation

To build NERDSS, you need:

1. A C++ compiler:
    - For macOS, install XCode
    - For Ubuntu, install a compiler through apt
2. GNU Scientific Library (GSL) version 1.0 or higher:
    - For macOS, install via Homebrew
    - For Ubuntu, install via apt
3. To compile using make:
    - Navigate to the main directory
    - Run *make serial*
    - The executable will be placed in the *./bin* directory

### Running Simulations

1. Example input files are located in the subdirectories within the *sample_inputs* folder. They can also be generated using the [python GUI program](./gui.py), which is also included in the ioNERDSS tool.

2. To start the simulation, use the command *./nerdss -f parms.inp*.

### Analyzing Results

1. Use the ioNERDSS PyPi library for visualizing simulation results.
2. Install ioNERDSS with *pip install ioNERDSS*.
3. Refer to the [ioNERDSS repository](https://github.com/mjohn218/io_nerdss) for more details.
