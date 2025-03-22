## NERDSS

Structure-Resolved Reaction-Diffusion Simulation Software by Johnson Lab, JHU

### [Website](https://mjohn218.github.io/NERDSS/)

### Version Information

- Current Version: 1.2.1
    - This version includes the latest stable features

- Parallel NERDSS
    - A parallelized version of NERDSS is in the [mpi](https://github.com/mjohn218/NERDSS/tree/mpi) branch.

#### Run a quick trial with our server

Go to the [NERDSS server](http://18.188.233.206:5000/).

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
3. Refer to the [ionerdss repository](https://github.com/JohnsonBiophysicsLab/ionerdss) for more details.
