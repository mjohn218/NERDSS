# Simulation of the Remodeling of HIV-1 Immature Lattice
## Constructing the initial lattices on the membrane
1. Forming a single spherical lattice by assembling the structure fully in solution using the inputs files in the folder named initialAssemblyInSolution.
2. Determining the time point that the lattice having the expected coverage and taking the restart file at that time point as the input file for the next step. I used the script in the folder named findSystemWithSpecifedMoleculeCopies.
3. Using the restart files from step 2 and put the assembled single lattice into a spherical system by linking the structure to the membrane using one lipid attachment per monomer. Meantime,modifying the reaction rates and other system parameters to the desired values. I used the cpp code in the folder named hackRestartFile.
4. Running the remodeling simulations using the restart file from step 3 in NERDSS.