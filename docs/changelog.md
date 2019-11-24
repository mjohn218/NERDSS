\page changelog Change Log
\brief Pretty self-explanatory

\section todo ToDo

- change PBC vs NOPBC to cmdline flag
  - In fact, it currently doesn't even support PBC
- Currently doesn't read exponentials from input files properly, i.e. `1E8` is read in as `1`.
- make proteins in the same complex able to bind at other interfaces (`evaluate_binding_pair_com.cpp`)
- VMD does not allow for trajectories that have changing atom numbers
  - Need to switch from `.xyz` format trajectories (maybe to `.dcd`?)
- Read/Write restart files
- Allow empty interface lists in reaction block (e.g. `A()`)
  - Would allow the following two reactions to be replaced with one, as long as they had the same rate

        PrA(a!*).A(a~X!*) -> PrA(a!*).A(a~X!*) + RNA(a)
        PrA(a) -> PrA(a) + RNA(a)
        // replaced by
        PrA() -> PrA() + RNA(a)

- Idea: store probabilities/factorials for zeroth and first order reactions. Might be faster.
- Idea: Quick and dirty parallelization through multithreading (split adjacent cell interaction search into two threads)?
- Idea: Move species copy numbers from `.mol` files to main `.inp` file, to more closely follow BNGL, e.g.
  
      start molecules
          PIP2 1E5
          AP2 200
      end molecules

- Change `start molecules` block to `start species` block, to more closely follow BNGL

---

\section change Change Log

Warning: There are gaps in this

\subsection feb19 February 2019

\subsubsection feb2819 02/28/2019

- Fixed various bugs parsing creation and state change reactions involving bound species
- Changed reactant checking for unimolecular reactions to check if bound state is correct by checking all interfaces of target Molecule
- Added rudimentary plotting tool for observables file (in `new_version/sample_inputs/analysis_scripts`)

\subsubsection feb122019 02/12/2019

- Added ability to perform unimolecular destruction reaction
- Added ability to perform  unimolecular creation reaction
- Added overlap check when performing creation
- Added restart file writing (working on reading functionality)

\subsubsection feb072019 02/07/2019

- Fixed issue where parser didn't read correct rates for zeroth order creation reactions.
- Added ability to perform zeroth order creation reactions (need to write overlap test)
- Added ability to perform unimolecular state changes
- Added unfinished restart file functionality and unfinished species tracking.
- Updated docs.

\subsubsection feb0519 02/05/2019

- Compiles and runs with GCC 4.9.0
  - Removed brace initialization for object references and arrays/vectors (see \ref guidelines).
  - Replaced `abs` with `std::abs` (this was giving horribly strange errors on GCC 4.9.0)
- Fixed an issue where I broke association due to deleting Interface::State::rxnPartners list population.
- Updated documents with dependencies and development guidelines.
- Added check for if the molecules are outside of the box during initial coordinate creation.
- Fixed issue where the initial coordinate xyz file did not have the correct number of atoms in the header.

\subsubsection feb0119 02/01/2019

- Added `<algorithm>` includes to `shared_reaction_functions.cpp`, `parser_functions.cpp`, and `dissociation.cpp` to fix issue with not finding `std::find_if` on Linux machines.


\subsection jan19 January 2019

\subsubsection jan3119 01/31/2019

- Updated documents with a quickstart guide

\subsubsection jan2919 01/29/2019

- Merged association/BNGL parser development branch with master branch.

\subsubsection jan2819 01/28/2019

- Added a global constants header, `globals.hpp`
- Fixed unimolecular reaction parsing

\subsubsection jan2519 01/25/2019

- Added bimolecular state change reaction parsing.
- Updated documentation with build instructions.

\subsubsection jan2419 01/24/2019

- Fixed bug in `shared_reaction_functions.cpp`, `evaluate_binding_pair()`, where it was looking at reactions for all states of an interface on proIndex1, instead of only its current state.

\subsubsection jan2319 01/23/2019

- Parser now properly reads in bimolecular state change reactions (e.g. kinase facilitated phosphorylation, K(k) + A(a~X) <-> K(k) + A(a~Y)). See `parser_functions.hpp` for information on limitations.
- Updated documentation in `parser_functions.hpp`, `parser_functions.cpp`, `class_bngl_parser.hpp`, and `class_bngl_parser_functions.cpp`.

\subsubsection jan2219 01/22/2019

- Changed `ParsedMol::rxnReactants` and `ParsedMol::rxnProducts` from `std::map<std::string, ParsedMol::IfaceInfo>` to `std::vector<ParsedMol::IfaceInfo>`, because there was no reason for it to be a map.

\subsubsection jan1419 01/14/2019

- Updated angle documentation. Now includes angle bounds and positive/negative angle definitions for torsion angles.

\subsubsection jan1119 01/11/2019

- Fixed bug in `angleSignIsCorrect()`, where it was errantly using an x-y projected vector instead of an x-z projected vector during 2D association.

\subsection dec18 December 2018

\subsubsection dec1918 12/19/2018

- New version works in 3D. Testing in 2D.

\subsection oct18 October 2018

\subsubsection oct3018 10/30/2018

- Fixed an issue with norms not being associated with the correct reactant. In BNGL reaction, reactant 1 is always the leftmost reactant, and should be associated with angles[0, 2] and norm1, with reactant two being rightmost and associated with angles[1, 3] and norm2.

\subsubsection oct0818 10/08/2018

- Added `<algorithm>` includes so it properly compiles on Linux machines

\subsection sept18 September 2018

\subsubsection sept1818 09/18/2018

- Data structure changes, combined parameter file and reaction file into one input file.

\subsection june18 June 2018

\subsubsection june1418 6/14/18

- function to choose which reaction to use now uses the correct rate state based on ancillary interfaces

\subsubsection june1318 6/13/18

- parser and association functions can now talk to each other
- fixed a bug in parser where it would incorrectly parse reactant wildcard bonds 
- function to choose which reaction to use still only chooses the first one, regardless of ancillary interfaces present on the reacting molecule (in the process of fixing)

\subsubsection june0818 6/8/18

- fixed bug in input parser: creation of a BackRxn would assign the on rate to its rateList, not the off rate

\subsubsection june0718 6/7/18

- updated docs
- moved RxnBase::stateChangeList to RxnBase::RateState
  - conditional rates now also applies to state changes of interfaces

\subsubsection june0618 6/6/18

- added Vector class to association and updated to reflect changes in class structure made while developing input parser
- fixed bug in calculate_omega() and calculate_phi() in association

\subsubsection june0418 6/4/18

- added docs

\subsection may18 May 2018

\subsubsection may3118 5/31/18

- updated parameter and added molecule file parsing
  - BNGL style parameter, molecule, and reaction blocks
- parsing of symmetric reactions now works

\subsubsection may2918 5/29/18

- creation and destruction reaction parsing now works
- added parsing of simulation parameters
- git now ignores test and cmake directories, .dat, .mol, .info, and .inp files

\subsubsection may2418 5/24/18

- reaction spliting based on implied states now works
  - if an interface lacks a state in a reation, but has named states in its MolTemplate, the reaction is split into multiple reactions, one for each state
  - seems to work with conditional rates

\subsubsection may2218 5/22/18

- added conditional rates
  - if reactants and products of two reactions are the same, but they have different ancillary interfaces, they are counted as the same reaction with different rates (as stored in std::vector<RateState> rateList)

\subsubsection may2018 5/20/18

- added detection of interface state changes
- calculation of angles from provided complex internal coordinates now works (mostly proof of concept for GUI)

\subsubsection may1518 5/15/18

- cleaned up directories
- committed in order to create feature_develop-bngl_parser branch

\subsubsection may0418 5/4/18

- adding cooperativity for states and wildcard bongs in reactants
  - not finished

\subsection jan18 January 2018

\subsubsection jan1418 1/4/18

- fairly certain NOPBC is free of bugs (multiple tests output the same as old version)
- starting on PBC version

\subsection dec17 December 2017

\subsubsection dec2217 12/22/17

- last known bug: sometimes has an error which manifests in pulling a lipid off the membrane in sweep_separation_complex_rot_memtest
- added command line flags:
    --rxnfile file  (reads the rxn.inp file)
    --parmfile file (reads the parms.inp file)
    --seed int(seed) (sets a custom seed, defaults to pseudorandom without)
    --extended-output 0,1 (writes out bases/complexes every configwrite timesteps, slows down ~x2)
    --debug_force_dissoc 0,1 (forces dissociation every time)
    --debug_force_assoc 0,1 (forces association every time)

\subsubsection dec1417 12/14/17

- still debugging

\subsubsection dec1317 12/13/17

- debugging still
- mostly looks good, testing different edge cases

\subsubsection dec1017 12/10/17

- produces identical output to old version, without association/dissociation
- debugging association

\subsubsection dec0718 12/7/17

- debugging association
- testing new Interaction class to take over bndlist, bndiface, bndpartner

\subsubsection dec0517 12/5/17

- fixed 12/4 bug
- debugging association

\subsubsection dec0417 12/4/17

- debugging issue where molecules with Dz = 0 will sometimes get traj[2] > 0
  - seems to only affect molecules which have ncross > 0

\subsubsection dec0117 12/1/17

- testing (in evaluate_binding_pair_com)

\subsection nov17 November 2017

\subsubsection nov3017 11/30/17

- full NOPBC compiles with no errors
  - in testing phase

\subsubsection nv2817 11/28/17

- up to checking for overlap in main

\subsubsection nov2117 11/21/17

- up to and including associate_freeleg in perform_association_sigma_com.cpp

\subsubsection nov1517 11/15/17

- evaluate_binding_pair_com.cpp and its children subroutines compile properly. still need to test

\subsubsection nov1417 11/13/17

- dissociate_sigma_com.cpp and its children subroutines compile properly. still need to test
- replacing c-style arrays to std::array where available

\subsubsection nov0817 11/8/17

- working on dissociate_sigma_com.cpp
- updated determine_which_complex.cpp

\subsubsection nov0617 11/6/17

- compiles and runs up to simulation start

\subsubsection nov0317 11/3/17

- working on get_bin2()

\subsection oct17 October 2017

\subsubsection oct3117 10/31/17

- compiles and runs up to and including gen_psf_system()
- renamed ind_com -> Complexlist

\subsubsection oct2717 10/27/17

- replaced ***coordcont with a Coord class within each InterfaceSite class object.

\subsubsection oct2617 10/26/17

- fixed the problem with read_inputs.cpp not handling Rin vector properly
