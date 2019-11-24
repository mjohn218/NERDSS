\page input Input File Introduction
\brief Information on how to read, use, and format the input files for the software.
\tableofcontents

\section params Input File Format
  - Blocks are started with the `start` keyword and ended with the `end` keyword, e.g.

        start molecules
        end molecules

    - Three block types are used, `parameters`, `molecules`, and `reactions`. `molecules` must be defined before `reactions`
  - Molecule (MolTemplate) information is stored in `*.mol` files, with the prefix being the molecule with named identically to that used in the reactions

\subsection paramKeys Parameter Keywords
\subsubsection sys System Parameters
  | Keyword    | Description |
  |------------|-------------|
  |NumMolTypes | total number of molecules in the system (integer)|
  |NumIfaces   | total number of interfaces and states in the system (integer) |

\subsubsection simul Simulation Parameters
  | Keyword    | Description |
  |------------|-------------|
  | NItr<sup>✝</sup> | requested number of iterations (timesteps) (integer) |
  | deltat<sup>✝</sup> | timestep length (integer)|
  | timeWrite | interval to write timestep information (integer, default: 10) |
  | trajWrite | interval to write to coordinates file (integer, default: 10) |
  | restartWrite<sup>1</sup> | interval to write to restart file (integer, default: trajWrite) |
  | pdbWrite | interval to write individual pdb files (mainly for nonequilibrium simulations (integer, default: 0) |
  | fromRestart | Does the simulation start from a restart file. Not yet implemented (01/25/2019) (boolean, default: false) |
  | waterBox<sup>✝</sup> | water box dimensions (array, [x, y, z]) |

✝ indicates a required keyword

<sup>1</sup> Must be set such that restartWrite % trajWrite = 0

---

\subsection rxnKeys Reaction Block Keywords

| Keyword | Description |
|---------|-------------|
| onRate<sup>✝</sup>  | On rate of the reaction (double) |
| offRate | Off rate of the reaction, if reversible (double) |
| sigma | distance between the two reacting interfaces for a bimolecular reaction (double or array, [x, y, z], default: 1.0 nm) |
| assocAngles<sup>1</sup> | Angles for bimolecular association (array, [theta1, theta2, phi1, phi2, omega], radians) |
| norm1, norm2<sup>2</sup> | vector used to calculate the phi1 and phi2 angles (and sometimes omega) |
| observeLabel | label for tracking the reaction product (string, see Species Tracking) |
 
✝ indicates a required keyword for all reaction

<sup>1</sup> If the angle does not exist, e.g. two rods-type molecules associating with \f$\theta\f$'s equal to \f$\pi\f$ will not not have an omega angle (see \ref association), the angle should read "nan" or "null"  

<sup>2</sup> numbering must match the order of reactants (e.g. norm1 goes with the leftmost reactant  

- Reaction specific required keywords:
  - Bimolecular association: OnRate, Sigma, Angles
  - Reversible reactions: OnRate, OffRate

- Rate units:
  - Unbinding: s<sup>-1</sup>
  - Creation: M/s
  - Association:
    - 3D: nm<sup>3</sup>/us
    - 2D: nm<sup>2</sup>/us

---

\subsection paramExample Example File

    start parameters
        nItr = 5000
        deltaT = 0.1
    
        waterBox=[400,400,400]
    
        timeWrite = 1000
        trajWrite = 1000
        restartWrite = 1000
        fromRestart = false
    end parameters
    
    start molecules
        pip2
        ap2
        clat
    end molecules
    
    start reactions
        # Reactions
        ap2(b2clat, m2muh) + pip2(head) <-> ap2(b2clat,m2muh!1).pip2(head!1)
        onRate = 50 #On-rate unit: uM-1/s
        offRate = 10.0 #Off-rate unit: 1/s
        sigma = 1.0 #Sigma unit: nm
        assocAngles = [1.5708, 1.5708, nan, nan, M_PI] #Angles unit: radians
    end reactions

---

\section mol Molecule Files

\subsection molKeys Mol File Keywords

| Keyword   | Description |
|-----------|-------------|
| Name<sup>✝</sup>     | Molecule name (must match that used in the provided reactions (string)|
| Copies<sup>✝</sup>    | Copy numbers of the molecule (integer)|
| IsRod     | The molecule is one dimensional (boolean, optional, default: false) |
| IsLipid   | The molecule is a lipid (boolean, optional, default: false) |
| isPoint | The molecule is a point, with no radius (boolean, optional, default: false)
| D<sup>✝</sup>      | Translational diffusion constants (array, [x, y, z]) |
| Dr<sup>✝</sup>     | Rotational diffusion constants (array, [x, y, z]) |
| COM<sup>✝</sup>    | Center of mass coordinates (x,y,z). Starts internal coordinates block (see below) |
| state     | Define interface states (string, see below) |
| bonds  | Predefined bonds for creating the PSF (two strings, see below) |

✝ indicates a required keyword

  - COM denotes the start of the internal coordinates block, which follows this format
        
        COM       x y z
        ifaceName x y z
        ...

  - `ifaceName` must be the name of the iface, identical to that used in the reactions block
  - States are declared with the following format:
        
        state = $ifaceName~X~Y

    where `$ifaceName` is the name of the interface, as used in the reactions block, and `X` and `Y` are the state identities, which must be alphabetic characters
  - Bonds for the PSF are declared as a pair of interface names separated by whitespace, preceded by a declaration line:
    
        bonds = 1 # 
        $ifaceName1 $ifaceName2

---

\subsection molFileExample Example File

    ##
    # pip2 molecule information file
    ##
    
    Name = pip2 
    copies = 100
    isRod = true
    isLipid = true
    
    # translational diffusion constants
    D = [8, 1, 0]
    
    # rotational diffusion constants
    Dr = [1, 8, 0]
    
    # Coordinates
    COM    0.0000    0.0000    0.0000
    head   0.0000    0.0000    1.0000

    bonds
    COM head

---

\section rxn Reaction Declarations

\subsection rxnSyn Syntax

- Reactions are defined as such: `A(a) + B(b) <-> A(a!1).B(b!1)`, where `A` and `B` are the reacting molecules, `a` and `b` are the reacting interfaces, "`!`" denotes an interaction with index `1`, and "`.`" indicates the two molecules are interacting
- Reversible reactions are denoted by a double-headed arrow `<->`, as opposed to a single-headed arrow `->`,  and the existence of a defined off rate
- Interfaces must be uniquely named, at least on each molecule type
- States are allowed, and are not required to be binary. Denoted with a tilde, e.g. `A~P`
- Allowed one state change and one interaction change per reaction. Note that these are mutually exclusive w.r.t. individual interfaces, i.e. an interface can only change its state or interaction, not both
- Ancillary interfaces are allowed. These can include interfaces with/without states/interactions that do not change their state or interaction in the particular reaction, but are required for the reaction to occur.
  - Example:  
    Let's say a molecule `A` has two interfaces, `a1`, and `a2`, with `a2` having two states, `P` and `U`. If an interaction between `a1` and some interface `b` on molecule `B` is dependent on `a2` being in the `P` state, it could be written as:

        A(a1,a2~P) + B(b) <-> A(a1!1,a2~P).B(b!1)
  - Ancillary interactions are allowed in the reactants, but only as a wildcard, not with an index, e.g.
        
        A(a1!*,a2) ...
    denoting `a1` is bound to something, it doesn't matter what.

- Wildcard states are allowed in the reactants/products by omitting the state of an interface which has states
    - Example:  
      Using the above defined molecules `A` and `B`, this is allowed:
          
          A(a1,a2) + B(b) <-> A(a1!1,a2).B(b!1)

      Indicating the reaction can occur regardless of the state of `a2`


\subsection supportedRxns Supported Reaction Types

- Symmetric reactions
  - Example: `(A(a) + A(a) <-> A(a!1).A(a!1)`
- Bimolecular association
- Bimolecular state change
  - Enzyme facilitated state changes
- Unimolecular state change
- Dissociation
- Creation from concentration
  - Example: `0 -> A(a)`
- Creation from molecule
  - Example: `A(a) -> A(a) + B(b)`
- Destruction
  - Example: `A(a) -> 0`

\subsection limits Current Limitations
- A Molecule must have at least one explicit interface listed, e.g. `A()` is not allowed.
- Unimolecular destruction reactions can only have one associated rate and will only act upon a Molecule (e.g. it will only consider wildcards, not explicitly stated bound molecules)
- Bimolecular state changes have a few limitations:
  - Bimolecular state changes can only occur when both interfaces are free
  - The facilitator can only have one interface explicitly listed (this will be changed)
