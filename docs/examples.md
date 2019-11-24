\page examplesPage Example Simulations
\brief Example simulations, with descriptions, parameters, and reaction networks.
\tableofcontents

\section zeroCreate Zeroth Order Creation

\subsection zeroCreateDesc Description

Demonstration of FPR reproducing theoretical zeroth order creation.

\subsection zeroCreatePros Molecules

| Name | FPR Name | Partners |
|------|----------|----------|
| A | A | None

\subsection zeroCreateRxns Reactions

Null -> A(a) 1E-5 Ms<sup>-1</sup>

\subsection zeroCreateResults Results

<!-- <img src="docs/examples/zero_create/zeroth_order_creation_indiv.png" width="300" height="300"/> -->
![Zeroth order creation, individual simulations](docs/examples/zero_create/zeroth_order_creation_indiv.png =500x) ![Zeroth order creation, average](ocs/examples/zero_create/zeroth_order_creation_avg.png =500px)

Plot of 100 individual trajectories (left, red) and the average of those trajectories (right, red) against the theoretical result (both, dashed black).

---

\section flatClat Clathrin Sheet Formation

\subsection flatClatIrrev Irreversible Ring Closure
\subsubsection flatClatIrrevDesc Description

Formation of clathrin sheets in solution (3D), with irreversible ring closure. Binding of reacting interfaces within sigma has probability unity, so only two dissociation events in one or two timesteps can result in a ring breaking apart.

\subsubsection flatClatIrrevPros Proteins
| Name | FPR Name | Partners |
|------|----------|----------|
| clathrin | clat(cd1, cd2, cd3) | clat(cd1, cd2, cd3) |

\subsubsection flatClatIrrevParams Simulation Parameters

Download the input and molecule files here.

\subsubsection flatClatIrrevRxns Reactions

clat(cd1) + clat(cd1) <-> clat(cd1!1).clat(cd1!1) onRate, offRate

clat(cd2) + clat(cd2) <-> clat(cd2!1).clat(cd2!1) onRate, offRate

clat(cd3) + clat(cd3) <-> clat(cd3!1).clat(cd3!1) onRate, offRate

clat(cd2) + clat(cd2) <-> clat(cd2!1).clat(cd2!1) onRate, offRate

clat(cd1) + clat(cd2) <-> clat(cd1!1).clat(cd2!1) onRate, offRate

clat(cd1) + clat(cd3) <-> clat(cd1!1).clat(cd3!1) onRate, offRate

| K<sub>D</sub> | on rate (M<sup>-1</sup>s<sup>-1</sup>) | off rate (s<sup>-1</sup> |
|--|-----|-----|
| 0.2 | 18.87 | 1.334 |
| 1 | | |
| 100 | | |

\subsubsection flatClatIrrevResults Results

Shown here is the number of clathrin-clathrin bonds (clathrin legs) over time, for each K<sub>D</sub> value.

<img src="docs/examples/irrev_flat_clat/species.png" width= 300px>

---

\section repress Genetic Oscillator Model

\subsection repressDesc Description

Oscillations in protein expression levels can exhibit sensitivity to copy number fluctuations that are only captured in stochastic solvers such as Gillespie or single-particle RD.
Based on a simple model of circadian oscillations ([Vilar et al, PNAS 2002](https://www.pnas.org/content/99/9/5988)), we have simulated the behavior of an activator protein A and repressor protein R that are produced from mRNA transcribed from a single copy of a gene (one for each protein).
Coupling of A and R expression is driven by positive feedback of the activator A, which binds to each gene’s promoters to enhance transcription.
Protein R also binds to A to effectively degrade it, and all proteins and mRNA are degraded.
If the degradation rate of protein R is slow, the oscillations will end in the deterministic model, but persist in the stochastic method. Input file can be downloaded here.

\subsection repressPros Proteins

| Name | FPR Name | Partners |
|------|----------|----------|
| Activator protein | A(a~X~Y) | R, PrmA, PrmR |
| Repressor protein | R(r) | A |
| mRNA<sub>A</sub> | RNA(a) | None |
| mRNA<sub>R</sub> | RNR(r) | None |
| Protein A | PrA(a) | A |
| Protein R | PrR(r) | A |

\subsection repressParams Simulation Parameters

- Timestep: 50 us
- Iterations: 4E6
- Simulation volume: 4.189 um<sup>3</sup>
  - Box dimensions: 1600 nm x 1600 nm x 1600 nm

\subsection repressRxns Reactions

A(a~X) + R(r) -> A(a~X!1).R(r!1) 2000 nm<sup>3</sup>/us

PrA(a) + A(x~X) <-> PrA(a!1).A(a~X) 1000 nm<sup>3</sup>/us, 50 s<sup>-1</sup>

PrR(a) + A(x~X) <-> PrR(a!1).A(a~X) 1000 nm<sup>3</sup>/us, 100 s<sup>-1</sup>

PrA(a) -> PrA(a) + RNA(a) 500 s<sup>-1</sup>

PrR(r) -> PrR(r) + RNR(r) 0.01 s<sup>-1</sup>

RNA(a) -> RNA(a) + A(a~X) 50 s<sup>-1</sup>

RNR(r) -> RNR(r) + R(r) 5 s<sup>-1</sup>

PrA(a!\*).A(a~X!\*) -> PrA(a!\*).A(a~X!\*) + RNA(a) 500 s<sup>-1</sup>

PrR(r!\*).R(r!\*) -> PrR(r!\*).R(r!\*) + RNR(r) 50 s<sup>-1</sup>

RNA(a) -> NULL 10 s<sup>-1</sup>

RNR(r) -> NULL 0.5 s<sup>-1</sup>

A(a~X) -> NULL 1 s<sup>-1</sup>

R(r) -> NULL 0.2 s<sup>-1</sup>

The following reactions take the place of A.R -> R

A(a~X!\*).R(r!\*) -> A(a~Y!\*).R(r!\*) 100000 s<sup>-1</sup> (force to happen)

A(a~Y!\*).R(r!\*) -> A(a~Y) + R(r) 1 s<sup>-1</sup> (This is the true rate of the A.R -> R reaction)

A(a~Y) -> NULL 100000 s<sup>-1</sup> (force to happen)
