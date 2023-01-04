READ ME:

This model for spherical Gag self-assembly is obtained from the 5l93.pdb and 5l93_work.pdb(fetched version, including 18 monomers). There are two sets of input files: one for radius 50 nm and another for radius 65 nm. Gag monomer information is obtained in A_R50.mol/A_R65.mol, each Gag has one COM and 7 binding sites: b—dimer site; c — trim1; d-hex1; p — hex2; q — trim2; m — membrane binding site(outside the gag sphere); r — RNA binding site(inside the gag sphere). The last two binding sites and the COM are set to be on the same line. The variation among 18 monomers is corrected to avoid monomers assembling into spherical fragments with distinct curvatures.

There are three kinds of intro-Gag interactions: the homo-dimer interaction, the hexamer interaction mediated by the two distinct hexamer sites, and the trimer interaction mediated by the two distinct trimer sites. The distance between two consecutive hexagons is 7.7nm(R50)/8.1nm(R65). Interaction info (orientations and angels) is included in the parm_R50.inp/parm_R65.inp.

FullSphere_R50.pdb is the simulation result based on A_R50.mol and parm_R50.inp. It includes a complete sphere. When R=50, a complete Gag sphere includes ~3700 monomers.
