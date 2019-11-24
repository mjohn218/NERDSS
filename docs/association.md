\page association Association Algorithm
\brief A description of the algorithm used when a binding event between two molecules occurs.
\tableofcontents

- See \ref Associate for relevant classes

\section algo Algorithm

  - See <a href="angle_documentation.pdf"><b>here</b></a> for a more complete (math-ey) description of the rotations.
  - Five angles are used to constrain the molecules into their bound conformation (see below)
  - Note that all translations/rotations are done proportional to each complex's total diffusion constants, i.e. a complex which a diffusion constant twice that of the other will move twice as far.
  
  1. The two reacting interfaces are placed at contact by translation along the vector between them.
  2. The complexes are pushed apart to sigma by translation along the interface-interface vector
  3. \f$\theta\f$ rotations are performed, by determining the difference between the current and desired \f$\theta\f$ angles and determining a rotation quaternion based on that difference and the rotation axis vector, the cross product between the (center of mass)-(reacting interface) and sigma vectors. \f$\theta_2\f$ is corrected first, with \f$\theta_1\f$ following.
  4. Next, the omega angle is corrected, using sigma as the axis of rotation.
      1. Determine the current omega angle:
          - The two (center of mass)-(reacting interface) vectors are projected onto a plane looking down the rotation axis (an orthographic projection)
          - The sign of the angle between them is determined by taking the cross product and looking at the sign of the z value (either pointing towards the camera or away from it)
      2. Determine the two rotation quaternions for the reactants (positive and negative rotations, based on their contributions to the total diffusion constant)
      3. Rotate and translate.
  5. Phi rotations are performed, using the (center of mass)-(reacting interface) as the axes of rotation.

  - All rotations are performed using the reacting interface as the origin, e.g. for theta 2, reactant 2 is the dominant interface, and rotations are done using the vector between it and whatever interface or center of mass is being rotated.
  - For all three rotation methods, if the first rotation fails to produce the desired angle, the rotation is reversed and the centers of mass and interfaces are returned to their original positions. The positive and negative half angles are then swapped and the rotation is performed with those new angles (i.e. the rotation is reversed in direction from the original). If this still produces the wrong angle, the association is cancelled. 

\image html docs/associate_figs/algo_flowchart_horiz.pdf

---

\section angles Angles

  - Five angles are used to orient the two associating molecules, \f$ \theta_1,~\theta_2,~\phi_1,~\phi_2,~\omega \f$
  - Using a toy model of two proteins interaction (red and blue) through an interface (iface1 and iface2), the below show the angles in practice.

\image html docs/associate_figs/angles/assoc_before.pdf

\subsection theta Theta

  - Angle between the (center of mass)-(reacting interface) and the \f$ \sigma \f$ vectors. Each association reaction has two \f$ \theta \f$ angles.

\image html docs/associate_figs/angles/theta.pdf

\subsection phi Phi

  - Dihedral between an arbitrary (by default the z-axis, or x- or y-axis if an interface to center of mass vector lies along it), center of mass to interface, and sigma vectors for each protein. Each association reaction has two \f$ \phi \f$ angles.

\image html docs/associate_figs/angles/phi.pdf

\subsection omega Omega

  - Dihedral between protein 1's center of mass to interface vector, the \f$ \sigma \f$ vector, and protein 2's center of mass to interface vector. Each association reaction has one \f$ \omega \f$ angle.

\image html docs/associate_figs/angles/omega.pdf

---

\section examples Examples

\subsection homo Homogenous Association

  - AP2 associationg with a clathrin-AP2 complex
\image html docs/associate_figs/ap2_example.pdf

\subsection hetero Heterogenous Association

  - Clathrins associating with each other
\image html docs/associate_figs/pucker_clat_example.pdf

---

\section Math

\subsection Quaternions
The algorithm uses quaternions for all rotations

\subsubsection intro Introduction
  > Quaternions came from Hamilton after his really good work had been done; and, though beautifully ingenious, have been an unmixed evil to those who have touched them in any way...
  > - Lord Kelvin, 1892

  - Number system which extends complex numbers
  - Equivalent to the quotient of two vectors, expressed in the form \f$ a + b\mathbf{i} + c\mathbf{j} + d\mathbf{k}\f$, where \f$ a,~b,~c,\f$ and \f$d\f$ are real numbers and \f$\mathbf{i},~\mathbf{j},\f$ and \f$\mathbf{k}\f$ are fundamental quaternion units.
    - Or, in the language of complex numbers, \f$b\mathbf{i} + c\mathbf{j} + d\mathbf{k}\f$ is the complex part, with \f$a\f$ being the real part of the quaternion \f$Q\f$
  - Relevant to rotations, quaternion multiplication is non-commutative.

\subsubsection why Why Quaternions?
  - Avoids gimbal lock, which can arise when using Euler angles (which may or may not be relevant in this algorithm, I'm not sure)
  - Faster (probably) and easier to read/more compact than matrices (definitely).

\subsubsection rotation Rotations with Quaternions
  - Axis-angle rotations, i.e. rotations through some angle \f$\theta\f$ using a rotation axis \f$\vec{u} = u_x\mathbf{i} + u_y\mathbf{j} + u_z\mathbf{k}\f$, can be represented by a quaternion through Euler's formula \f$\left( e^{ix} = \cos x + i\sin x \right)\f$:
    \f[
      Q = e^{\frac{\theta}{2}{(u_x \mathbf{i} + u_y \mathbf{j} + u_z \mathbf{k})}} = \cos\frac{\theta}{2} + (u_x \mathbf{i} + u_y \mathbf{j} + u_z \mathbf{k})\sin\frac{\theta}{2}
    \f]
  - Must use half angles for association, e.g. \f$\frac{\theta}{2}\f$ instead of \f$\theta\f$.
  - \f$Q (vec) Q^{-1}\f$ --> rotation of vector using quaternion \f$Q\f$. \f$Q^{-1}\f$ is the inverse of \f$Q\f$,
    \f$Q = (C, w)\f$ where C is the rotation angle and w is the axis of rotation \f$(a\mathbf{i} + b\mathbf{j} + c\mathbf{k})\f$.
  - Q must be normalized for rotation. Normalizes entire thing, including C
