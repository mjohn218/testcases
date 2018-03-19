# Irreversible Binding 2D

This test case examines the irreversible interaction between two, not necessarily identical, species. The 2D reversible binding test case assesses whether the simulation method depicts the proper kinetics and equilibrium for an irreversible reaction in two-dimensions.        

<pre>
Species and Diffusion Coefficients:
A       1 nm^2.us-1
B       1 nm^2.us-1
C       0.5 nm^2.us-1

Initial Conditions (Species, Number of Molecules):
A       1000 molecules
B       1000 molecules

Reaction (Equation, Kf, Kr):
A + B -> C   1 and 10 um^2.s-1    

MCell Parameters
Geometry:
Plane
Surface Area - 1 um^2
Side Length - 1 um

Time Step: 1e-6 s
Length of Simulation: 1 s
Number of Seeds (MCell): 10


VCell Parameters
Geometry:
Spherical Box
Cell Effective Radius - 0.080 um

Time Step: 1e-7 s
End Time: 1 s
Number of Trajectories Averaged: 10
</pre>
