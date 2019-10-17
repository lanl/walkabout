# Amanzi - Walkabout

## particleDriver V1 – works with Amanzi versions before September 2017 

particleDriver is a C++ program for a particle advection system using Amanzi, LaGriT, and Walkabout. It automates the process of producing input and output files for LaGrit and Walkabout from an Amanzi checkpoint file and creates data for visualization that can be used in Paraview, PlumeCalc, and Meshlab. 

- This works ok for a box shaped (convex) mesh, but does not work for non-convex, non-Delaunay mesh.
- Any work around will result in missing tet elements, or inverted tet elements along the boundary.
- Walkabout Particles stop or have undefined behavior at these locations.


## particleDriver V2 – Amanzi(V.88+)

This new version particleDriver2 reads the h5 output from Amanzi Version 0.88 or newer. The model velocity fields are calculated on the input mesh vertices. It requires the tet5.inp mesh file used to write the Amanzi exodus file so both the Amanzi input mesh and output mesh are the same. 

- Amanzi has been modified to output Vxyz at original mesh vertices following the least squares framework from Painter, Gable, Kelkar. However, there can be ongoing issues with particles leaving the domain on boundaries. This behavior can be controlled by using the noflow boundary condition zone built into Walkabout.
- Uniform Material paraview and walkabout results nearly the same.
- Non-uniform Material paraview and walkabout differ for some particles.


## Workflow Using Amanzi plot mesh files

This method uses LaGriT to create a tet mesh from the hex plot mesh by connecting the plot mesh cell centers. The cell centers become the tet vertices. The original mesh boundary nodes are lost in this conversion. This method has the best results.



## New Workflow Development

particleDriver and workflow development uses Walkabout tests plus some new curved boundary models.

See test repo at https://github.com/amanzi/amanzi-walkabout

