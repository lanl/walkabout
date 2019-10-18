# Amanzi - Walkabout

## particleDriver V1 works with Amanzi versions before September 2017 

particleDriver is a C++ program for a particle advection system using Amanzi, LaGriT, and Walkabout. It automates the process of producing input and output files for LaGrit and Walkabout from an Amanzi checkpoint file and creates data for visualization that can be used in Paraview, PlumeCalc, and Meshlab. 

- This works ok for a box shaped (convex) mesh, but does not work for non-convex, non-Delaunay mesh.
- Any work around will result in missing tet elements, or inverted tet elements along the boundary.
- Walkabout Particles stop or have undefined behavior at these locations.


## particleDriver V2 Amanzi(V.88+) for non-convex meshes

This new version particleDriver2 reads the h5 output from Amanzi Version 0.88 or newer. The model velocity fields are calculated on the input mesh vertices. It requires the tet5.inp mesh file used to write the Amanzi exodus file so both the Amanzi input mesh and output mesh are the same. 

- Amanzi has been modified to output Vxyz at original mesh vertices following the least squares framework from Painter, Gable, Kelkar. However, there can be ongoing issues with particles leaving the domain on boundaries. This behavior can be controlled by using the noflow boundary condition zone built into Walkabout.
- Uniform Material paraview and walkabout results nearly the same.
- Non-uniform Material paraview and walkabout differ for some particles.
- Non-convex non-Delaunay mesh results in many negative coefficients that can impact results. Use a median .stor file as a work around, all volumes and coefficients are positive, but the geometry is not guaranteed to be orthogonal.


## plot_to_tet5 uses Amanzi hex mesh copied on to cell center median tet mesh 

This method uses LaGriT to create a tet mesh from the hex plot mesh by connecting the plot mesh cell centers. The cell centers become the tet vertices. The original mesh boundary nodes are lost in this conversion. 

- This method avoid interpolation of Amanzi values from cell center to cell vertices.
- This method has the best results.
- Non-convex non-Delaunay mesh results in many negative coefficients that can impact results. Use a median .stor file as a work around, all volumes and coefficients are positive, but the geometry is not guaranteed to be orthogonal.


See particleDriver and workflow development for Amanzi-Walkabout models:

   test repo at https://github.com/amanzi/amanzi-walkabout

