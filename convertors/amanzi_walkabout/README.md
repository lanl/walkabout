# Amanzi - Walkabout


The goal is to add capability to Amanzi software to write Walkabout compatible files for modeling particle tracking. The main issue is that Amanzi is computing on general polyhedral elements and boundary faces using cell centers. Walkabout uses a tetrahedral mesh with properties defined on the cell vertices. 

LaGriT is used as part of this workflow to connect the Amanzi center points into tetrahedral elements and from this write files required by Walkabout. These include the the sparse matrix of geometric coefficients (FEHM .stor file), an adjacency face file, and other files as described in the [Walkabout Manual](https://lanl.github.io/walkabout/index.html).


## METHOD 1 
### particleDriver V1 works with Amanzi versions before September 2017 

particleDriver is a C++ program for a particle advection system using Amanzi, LaGriT, and Walkabout. It automates the process of producing input and output files for LaGrit and Walkabout from an Amanzi checkpoint file and creates data for visualization that can be used in Paraview, PlumeCalc, and Meshlab. 

- This works ok for a box shaped (convex) mesh, but does not work for non-convex, non-Delaunay mesh.
- Any work around will result in missing tet elements, or inverted tet elements along the boundary.
- Walkabout Particles stop or have undefined behavior at these locations.


## METHOD 2 
### particleDriver V2 Amanzi(V.88+) for non-convex meshes

This new version particleDriver2 reads the h5 output from Amanzi Version 0.88 or newer. The model velocity fields are calculated on the input mesh vertices. From Konstantin; Amanzi supports multiple conservative formulations (FV, monotone FV, mixed MFD). All of them provide normal component of the Darcy velocity on a mesh face.The nodal velocity is commuted using a constrained least square algorithm, where constraints are the Neumann boundary conditions.

*This requires the tet5.inp mesh file used to write the Amanzi exodus file so both the Amanzi input mesh and output mesh are the same.* 

- Amanzi has been modified to output Vxyz at original mesh vertices following the least squares framework from Painter, Gable, Kelkar. However, there can be ongoing issues with particles leaving the domain on boundaries. This behavior can be controlled by using the noflow boundary condition zone built into Walkabout.
- Uniform Material paraview and walkabout results nearly the same.
- Non-uniform Material paraview and walkabout differ because the vertices of a tetrahedron can lie on a material interface between materials with strongly contrasting hydraulic properties. The cell centers are within a single material, 
velocities interpolated onto vertices that lie along an interface are ill-defined. Doubly defined nodes at interfaces may be a solution to this issue, but requires code development and testing.
- A mesh that is non-convex non-Delaunay will have many negative coefficients that can impact results. Use a median .stor coefficient file as a work around.  Median stor volumes and coefficients are always positive, but the geometry is not guaranteed to be orthogonal. For median version of .stor file, see [Voronoi](https://github.com/lanl/voronoi).


## METHOD 3 
### plot_to_tet5 uses Amanzi hex mesh copied on to cell center median tet mesh 

This method uses LaGriT to create a tet mesh from the hex plot mesh by connecting the plot mesh cell centers. The cell centers become the tet vertices. The original mesh boundary nodes are lost in this conversion. 

- This method avoid interpolation of Amanzi values from cell center to cell vertices.
- This method has the best results.
- Non-convex non-Delaunay mesh results in many negative coefficients that can impact results. Use a median .stor file as a work around, all volumes and coefficients are positive, but the geometry is not guaranteed to be orthogonal.


See particleDriver and workflow development for Amanzi-Walkabout models:

   test repo at https://github.com/amanzi/amanzi-walkabout

