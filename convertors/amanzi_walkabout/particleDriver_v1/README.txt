METHOD 1 particleDriver V1 
Originally developed for Amanzi to Walkabout conversions using Amanzi HDF5 data file.
Output includes data and points at Amanzi cell centers and the boundary nodes.
This particleDriver uses LaGriT to connect the points into a tetrahedral mesh for Walkabout.
 
This works ok for a box shaped (convex) mesh, but does not work for non-convex, non-Delaunay mesh.
Any work around will result in missing tet elements, or inverted tet elements along the boundary.
Walkabout Particles stop or have undefined behavior at these locations.

Other Methods are under development.

