# particleDriver V1 – Amanzi - Walkabout conversion (Sept 2017 and earlier)

particleDriver is a C++ program for a particle advection system using Amanzi, LaGriT, and Walkabout. It automates the process of producing input and output files for LaGrit and Walkabout from an Amanzi checkpoint file and creates data for visualization that can be used in Paraview, PlumeCalc, and Meshlab.

- Originally developed for Amanzi to Walkabout conversions using Amanzi HDF5 data file.
- Amanzi output includes data and points at Amanzi cell centers and the boundary nodes.
- This particleDriver uses LaGriT to connect the points into a tetrahedral mesh for Walkabout.
 
This works ok for a box shaped (convex) mesh, but does not work for non-convex, non-Delaunay mesh.
Any work around will result in missing tet elements, or inverted tet elements along the boundary.
Walkabout Particles stop or have undefined behavior at these locations.


## Workflow Development 

New particleDriver and workflows with Amanzi to Walkabout is under development.
See test repo at https://github.com/amanzi/amanzi-walkabout

