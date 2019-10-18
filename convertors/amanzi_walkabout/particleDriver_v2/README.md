# particleDriver V2 Amanzi(V.88+) - Walkabout conversion  

This is a modification of particleDriver V1 using Amanzi HDF5 data file.
The Amanzi mesh and the walkabout mesh are the the hex cells converted into 5 tets using LaGriT.
Use new July 2018 version Amanzi to output velocities to cell 5 tet vertices instead of cell centers.

- This is a workaround to avoid the issues from non-convex mesh, non-Delaunay mesh.
- Amanzi will not output a mesh, the input mesh tet5.inp is used instead.
- Amanzi will write the velocities interpolated from hex to 5 tet vertices.
- The hex to tet5 mesh is not quaranteed to be Delaunay and may have many negative coefficients in the stor file.
- The solutions from a mesh with negative coefficients may be stable but inaccurate. Particles encountering negative coefficients may have poorly defined behavior. A median stor file will ensure there are no negative coefficients, but the geometry may not be orthogonal. 

LaGriT lagrit.lanl.gov and https://github.com/lanl/LaGriT
Voronoi https://lanl.github.io/voronoi/
Walkabout https://github.com/lanl/walkabout


## Workflow Development 

New particleDriver and workflows with Amanzi to Walkabout is under development.
See test repo at https://github.com/amanzi/amanzi-walkabout

