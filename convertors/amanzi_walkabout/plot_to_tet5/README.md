# This Workflow uses Amanzi h5 hex plot and data files instead of interpolation methods  

<b>Input:</b><br>
- h5_plot_output.cpp (code to read plot data and write files)
- plot_mesh.h5 (hex mesh from current Amanzi run)
- plot_data.h5 (hex mesh data from current Amanzi run)
- tet5_median_clrs.inp (hextotet5 median mesh for all runs)

<b>Output:</b><br>
- Amanzi hex mesh connected into 5 tets with vertices at cell center points.
- Walkabout files with data copied from Amanzi hex centers to tet vertices
- Walkabout .stor median version instead of voronoi


This work is a continuation of efforts to integrate Amanzi cell center modeling with Walkabout using tet vertices.
This method is used for non-convex hex meshes. Non-convex point distributions results in poorly formed tets across boundaries. Removing bad tets did not work as the removal created poor behavior on the boundary where divots formed.


## Description

We use the Amanzi h5 plot and data files written on hex cell centers, these same centers used in Amanzi runs.
These are copied on to the tet5 median mesh vertices, there is therefore no interpolation of values.
Since this tet5 mesh has tet vertices formed from the Amanzi hex centers the hex boundary nodes are not included.
Therefore the Walkabout mesh is smaller on all sides by half a cell width.

- This is a workaround to avoid the issues from a non-convex mesh, non-Delaunay mesh.
- There is no interpolation of Amanzi values, plot data are copied directly the to tet5 vertices formed by the Amanzi hex cell centers (median).
- The hex to tet5 mesh is not quaranteed to be Delaunay and may have many negative coefficients in the voronoi stor file.
- The solutions from a mesh with negative coefficients may be stable but inaccurate. Particles encountering negative coefficients may have poorly defined behavior. 
- A workaround for the neg coefficients is to use median stor file. This will ensure there are no negative coefficients, but the geometry may not be orthogonal. 

## Software for this Method

LaGriT lagrit.lanl.gov and https://github.com/lanl/LaGriT

Voronoi https://lanl.github.io/voronoi/

Walkabout https://github.com/lanl/walkabout


## Workflow Development 

New particleDriver and workflows with Amanzi to Walkabout is under development.
See test repo at https://github.com/amanzi/amanzi-walkabout

