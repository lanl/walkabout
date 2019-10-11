
# Walkabout and FEHM

Walkabout Version 1.0 (LA-CC-11-033) performs random walk particle tracking simulations of solute transport based on groundwater flow solutions that use fully unstructured control volume grids. Walkabout V1.0 is designed to work within the [FEHM](https://fehm.lanl.gov) code system, and accepts groundwater flow solutions from FEHM and computational mesh descriptions from LaGriT (Los Alamos Grid Toolbox). Walkabout also provides input to the [PLUMECALC](https://plumecalc.lanl.gov) particle tracking post processor. Although designed to work within the FEHM system, other control-volume flow codes and mesh generators could, with appropriate reformatting of the output, provide the flow fields and mesh descriptions.

A typical workflow for Walkabout within the FEHM system would use LaGrit to generate unstructured grids. FEHM then provides a discretized representation of the steady-state flow field to Walkabout. Given this discrete solution, Walkabout then reconstructs a groundwater flow field, and performs the random walk particle tracking calculation. 

The following FEHM files are used by Walkabout, see the User Manual for a full description of input and output files. 


## Input Mesh Geometry

These files describing the mesh geometry are required. LaGrit can be used to write these files, although any grid generation software will work provided the output is converted to the required format. [See LaGriT Documentation](http://lanl.github.io/LaGriT)


* **fehmn** - mesh geometry file with node locations and element connectivity list.

* **stor** - sparse matrix geometric coefficient file in FEHM file format.
  1) Sparse matrix structure of geometric coefficients is defined by the face graph of the Voronoi tessellation or equivalently the edge graph of the Delaunay dual of the Voronoi tessellation.
  2) Volume of each Voronoi polygon - to fill the Diagonal entries of the sparse matrix
  3) Areas of each Voronoi polygon face - to fill the off diagonal entries of the sparse matrix where area_ij = area_ji so only half of the matrix is written.
  
*  **graph** or **ealist** - mesh element adjacency list

* **cbound** file is optional and provides a list of nodes on boundaries that are closed to transport using zone format.


## Input FEHM Result Files
  
* **fin** FEHM restart file with internodal liquid fluxes

* **avs** - AVS format file written by FEHM with mesh node properties porosity, liquid saturation and liquid density.


## Output Files

Two output files are always produced, a log file that echoes back input parameters and a reduced SPTR2 file for PLUMECALC. The latter may be read directly by PLUMECALC.
 
In addition, to the two files always produced, a trajectory file (trajout) may be produced if the toutfreq parameter is set greater than 0.

 
  
  
