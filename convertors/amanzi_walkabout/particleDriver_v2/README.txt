METHOD 2 particleDriver V2 
This is a modification of particleDriver V1 using Amanzi HDF5 data file.
Amanzi mesh and walkabout mesh the hex cells converted into 5 tets using LaGriT.
This is a work around to avoid for a non-convex mesh and Delaunay is not guaranteed.
Use new July 2018 version Amanzi to write velocities to cell nodes instead of cell centers.

Amanzi will write walkabout files that are processed by particleDriver
Amanzi will not output a mesh, tet5.inp is used instead.
Amanzi will write the velocities interpolated to hex to tet5 nodes.

LaGriT lagrit.lanl.gov and https://github.com/lanl/LaGriT
Voronoi https://lanl.github.io/voronoi/
Walkabout https://github.com/lanl/walkabout

METHOD -----------------------------------------------

IMPORTANT: tet5.inp (required) must have same node order as sorted Exodus mesh.

NOTE: This requirement is no longer necessary with new lagrit interpolation speed-up.
This method as is, required tet5.inp and Amanzi output to have the same order.
This is still a good practice, but not necessary if lagrit is used for interpolation.

Check visually or with lagrit:

tet5_presort.inp (before sort)
07055000   5.450000000000E+05  4.135000000000E+06  1.418378897951E+03
07055001   5.451000000000E+05  4.135000000000E+06  1.417707782742E+03
07055002   5.452000000000E+05  4.135000000000E+06  1.417019598995E+03
07055003   5.453000000000E+05  4.135000000000E+06  1.416318913546E+03
07055004   5.454000000000E+05  4.135000000000E+06  1.415593589167E+03
07055005   5.455000000000E+05  4.135000000000E+06  1.414905005206E+03

tet5_flts_wtrfit_exo.inp (sorted)
07055000   5.359000000000E+05  4.134500000000E+06  1.462898881611E+03
07055001   5.365000000000E+05  4.134900000000E+06  1.462932470123E+03
07055002   5.362000000000E+05  4.134700000000E+06  1.462967944515E+03
07055003   5.357000000000E+05  4.134400000000E+06  1.463008127633E+03
07055004   5.355000000000E+05  4.134300000000E+06  1.463088145057E+03
07055005   5.366000000000E+05  4.135000000000E+06  1.463123825919E+03


Step 1)
Write tet5 version of the hex mesh used for Amanzi
Output: tet5_presort.inp and tmp.stor
  lagrit < hextotet.lgi
  cp outx3dgen hextotet.out.txt

Step 2)
Color the tet5 mesh
Output: tet5_mat_zones.inp
  lagrit < set_materials.lgi
  cp outx3dgen set_materials.out.txt


Step 3)
Write Amanzi tet5 mesh with facesets
Output: tet5_flts_wtr.exo and sorted AVS file tet5_flts_wtr_exo.inp
  lagrit < write_exo.lgi 
  cp outx3dgen write_exo.out.txt


Step 3)
Write FEHM files for Walkabout (Used after Amanzi runs)
Output: walkabout.tet5.fehmn walkabout.tet5.stor walkabout.tet5.graph walkabout.tet5_material.zone 
Output: tet5.inp (sorted mesh without attributes for walkabout from Amanzi)
Note these must have same node order as the exo file.
Check for postive element and voronoi volumes, there will be neg ccoefs on outside nodes.
  lagrit < tet_to_fehm.lgi

  outx3dgen is very large with 12598917 neg coef lines, remove these and save.
  awk '{if ($4 != "row" && $1 != "Negative") print $0}' outx3dgen > tet_to_fehm.out.txt


Step 4)
Write user defined outside zones for Walkabout
where keyword top = noflow zones
Output: outside_top.zone and outside_top_back.zone
  lagrit < get_outside_zones.lgi


Step 5)
The voronoi stor file has many neg coefficients as expected from tet5 mesh with non-planar boundaries and interfaces.
Write new stor file with median volumes more appropriate for non-Delaunay mesh:
  voronoi -avs tet5.inp -cv median -o tet5_median.stor -d

The output report for this voronoi stor file from lagrit:

*** Construct and Compress Sparse Matrix:3D ***                                 
AMatbld3d_stor: Total Number of Negative Coefficients  12730985                 
AMatbld3d_stor: Number of Significant Negative Coefs  12730985                  
AMatbld3d_stor: npoints =  7055100  ncoefs =   88523102                         
AMatbld3d_stor: Volume min =   2.8171515E+04                                    
AMatbld3d_stor: Volume max =   2.4928255E+05                                    
AMatbld3d_stor: Total Volume:   1.6217985E+12                                   

The output report for this median stor file from Voronoi:

 COEFFICIENTS
 =================================================================
    Min. Voronoi area:                               0.2077724E+01
    Max. Voronoi area:                               0.4154585E+06
    --------------------------------------------------------------
    Min. Voronoi edge length:                        0.1126861E+02
    Max. Voronoi edge length:                        0.1611800E+04
    --------------------------------------------------------------
    Min. area/len coefficient:                       0.0000000E+00
    Max. area/len coefficient:                       0.1478468E+03
    --------------------------------------------------------------
    Total Voronoi area:                              0.1721339E+13
    Total triangle face area:                        0.1495281E+12




