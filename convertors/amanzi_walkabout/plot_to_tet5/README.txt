METHOD 3 Workflow uses Amanzi hex plot and data files instead of HDF5 data.

This method uses the median points of the Amanzi hex mesh, connected by hex to tet5.
We use the Amanzi h5 plot and data files written on hex cell centers, these same centers used in Amanzi runs.
These are copied on to the tet5 median mesh, there is therefore no interpolation of values.
Since this tet5 mesh has tet vertices formed from the Amanzi hex centers the hex boundary nodes are not included.
Therefore the Walkabout mesh is smaller on all sides by half a cell width.

LaGriT lagrit.lanl.gov and https://github.com/lanl/LaGriT
Voronoi https://lanl.github.io/voronoi/
Walkabout https://github.com/lanl/walkabout

Use run_plot.scr or follow steps:

Write script for median mesh from Amanzi plot mesh cell center data.
This assumes the hex and tet median mesh have been created.
See directory lagrit_files if this needs to be done.

STEP 1
Copy from Amanzi directory:
cp /scratch/fwo/zhiming/PM_ASCEM/multiwells/inverse_results_from_wolf/newFluxConstrain_grizzly/iteration4_hex_uniformK/plot_data.h5 .
cp /scratch/fwo/zhiming/PM_ASCEM/multiwells/inverse_results_from_wolf/newFluxConstrain_grizzly/iteration4_hex_uniformK/plot_mesh.h5 .

STEP 1a
If median tet5 mesh has been created:
  ln -s ../hex_to_median/walkabout.h5_med_tet.inp tet5_median.inp

STEP 1b
If median mesh does not have materials, add them.
Interpolate Amanzi hex mesh materials to median nodes.
  lagrit < hex_materials_to_median.lgi

STEP 2
Read h5 data and interpolate on to tet5 mesh.
  ./build.sh
  ./h5_output

STEP 3 (final files)
Read hex_plot_mesh.inp and tet5 median mesh to write walkabout files.
  lagrit < plot_to_walkabout.lgi

STEP 3a
  # write .ama file for walkabout
  # edit walkabout.h5_med.ama header to be number of nodes
  # Note the first column with node number is ignored
      6930000

  # write .avs file for walkabout
  # header should look like this
  # 03 1   1   1
  # Pressure (MPa), (MPa)
  # Saturation, (no dim)
  # Porosity, (no dim)


STEP 3b (stor file)
If this median mesh needs a median stor file instead of voronoi stor file to avoid large number of neg coefficients
that result from a non-delaunay connected mesh.

   voronoi -avs tet5_at_zero.inp  -cv median -o tet5_median.stor -d


NOTES:

To create c++ text for writing the input .lgi file:
awk -f myfile.awk input.lgi > input_lgi.txt




