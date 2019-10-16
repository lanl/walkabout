# Appendix A. Input and Output Formats

## A.1 The control file

The control file contains two keyword blocks that allow the user to specify dispersion tensors and particle start locations. In addition, it contains numerical parameters that control the functioning of the particle tracking simulation.

### A.1.1   *INITIAL* keyword block

Particle starting positions are specified in the *INITIAL* keyword block. The keyword block may be placed anywhere in the **control** file. The form of the *INITIAL* keyword block is as follows

```
INITIAL ! case sensitive
dist_type
num_part ! conditional input depends on dist_type
location_info ! conditional input depends on dist_type
:   ! repeat as needed
```

The dist_type keyword specifies how particles are to be distributed. Allowed values are *MANUAL*, *UNIFORM*, and *RANDOM* (case sensitive).

If *MANUAL* is specified then a particle location must be specified for each particle. In this case, the number of particles should be entered in the **num_part** field and start coordinates for each particle should be entered (**location_info** field, one x,y,z triplet on each line).

If *RANDOM* is specified, the particles are to be distributed randomly in a user- specified rectangular box. In this case, the number of particles should be entered in the **num_part** field and the **location_info** has two lines representing locations of minimum (first line) and maximum (second line) coordinates of the bounding box for the start region.

The *UNIFORM* keyword is similar to the *RANDOM* keyword except that particles will be distributed uniformly in the region of interest. In this case, three integers should be entered in the **num_part** representing spacing in the x, y, and z directions. The **location_info** lines are identical to the *RANDOM* case.

### A.1.2   DTENSOR keyword block

Dispersion tensor information is given in the *DTENSOR* keyword block. The keyword block may be placed anywhere in the **control** file. The form of the *INITIAL* keyword block is as follows

```
DTENSOR
region_specifier !
type ! tensor type
dispersivities ! conditional on type
  : ! repeat above three lines as many times as needed
END ! terminates the keyword block
```

The **region_specifier** comes in two forms. If three integers (**min max stride**) are given, then the region comprises the nodes between **min** and **max** (inclusive) with the specified stride. The value 0 for **max** is interpreted as the last node in the grid. This format is identical to FEHM’s format for specifying nodes in a region. If **region_specifier** starts with a nonnumeric character, then it is interpreted as a filename that contains zone information. In this case, the first field of the **region_specifier** is the filename and the second field of the **region_specifier** is the zone number in the specified zone file. The zone file should be in LaGrit’s format.

The **type** keyword is specifies the type of the dispersion tensor. In Version 1.0, the only allowed value is *BF* denoting the Burnett and Frind tensor. The dispersivity_values required for the Burnett and Frind tensor are (in order) longitudinal dispersivity, horizontal transverse dispersivity, vertical transverse dispersivity, and molecular diffusion coefficient. The dispersivities have units of
m. The molecular diffusion coefficient has units of \(m^2/\mathrm{day}\).

### A.1.3   Numerical control parameters

Several numerical control parameters may also be specified in the control file. Defaults exist for each parameter. The parameters take scalar values may come in any order. The form is

**`parameter_name parameter_value`**

where **parameter_name** is one of the following

* *dtmax* – maximum allowed value of time step in days. Optional. Defaults to 1000

* *dt0* – initial time step in days. Optional. Defaults to 0.01

* *maxstretch* – maximum allowed time step relative to previous time step (dimensionless). Must be greater than 1. Optional. Defaults to 1.2

* *maxsteps* – maximum number of time steps in the simulation. A particle is terminated if maxstep timesteps are taken. Optional defaults to 100000

* *dxtarget* – Courant factor, maximum allowed time step relative to time required to advect across the cell. Optional. Defaults to 0.1

* *dttarget* – maximum allowed time step as fraction of characteristic time to disperse across cell. Optional. Defaults to 0.1

* *toutfreq* – number of time steps between output of trajectory data. If 0 is specified, no trajectory output is written. Optional. Defaults to 0

### A.1.4   Example control file

Examples of control files are given in Figures A-x and A-y. In the example in Figure A-x, 10000 particles are released randomly in the region between (10,-50,-50) and (10,50,50). The dispersion coefficient in this example is

```
Title – the title goes here
dtmax 100   !maximum step size days
dxtarget 0.1    !relative to grid size
dttarget 0.1    ! relative to characteristic dispersion time 
maxstretch 1.3
maxsteps 100000 ! maximum number of steps allowed
dt0 0.1   ! initial step days
toutfreq 0  !no trajectory output

INITIAL
 RANDOM
 10000
 10.0 -50.0 -50.0
 10.0    50.0    50.0

DTENSOR
1 0 0
BF  !burnett and frind tensor
 40. 0.0 0.0 0.0 !dispersivity values
END
```

*Figure A-1. Example control file.*

```
Title Goes Here
dtmax 100   !maximum step size days
dxtarget 0.1    !relative to grid size
dttarget 0.1    ! relative to characteristic dispersion time 
maxstretch 1.3
maxsteps 100000 ! maximum number of steps allowed
dt0 0.1   ! initial step days
toutfreq 0  !no trajectory output

INITIAL
 RANDOM
 10000
 10.0 -50.0 -50.0
 10.0    50.0    50.0

DTENSOR
1 0 0
BF  !burnett and frind tensor
 40. 0.0 0.0 0.0 !dispersivity values
END
```

*Figure A-2. Example control file.*

## A.2 Files describing geometry of the mesh

Three files describing geometry of the computational grid are required, and a fourth is optional. It is anticipated that these files will be produced by the LaGrit software, although any grid generation software could be used provided the output is converted to the required format.

### A.2.1  The fehmn file

The **fehmn** file provides geometry information (location of nodes and lists of nodes that compose each element). See the LaGrit and FEHM manuals for details. In Walkabout Version 1.0, the **fehmn** file must be provided in ASCII format.

### A.2.2   The stor file

The **stor** file provides information about nodal connectivity and interface areas. See the LaGrit and FEHM manuals for details. LaGrit options exist to produce vector areas, scalar areas, or ratio of scalar area to distance for each node-to- node connection. The latter option is required by Walkabout. In Walkabout Version 1.0, the **stor** file must be provided in ASCII format.

### A.2.3   The ealist file

The **ealist** file provides information about element adjacency, which should not be confused with nodal connectivity. Element adjacency is not needed by FEHM, but is required by Walkabout. See the LaGrit manual for details on how to produce the element adjacency lists.

### A.2.4   The cbound file

The **cbound** (closed boundary) file provides a list of nodes on boundaries that are closed to transport. It has the same format as the LaGrit outside zone file, but excludes outflow boundaries. A simple strategy for producing the **cbound** file is to use LaGrit to produce a list of all outside nodes, and then remove those nodes associated with outflow boundaries.

In Version 1, external faces of cells on a boundary must aligned with the principal directions in the coordinate system. That is, boundary faces must be top, bottom, left, right, back or front. See the LaGrit manual. It is important to recognize, that this restriction only applies to cell faces on boundaries. Cell faces internal to the model have no such restrictions. A node/cell may have more than one no- transport boundary face associated with it, in which case it would appear more than once in the **cbound** file.

For nodes on a no-transport boundary, Walkabout first attempts to reconstruct the nodal velocity using the unconstrained algorithm, Eq. 5. If this procedure results in inflow into the domain at the boundary node, the reconstructed velocity is used as is. If the unconstrained procedure results in outflow on a no-transport boundary, then the velocity reconstruction is repeated using the constrained procedure Eq. 7 to enforce the no-flow condition on the cell’s boundary face.
Particles are not allowed to disperse across boundaries that are closed to transport.

## A.3 Files containing FEHM results

Two files containing required FEHM results are required.

Internodal liquid fluxes are read from the **fin** file, an FEHM restart file in ASCII format. If the liquid fluxes are not found in the fin file, Walkabout will terminate.

Porosity, liquid saturation and liquid density are read from an ASCII **avs** file produced by FEHM. The **avs** file is produced by the FEHM contour macro. See the FEHM manual for details.

## A.4 Output files

Two output files are always produced, a log file that echoes back input parameters and a reduced SPTR2 file for PLUMECALC. The latter may be read directly by PLUMECALC.

In addition, to the two files always produced, a trajectory file (**trajout**) may be produced if the *toutfreq* parameter is set greater than 0 (see section A.1). The **trajout** file has format

```
!header line
!header line
npart !number particles
 ntimes !number time steps reported for this particle
    t1 x1 y1 z1 ! spatial position at time t1
    : ! repeat above line for a total of ntimes
 : ! repeat each particle block for a total for npart particles
```

