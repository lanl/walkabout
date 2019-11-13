# Appendix A.1  Control

The **`control`** file contains two keyword blocks that allow the user to specify dispersion tensors and particle start locations. In addition, it contains numerical parameters that control the functioning of the particle tracking simulation.

## A.1.1   INITIAL keyword block

Particle starting positions are specified in the **INITIAL** keyword block. The keyword block may be placed anywhere in the **`control`** file. The form of the **INITIAL** keyword block is as follows

```
INITIAL ! case sensitive
dist_type
num_part ! conditional input depends on dist_type
location_info ! conditional input depends on dist_type
:   ! repeat as needed
```

The *dist_type* keyword specifies how particles are to be distributed. Allowed values are **MANUAL**, **UNIFORM**, and **RANDOM** (case sensitive).

If **MANUAL** is specified then a particle location must be specified for each particle. In this case, the number of particles should be entered in the **num_part** field and start coordinates for each particle should be entered (**location_info** field, one x,y,z triplet on each line).

If **RANDOM** is specified, the particles are to be distributed randomly in a user- specified rectangular box. In this case, the number of particles should be entered in the *num_part* field and the *location_info* has two lines representing locations of minimum (first line) and maximum (second line) coordinates of the bounding box for the start region.

The **UNIFORM** keyword is similar to the **RANDOM** keyword except that particles will be distributed uniformly in the region of interest. In this case, three integers should be entered in the *num_part* representing spacing in the x, y, and z directions. The *location_info* lines are identical to the **RANDOM** case.

## A.1.2   DTENSOR keyword block

Dispersion tensor information is given in the **DTENSOR** keyword block. The keyword block may be placed anywhere in the **`control`** file. The form of the **INITIAL** keyword block is as follows

```
DTENSOR
region_specifier !
type ! tensor type
dispersivities ! conditional on type
  : ! repeat above three lines as many times as needed
END ! terminates the keyword block
```

The *region_specifier* comes in two forms. If three integers (*min max stride*) are given, then the region comprises the nodes between *min* and *max* (inclusive) with the specified stride. The value 0 for *max* is interpreted as the last node in the grid. This format is identical to FEHM’s format for specifying nodes in a region. If *region_specifier* starts with a nonnumeric character, then it is interpreted as a filename that contains zone information. In this case, the first field of the *region_specifier* is the filename and the second field of the *region_specifier* is the zone number in the specified zone file. The zone file should be in FEHM zone file format (See Appendix A.2 Mesh Geometry Files).

The *type* holds a keyword that specifies the type of the dispersion tensor. In Version 1.0, the only allowed value is **BF** denoting the Burnett and Frind tensor. The dispersivity_values required for the Burnett and Frind tensor are (in order) longitudinal dispersivity, horizontal transverse dispersivity, vertical transverse dispersivity, and molecular diffusion coefficient. The dispersivities have units of
m. The molecular diffusion coefficient has units of \(m^2/\mathrm{day}\).


## A.1.3   Numerical control parameters

Several numerical control parameters may also be specified in the **`control`** file. If not defined in the **`control`** file, defaults will apply. 
The parameters are scalar values and may be written in any order. The form is

*`parameter_name parameter_value`*

where *`parameter_name`* is one of the following



* **dtmax** – maximum allowed value of time step in days. Optional. Defaults to 1000

* **dt0** – initial time step in days. Optional. Defaults to 0.01

* **maxstretch** – maximum allowed time step relative to previous time step (dimensionless). Must be greater than 1. Optional. Defaults to 1.2

* **maxsteps** – maximum number of time steps in the simulation. A particle is terminated if maxstep timesteps are taken. Optional defaults to 100000

* **dxtarget** – Courant factor, maximum allowed time step relative to time required to advect across the cell. Optional. Defaults to 0.1

* **dttarget** – maximum allowed time step as fraction of characteristic time to disperse across cell. Optional. Defaults to 0.1

* **toutfreq** – number of time steps between output of trajectory data. If 0 is specified, no trajectory output is written. Optional. Defaults to 0



## A.1.4   Example **`control`** files



In this example shown in Figure A1-1, 10000 particles are released randomly in the region between (10,-50,-50) and (10,50,50). The region is indicated by 1 0 0 which is all grid nodes.  The dispersivity_values are 40. for longitudinal dispersivity and 0. for horizontal transverse dispersivity, vertical transverse dispersivity, and molecular diffusion coefficient. 

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
*Figure A1-1. Example **`control`** file.*

In this example shown in Figure A1-2,  different transport parameters are applied to different regions of the grid using FEHM style zone files.
In this example,  there are 5 node regions (zone numbers 4, 9, 11, 12, and 14) in the model domain and are defined in walk06337.h5_material.zone. 
The zone file format is described under **`cbound`** in [Appendix A.2 Geometry](appendix-A2.md). 

```
dtmax 365.25
dt0 0.10
maxstretch 1.3
maxsteps 100000
dxtarget 0.1
dttarget 0.1
seed 7127
toutfreq 0
INITIAL
MANUAL
5
   8200.2050        30.0000      1812.0000
   8200.2050        30.0000      1826.0000
   8200.6150        30.0000      1840.0000
   8200.6150        30.0000      1854.0000
   8200.6150        30.0000      1868.0000
 
DTENSOR
walk06337.h5_material.zone 4
BF
1.0  0.1  0.1  0.0
walk06337.h5_material.zone 9
BF
3.0  0.3  0.1  0.0
walk06337.h5_material.zone 11
BF
10.0  1.0  1.0  0.0
walk06337.h5_material.zone 12
BF
2.0  0.1  0.1  0.0
walk06337.h5_material.zone 14
BF
5.0  0.5  0.5  0.0
END
```
*Figure A1-2. Example **`control`** file with multiple regions.*

