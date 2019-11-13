# Appendix A.4  Output


Two output files are always produced, a log file that echoes back input parameters and a reduced **`sptr2`** file for PLUMECALC. The latter may be read directly by PLUMECALC.

In addition, to the two files always produced, a trajectory file **`trajout`** may be produced if the **toutfreq** parameter is set greater than 0.
See [`Appendix A.1 Control`](appendix-A1.md).


Walkabout and SPTR use different streamline tracing algorithms. The good agreement between the two thus helps build confidence in bothe codes.
This is done in Verification Test 3 comparing **`sptr2`** to **`trajout`** results, see [`Verification Test 3, section 4.3`](4-verification-tests.md)


## A.4.1 Plumecalc particle tracks **`sptr2`**

This streamline particle tracking file is written in FEHM sptr2 ASCII format and can be read by PLUMECALC.
The number of particles used in the simulation, followed by the particle number, time that the particle is leaving a cell (days), and the cell that the particle is leaving, for each travel segment of each particle. 
The example in Figure A4-1 is walkabout.sptr2 written for Test 3 (file walkabout.sptr2.ans is a copy).

```
      walkabout 1.2        20191008170857.604 
 dtmax 10.0                                                                      
          16
  Part_no    time_days     cell_leaving
        1    15.042780           7022
        1    73.054332           7023
        1    105.27757           7024
        1    129.20500           7050
        1    185.07995           7051
        1    206.95252           7052
(...)
       16    341.17299          10471
       16    347.76572          10472
       16    359.72616          10446
       16    380.07722          10447
       16    401.86738          10448
       16    425.39493          10449
       16    449.87253          10450
       16    475.15073          10425
       16    487.42966          10426
```
*Figure A4-1. This example is the file walkabout.sptr2 written for Test 3.*

## A.4.2 Walkabout particle tracks **`trajout`**


The walkabout particle file format is as follows.

```
walkabout title ! header line
[dtmax or tmax] max_step       ! max time step
npart           ! number particle tracks
nsteps          ! number of steps for this track
    t1 x1 y1 z1 ! spatial position at time t1
    : ! repeat above line for a total of nsteps
 : ! repeat each particle block for a total for npart particles
```

The example in Figure A4-2 is traj.out written for Test 3 (file traj.out.ans is a copy).
The max time steps is 10. The file contains 16 particle tracks.
The first particle track has 1355 lines with time, x, y, z
The  next particle track has 1284 lines with time, x, y, z
This is repeated for each of the 16 particle tracks.


```
      walkabout 1.2            20160727132144.729 
 dtmax 10.0                                                                      
          16
        1355
  0.0000000      5.0000000      40.000000      40.000000    
 0.23100000E-01  5.0015968      40.000247      40.000000    
 0.51051000E-01  5.0035288      40.000545      40.000000    
 0.84871710E-01  5.0058666      40.000907      40.000000    
(...)
 918.26953      99.779548      16.224802      39.999271
 919.48875      99.908511      16.224545      39.999271
 920.09356      99.972493      16.224516      39.999271
        1284
 0.0000000      5.0000000      60.000000      40.000000
 0.23100000E-01  5.0027597      60.000282      40.000000
 0.51051000E-01  5.0060991      60.000623      40.000000
 0.84871710E-01  5.0101397      60.001036      40.000000
 0.12579477      5.0150290      60.001537      40.000000
(...)

```
*Figure A4-2. This example is the file traj.out written for Test 3.*




