# Appendix A.4  Output Files 


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


