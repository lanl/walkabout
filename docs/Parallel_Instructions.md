# Intructions for using openMP in walkabout

Walkbout will use the maximum number of available threads by default. To change the number of availble threads:

export OMP_NUM_THREADS=1

This value can be between 1 and max threads on your machine; 1 will run the code in serial, yet still load openMP libraries.
