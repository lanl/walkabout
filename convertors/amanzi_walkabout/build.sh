g++ -I/usr/include/hdf5/serial -c h5output.cpp
g++ -I/usr/include/hdf5/serial -o particleDriver h5output.o -L/usr/lib/x86_64-linux-gnu/hdf5/serial /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_cpp.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5_hl.a /usr/lib/x86_64-linux-gnu/hdf5/serial/libhdf5.a -lstdc++ -lpthread -lhdf5 -ldl -lz -lsz

