F95=gfortran 
FFLAGS = -O3  -static 

OBJS = precision.o io.o ll.o cell_class.o mesh_class.o tracktools.o bftensor.o isotropictensor.o lichtnertensor.o disptensor.o flowfield_class.o walkabout.o  

walkabout: $(OBJS) 
	$(F95) $(FFLAGS) $(OBJS) -o walkabout  

precision.o: precision.f90  
	$(F95) $(FFLAGS) -c precision.f90   

ll.o: ll.f90 precision.o 
	$(F95) $(FFLAGS) -c ll.f90 

io.o: io.f90 
	$(F95) $(FFLAGS) -c io.f90 
 
cell_class.o: cell_class.f90 precision.o 
	$(F95) $(FFLAGS) -c cell_class.f90 

mesh_class.o: mesh_class.f90 cell_class.o io.o precision.o 
	$(F95) $(FFLAGS) -c mesh_class.f90 

tracktools.o: tracktools.f90 mesh_class.o cell_class.o precision.o 
	$(F95) $(FFLAGS) -c tracktools.f90 

bftensor.o: bftensor.f90 precision.o   
	$(F95) $(FFLAGS) -c bftensor.f90  

isotropictensor.o: isotropictensor.f90 precision.o   
	$(F95) $(FFLAGS) -c isotropictensor.f90  

lichtnertensor.o: lichtnertensor.f90 precision.o  
	$(F95) $(FFLAGS) -c lichtnertensor.f90  

disptensor.o: disptensor.f90 bftensor.o lichtnertensor.o isotropictensor.o precision.o  
	$(F95) $(FFLAGS) -c disptensor.f90  

flowfield_class.o: flowfield_class.f90 disptensor.o tracktools.o mesh_class.o io.o precision.o 
	$(F95) $(FFLAGS) -c flowfield_class.f90 

walkabout.o: walkabout.f90 tracktools.o flowfield_class.o mesh_class.o cell_class.o ll.o io.o precision.o  
	$(F95) $(FFLAGS) -c walkabout.f90 

pclite: pclite.o 
	$(F95) $(FFLAGS) pclite.o -o pclite  

pclite.o: pclite.f90   
	$(F95) $(FFLAGS) -c pclite.f90  

clean:
	rm *.o
	rm *.mod
	rm walkabout  

