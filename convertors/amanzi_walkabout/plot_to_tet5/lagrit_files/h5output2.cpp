// https://github.com/amanzi/amanzi-walkabout/blob/master/bin/h5output.cpp
//
// tamiller Jun 2018
// modified to read from mesh vertices, no cell centers
// assume the hextotet converted file is available
// Add error checks and screen output for user information
// Zhiming
// modified h5output.cpp for new hdf5 libs May 2018
// From  /scratch/nts/zhiming/RM_SM/walkabout_workflow/h5output.cpp May 2017
// This version does update include string to cstring

// TODO
// Need to check for Material exist or Material -1 so imt is positive in output mesh
// Note Walkabout does not read materials, they are folded into property values. This is for viewing purposes
// Perhaps use interpolate instead of copyatt to ensure node order Amanzi correlates to Walkabout mesh
 
#include <hdf5.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cstdlib>

using namespace std;

//dataStruct is a derived type that contains data and metadata from the .h5 checkpoint file
//num_properties are number of properties to write to the AVS grid files as attributes.
// current num_properties = 7 which are xv, yv, zv, por, pres, sat, material

typedef struct {
	double *x, *y, *z;
	double *xv, *yv, *zv;
	double minx, miny, minz;
	double maxx, maxy, maxz;
	double *por, *pres, *sat, *den, *por2;
        int *material;
        int num_properties;
	unsigned int size, size2;
	char *filename;
} dataStruct;

//parameterStruct is a derived type that contains information about what to run and where to find it
typedef struct {
	string lagrit, walkabout, amanzi;
	bool pre, runLagrit, runWalkabout, post;
} parameterStruct;

//function prototypes
double* readDoubleArray(hid_t file_id, const char* dset, int &size);
int* readIntegerArray(hid_t file_id, const char* dset, int &size);
dataStruct getData(char *filename);
void writePly(dataStruct data);
void writeInp(dataStruct data);
void writeAma(dataStruct data);
void writeAvs(dataStruct data);
void writeLgi(dataStruct data);
void writeControl(dataStruct data);
void writeWalkaboutFile(dataStruct data);
void writePlumeCalcFiles(dataStruct data);
void writeObj(dataStruct data);
void writePosVel(dataStruct data);
void runLagrit(dataStruct data, parameterStruct param);
void runWalkabout(dataStruct data, parameterStruct param);
void runAmanzi(dataStruct data, parameterStruct param);
bool readTets(dataStruct data, double** verts, int** faces, int &vCount, int& fCount);
void writeVTKs(dataStruct data);
void setOperations(parameterStruct &param, int argc, char* argv[]);
parameterStruct readConfig();

// MAIN main
// Controls the flow of the program by calling functions in the necessary order

int main(int argc, char* argv[]) {
	//get the locations of executables
	parameterStruct param = readConfig();
	//determine which operations the user wants to run
	setOperations(param, argc, argv);
	//get the data from the h5 file
	dataStruct data = getData(argv[1]);
	//do preprocessing steps if specified to do so
	if(param.pre){
		cout << "Running pre processing steps.\n";
		//create files needed to run LaGriT and Walkabout
		writePly(data);
		writeInp(data);
		writeAvs(data);
		writeLgi(data);
		writeAma(data);
		writeWalkaboutFile(data);
		writePlumeCalcFiles(data);
		writeControl(data);
	}
	//run LaGriT if specified to do so
	if(param.runLagrit){
		cout << "Running Lagrit.\n";
		runLagrit(data, param);
	}
	//run Walkabout if specified to do so
	if(param.runWalkabout){
		cout << "Running Walkabout\n";
		runWalkabout(data, param);
	}
	//do postprocessing steps if specified to do so
	if(param.post){
		cout << "Running post processing steps.\n";
		writeObj(data);
		writePosVel(data);
		writeVTKs(data);
	}
	//end the program
	return 0;
}

//getData creates a dataStruct populated with data from a specified h5 file
dataStruct getData(char *filename){
	int size;
	dataStruct data;
	hid_t file_id; 
	herr_t status;
	
	//open the h5 file
	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	
	//read in the walkabout data.  Default xyz are cell centered + vertices on boundaries
	//             Tet xyz are cell vertices, these properties will be copied to the grid 
	data.x 		= readDoubleArray(file_id, "x", size);
	data.y 		= readDoubleArray(file_id, "y", size);
	data.z 		= readDoubleArray(file_id, "z", size);
	data.xv 	= readDoubleArray(file_id, "pore velocity x", size);
	data.yv 	= readDoubleArray(file_id, "pore velocity y", size);
	data.zv 	= readDoubleArray(file_id, "pore velocity z", size);
	data.por 	= readDoubleArray(file_id, "porosity", size);
	
	//all the above arrays are the same size, so we save the size for future reference
	data.size 	= size;
	
	//initialize min/max locations in each dimension
	data.minx = data.x[0];
	data.maxx = data.x[0];
	data.miny = data.y[0];
	data.maxy = data.y[0];
	data.minz = data.z[0];
	data.maxz = data.z[0];
	
	//then loop through the data, finding the actual min and max in the x, y, and z dimensions
	for(int a = 1; a < size; a++){
		if		(data.x[a] < data.minx){data.minx = data.x[a];}
		else if	(data.x[a] > data.maxx){data.maxx = data.x[a];}
		if		(data.y[a] < data.miny){data.miny = data.y[a];}
		else if	(data.y[a] > data.maxy){data.maxy = data.y[a];}
		if		(data.z[a] < data.minz){data.minz = data.z[a];}
		else if	(data.z[a] > data.maxz){data.maxz = data.z[a];}
	}
	
	//then read in other datasets properties for LaGriT
	data.pres 	= readDoubleArray(file_id, "pressure", size);
	data.sat 	= readDoubleArray(file_id, "saturation", size);
	data.por2 	= readDoubleArray(file_id, "porosity", size);
        data.num_properties = 6;

        // material id is now included in the Amanzi h5 files
    	data.material 	= readIntegerArray(file_id, "material ids", size);
        data.num_properties = 7;
	
	//all these arrays are the same size, so we save this size too for future reference
	data.size2 = size;

	
	//save the filename for future reference
	data.filename = filename;
	
	//then close the h5 file
	status = H5Fclose(file_id);
	
	//and return the data
	return data;
}

//readDoubleArray is a helper function that gets a specified dataset from a specified h5 file
double* readDoubleArray(hid_t file_id, const char* dset, int &size){
	double *dset_data;
	herr_t status;
	hid_t dataset_id;
	//get the dataset ID
	dataset_id = H5Dopen(file_id, dset, H5P_DEFAULT);

	//get the number of values.  dividing by 8 because the size returned is in bytes, and a double is 8 bytes
	size = H5Dget_storage_size(dataset_id)/8;
	//allocate memory for the dataset
	dset_data = new double[size];
	//read data as doubles into the allocated memory

	status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
	//close the dataset
	status = H5Dclose(dataset_id);
	//return the data
	return dset_data;
}

//readIntegerArray is a helper function that gets a specified dataset from a specified h5 file
int* readIntegerArray(hid_t file_id, const char* dset, int &size){
        int *dint_data;
        int isize;
        herr_t status;
        hid_t dataset_id;
        int j;

        //get the dataset ID
        dataset_id = H5Dopen(file_id, dset, H5P_DEFAULT);

        //get the number of values.  dividing by 8 because the size returned is in bytes, and a double is 8 bytes
        isize = H5Dget_storage_size(dataset_id)/8;

	dint_data = new int[isize];

        status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, dint_data);

// DEBUG
// printf("Integer Dataset Size: %d\n",isize);
// for (j = 0; j < isize; j++) {
//    printf("%d %d\n ",j, dint_data[j]); }

        //close the dataset and return data
        status = H5Dclose(dataset_id);
        return dint_data;
}


//writePly writes a .ply file representing the points that make up the domain, using their velocities
//to specify the color.  These files can be opened in Meshlab and other point cloud software
void writePly(dataStruct data){
	ofstream myfile;
	string filename = data.filename;
	filename += ".ply";
	myfile.open (filename.c_str());
	
	//write out header information
	myfile << "ply\n";
	myfile << "format ascii 1.0\n";
	myfile << "element vertex " << data.size << "\n";
	myfile << "property float x\n";
	myfile << "property float y\n";
	myfile << "property float z\n";
	myfile << "property uchar diffuse_red\n";
	myfile << "property uchar diffuse_green\n";
	myfile << "property uchar diffuse_blue\n";
	myfile << "end_header\n";
	
	//initialize min/max variables
	double minxv, maxxv, minyv, maxyv, minzv, maxzv;
	minxv = maxxv = data.xv[0];
	minyv = maxyv = data.yv[0];
	minzv = maxzv = data.zv[0];
	
	//loop through velocities, finding the min and max in the x, y, and z dimensions
	for(int a = 0; a < data.size; a++){
		if(data.xv[a] < minxv)
			minxv = data.xv[a];
		else if(data.xv[a] > maxxv)
			maxxv = data.xv[a];
		if(data.yv[a] < minyv)
			minyv = data.yv[a];
		else if(data.yv[a] > maxyv)
			maxyv = data.yv[a];
		if(data.zv[a] < minzv)
			minzv = data.zv[a];
		else if(data.zv[a] > maxzv)
			maxzv = data.zv[a];
	}
	
	int r,g,b;
	//for each position, convert the velocities into RGB values relative to the min/max values for the data
	//the x range corresponds to the Red value, y to the Green value, and z to the Blue value.
	//then write out the position and corresponding RGB value
	for(int a = 0; a < data.size; a++){
		if(minxv != maxxv){
			r = int((data.xv[a] - minxv)/(maxxv - minxv) * 255.0);
		}else{
			r = 255;
		}
		if(minyv != maxyv){
			g = int((data.yv[a] - minyv)/(maxyv - minyv) * 255.0);
		}else{
			g = 255;
		}
		if(minzv != maxzv){
			b = int((data.zv[a] - minzv)/(maxzv - minzv) * 255.0);
		}else{
			b = 255;
		}
		myfile << data.x[a] <<  " " << data.y[a] << " " << data.z[a] << " " << r << " " << g << " " << b << "\n";
	}
	//then close the file
	myfile.close();
}

//writeInp writes AVS point .inp file, which is an input file for LaGriT to set tet vertices 
// vertices can be connected or copied into hextotet mesh
// add attributes to this AVS node file

void writeInp(dataStruct data){
        double xv, yv, zv;
	ofstream myfile;

	string filename = "walkabout.h5_node.inp";
	myfile.open (filename.c_str());
	
	//write out header information with number of attributes for the nodes
	myfile << data.size << " " << " 0 " << data.num_properties << "  0  0\n";
	
	//write node x y z positions 
	for(int a = 0; a < data.size; a++){
		myfile << (a+1) << " " << data.x[a] << " " << data.y[a] << " " << data.z[a] << "\n";
	}

	//write node attributes which are the arrays read from the h5 file 

        if(data.num_properties == 6) {
	  myfile << "06  1 1 1 1 1 1\n";
  	  myfile << "Pressure (MPa), (MPa)\n";
	  myfile << "Saturation, (no dim)\n";
	  myfile << "Porosity, (no dim)\n";
	  myfile << "Xv, (no dim)\n";
	  myfile << "Yv, (no dim)\n";
	  myfile << "Zv, (no dim)\n";
	
	  for(int a = 0; a < data.size2; a++){
                xv = data.xv[a];
                yv = data.yv[a];
                zv = data.zv[a];
myfile << (a+1) << "\t" << data.pres[a] << "\t" << data.sat[a] << "\t" << data.por2[a] << "\t" << xv << "\t" << yv << "\t" << zv << "\n";
	}}

        if(data.num_properties == 7) {
          myfile << "07 1 1 1 1 1 1 1\n";
          myfile << "Material, integer\n";
          myfile << "Pressure (MPa), (MPa)\n";
          myfile << "Saturation, (no dim)\n";
          myfile << "Porosity, (no dim)\n";
          myfile << "Xv, (no dim)\n";
          myfile << "Yv, (no dim)\n";
          myfile << "Zv, (no dim)\n";

          for(int a = 0; a < data.size2; a++){
                xv = data.xv[a];
                yv = data.yv[a];
                zv = data.zv[a];
myfile << (a+1) << "\t" << data.material[a] <<  "\t" << data.pres[a] << "\t" << data.sat[a] << "\t" << data.por2[a] << "\t" << xv << "\t" << yv << "\t" << zv << "\n";
        }}
	
	myfile.close();
}

//writeAma writes an ama file, which is an input for Walkabout containing cell velocities
void writeAma(dataStruct data){
	ofstream myfile;
//	string filename = data.filename;
	string filename = "walkabout.h5.ama";
	myfile.open (filename.c_str());
	
	//write out the number of data points
	myfile << data.size << "\n";
	
	double xv, yv, zv;
	//loop through points
	for(int a = 0; a < data.size; a++){
		//calculate a corrected velocity
		xv = data.xv[a];
		yv = data.yv[a];
		zv = data.zv[a];
		//write out corrected velocities
		myfile << a << "\t" << xv << "\t" << yv << "\t" << zv << "\n";
	}
	
	myfile.close();
}

//writeAvs writes an avs file, which is an input file for Walkabout containing 
//pressure, saturation, porosity
void writeAvs(dataStruct data){
	ofstream myfile;
//	string filename = data.filename;
	string filename = "walkabout.h5.avs";
	myfile.open (filename.c_str());
	
	//write out header information
	myfile << "03  1   1   1\n";
	myfile << "Pressure (MPa), (MPa)\n";
	myfile << "Saturation, (no dim)\n";
	myfile << "Porosity, (no dim)\n";
	
	//write out data
	for(int a = 0; a < data.size2; a++){
//		myfile << (a+1) << "\t" <<  data.pres[a] << "\t" << data.sat[a] << "\t" << data.por2[a] << "\t" << data.den[a] << "\n";
		myfile << (a+1) << "\t" <<  data.pres[a] << "\t" << data.sat[a] << "\t" << data.por2[a] << "\n";
	}
	
	myfile.close();
}

//writeLgi writes an lgi control file for LaGriT to use a hextotet connected tet mesh 
//mesh properties are added to this tet mesh
//This controls the input and output file names and controls what operations LaGritT should perform.
void writeLgi(dataStruct data){
	ofstream myfile;
	string filename = data.filename;
	filename += ".lgi";
	myfile.open (filename.c_str());
	
	//write out standard LaGriT lgi file
	myfile << "# Define input and output file names\n";
	myfile << "#\n";
	myfile << "define / INPUT /" <<  "walkabout.h5_node" << ".inp\n";
	myfile << "define / INPUT_TET /" <<  "tet5.inp\n";
	myfile << "define / OUT_AVS / " <<  "walkabout.h5_tet" << ".inp\n";
	myfile << "define / OUT_GMV / " <<  "walkabout.h5_tet" << ".gmv\n";
	myfile << "define / OUT_FEHM / " <<  "walkabout.h5" << "\n";
	myfile << "define / OUT_ADJ / " <<  "walkabout.h5" << ".graph\n";
	myfile << "#\n";
	myfile << "# Read xyz nodes and attributes\n";
        myfile << "read avs INPUT mopts\n";
        myfile << "cmo/setatt/mopts/imt/1\n";
        myfile << "read avs INPUT_TET motmp\n";
        myfile << "cmo/status/motmp \n";
        myfile << "# remove extra attributes and write\n";
        myfile << "dump avs2  tmp.inp motmp 1 1 0 0 0\n";
        myfile << "read avs tmp.inp motet\n";
        myfile << "cmo/delete/motmp\n";
        myfile << "cmo/status/motet\n";
	myfile << "#\n";

        // if material property exists, copy to imt for the tet mesh
	// DEBUG quick fix for -1 Materials from Amanzi
        if (data.num_properties == 7) {
          myfile << "cmo / setatt / motet / imt / 1 0 0 / 1\n";
          // myfile << "cmo/copyatt/ motet mopts/ imt Material\n";
        }
        else {
          myfile << "# Set default materials to 1\n";
          myfile << "cmo / setatt / motet / imt / 1 0 0 / 1\n";
        }
        myfile << "cmo / setatt / motet / itetclr / 1 0 0 / 1\n";

        myfile << "resetpts/itp\n";
        myfile << "# Replace copyatt with interpolate if node ordering differs from tet mesh\n";
        myfile << "cmo/addatt/motet/Pressure/VDOUBLE/scalar/nnodes/linear/permanent/gxaf/0\n";
        myfile << "cmo/addatt/motet/Saturation/VDOUBLE/scalar/nnodes/linear/permanent/gxaf/0\n";     
        myfile << "cmo/addatt/motet/Porosity/VDOUBLE/scalar/nnodes/linear/permanent/gxaf/0\n";
        myfile << "cmo/addatt/motet/Xv/VDOUBLE/scalar/nnodes/linear/permanent/gxaf/0\n";
        myfile << "cmo/addatt/motet/Yv/VDOUBLE/scalar/nnodes/linear/permanent/gxaf/0\n";     
        myfile << "cmo/addatt/motet/Zv/VDOUBLE/scalar/nnodes/linear/permanent/gxaf/0\n";
        myfile << "cmo/copyatt/ motet mopts/ Pressure Pressure\n";
        myfile << "cmo/copyatt/ motet mopts/ Saturation Saturation\n";
        myfile << "cmo/copyatt/ motet mopts/ Porosity Porosity\n";
        myfile << "cmo/copyatt/ motet mopts/ Xv Xv\n";
        myfile << "cmo/copyatt/ motet mopts/ Yv Yv\n";
        myfile << "cmo/copyatt/ motet mopts/ Zv Zv\n";
        myfile << "cmo/modatt/motet/icr1/ioflag/l\n";
        myfile << "cmo/modatt/motet/isn1/ioflag/l\n";
	myfile << "#\n";
        myfile << "# Output tet mesh file\n";
        myfile << "dump / avs / OUT_AVS / motet\n"; 
        myfile << "dump / gmv / OUT_GMV / motet\n";
        myfile << "# Element adjacency information.\n";
        myfile << "dump / elem_adj_elem / OUT_ADJ / motet\n";
        myfile << "# Output FEHM files:\n";
        myfile << "# *.fehmn\n";
        myfile << "# *_material.zone\n"; 
        myfile << "# *_outside.zone (top, bottom, left, right, front, back)\n";
        myfile << "# *_outside_vor.area (you may not need this)\n";
        myfile << "# *.stor (area, volume coefficients move to zero for accuracy)\n";
        myfile << "dump / coord / OUT_FEHM / motet\n";
        myfile << "dump / zone_imt / OUT_FEHM / motet\n";
        myfile << "dump / zone_outside / OUT_FEHM / motet\n";
        myfile << "cmo/copy/mozero/motet\n";
        myfile << "cmo/select/mozero\n";
        myfile << "trans/1,0,0/zero\n";
        myfile << "dump / stor / OUT_FEHM / mozero\n";
        myfile << "cmo/printatt/motet/-all- /minmax\n";
	myfile << "\n";
	myfile << "finish\n";
	
	myfile.close();
}

//writeLgi_connect writes an lgi control file for LaGriT to connect grid points into a tet mesh
//This controls the input and output file names and controls what operations LaGritT should perform.
void writeLgi_connect(dataStruct data){
	ofstream myfile;
	string filename = data.filename;
	filename += ".lgi";
	myfile.open (filename.c_str());
	
	//write out standard LaGriT lgi file
	myfile << "# Define input and output file names\n";
	myfile << "#\n";
	myfile << "define / INPUT /" <<  "walkabout.h5" << ".inp\n";
	myfile << "define / OUT_AVS / " <<  "walkabout_tet.h5" << ".inp\n";
	myfile << "define / OUT_GMV / " <<  "walkabout.h5" << ".gmv\n";
	myfile << "define / OUT_FEHM / " <<  "walkabout.h5" << "\n";
	myfile << "define / OUT_ADJ / " <<  "walkabout.h5" << ".graph\n";
	myfile << "# Read in point set\n";
	myfile << "read / avs / INPUT / mo\n";
	myfile << "# Set some node attributes to default values\n";
	myfile << "cmo / setatt / mo / imt / 1 0 0 / 1\n";
	myfile << "cmo / setatt / mo / itp / 1 0 0 / 0\n";
	myfile << "#\n";
	myfile << "# Compute Delaunay tet connectivity of point set\n";
	myfile << "connect\n";
	myfile << "cmo / setatt / mo / itetclr / 1 0 0 / 1\n";
	myfile << "resetpts / itp\n";
	myfile << "# Output AVS file\n";
	myfile << "dump / avs / OUT_AVS / mo\n";
	myfile << "# Output GMV file\n";
	myfile << "dump / gmv / OUT_GMV / mo\n";
	myfile << "# Output FEHM files:\n";
	myfile << "# *.fehmn\n";
	myfile << "# *_material.zone (you may not need this)\n";
	myfile << "# *_interface.zone (you don't need this\n";
	myfile << "# *_multi_mat.zone (you don't need this)\n";
	myfile << "# *_outside.zone (top, bottom, left, right, front, back)\n";
	myfile << "# *_outside_vor.area (you may not need this)\n";
	myfile << "# *.stor (area, volume coefficients)\n";
	myfile << "dump / fehm / OUT_FEHM / mo\n";
	myfile << "#\n";
	myfile << "# Element adjacency information.\n";
	myfile << "dump / elem_adj_elem / OUT_ADJ / mo\n";
	myfile << "\n";
	myfile << "finish\n";
	
	myfile.close();
}


//writeControl writes a Walkabout control file (control.dat) which is used to specify parameters for 
//Walkabout
void writeControl(dataStruct data){
	ifstream ifile("control.dat");
	
	//check to see if a control file exists or not
	if(!ifile){
		//if not, write out data to the file
		ifile.close();
		ofstream myfile;
		string filename("control.dat");
		myfile.open (filename.c_str());
	
		cout << "Writing default Walkabout control.dat file.  Please modify control.dat to suit your needs.\n";
	
		myfile << "tmax 100\n";
		myfile << "seed 7127\n";
		myfile << "dxtarget 0.1\n";
		myfile << "dttarget 0.1\n";
		myfile << "maxstretch 1.3\n";
		myfile << "maxsteps 50000\n";
		myfile << "dt0 0.01\n";
		myfile << "toutfreq 2\n";
		myfile << "INITIAL\n";
		myfile << " RANDOM\n";
		myfile << " 100\n";
		myfile << data.minx << " " << data.miny << " " << data.minz << "\n";
		myfile << data.maxx << " " << data.maxy << " " << data.maxz << "\n";
		myfile << "\n";
		myfile << "DTENSOR\n";
		myfile << "1 0 0\n";
		myfile << "BF\n";
		myfile << " 1. 0.1 0.01 0.0\n";  
		myfile << "END\n   \n"; 
	
		myfile.close();
	}else{
		//otherwise, don't do anything
		ifile.close();
		cout << "Not writing Walkabout control.dat because it already exists.\n";
	}
}

//writeWalkaboutFile writes out a file (walkabout.files) which specifies what other files Walkabout should
//use as input and output files.
void writeWalkaboutFile(dataStruct data){
	ifstream ifile("walkabout.files");
	ofstream myfile;
	
	string filename = "walkabout.files";
	myfile.open (filename.c_str());
	
	//write out walkabout file data
	myfile << "fehmn:" << "walkabout.h5" << ".fehmn\n";
	myfile << "stor:" << "walkabout.h5" << ".stor\n";
	myfile << "ealist:" << "walkabout.h5" << ".graph\n";
	myfile << "ama:" << "walkabout.h5" << ".ama\n";
	myfile << "avs:" << "walkabout.h5" << ".avs\n"; 
	myfile << "control:control.dat\n"; 
	myfile << "trajout:traj.out\n";
	
	//close output file
	myfile.close();
}
//writePlumeCalcFiles writes out a set of files for running PlumeCalc
void writePlumeCalcFiles(dataStruct data){
	ifstream ifile("plumecalc.files");
	ofstream myfile;
	
	string filename = "plumecalc.files";
	myfile.open (filename.c_str());
	
	//write out plumecalc file data
	myfile << "walkabout.h5" << ".fehmn\n";  // grid
	myfile << "walkabout.h5" << ".stor\n";   // stor file
	myfile << "walkabout" << ".sptr2\n";  //sptr2 particle tracking file 
	myfile << "plumecalc" << ".rock\n";   //rock file
	myfile << "plumecalc" << ".sim\n";    //simulation control file
	myfile << "plumecalc" << ".out\n";    //output file
	
	//close output file
	myfile.close();

        //grid, stor, and sptr2 files are output from codes, while rock, sim files
        //have to be created here
        filename = "plumecalc.rock";
        myfile.open (filename.c_str());

	myfile << "rock\n";  
        //loop through points
        for(int a = 0; a < data.size; a++){
 //             write out density and porosity. 
 //             for now we use 1870 kg/m^3. Will modify this when the bulk density is available in the Amanzi output.
 //             myfile << a + 1 << "\t" << a + 1 << "\t" << "1" << "\t" << data.den[a] <<  "\t" << "0.00" << "\t"  <<  data.por[a] << "\n";
                 myfile << a + 1 << "\t" << a + 1 << "\t" << "1" << "\t" << 1875.00 <<  "\t" << "0.00" << "\t"  <<  data.por[a] << "\n";

        }
	myfile << "\n";  
	//close output file
	myfile.close();

/*   *.sim file will be prepared separately
        string filename = data.filename;
        filename += ".sim";
        myfile.open (filename.c_str());
        
	//close output file
	myfile.close();

*/

}

//writeObj writes an obj file, which is a mesh of tets that can be displayed in Paraview and other mesh
//rendering software.  The format is a standard Wavefront OBJ file.
void writeObj(dataStruct data){
	ofstream myfile;
	string filename = data.filename;
	
	double* verts;
	int* faces;
	int vCount, fCount;
	
	//check to see if reading the inp file was sucessful, and get the data if it was
	if(!readTets(data, &verts, &faces, vCount, fCount)){
		cout << "INP file not found.  OBJ file not created.\n";
		return;
	}
	
	//open output obj file
	filename += ".obj";
	myfile.open (filename.c_str());
	
	myfile << "# OBJ file created by h5output.py\n";
	myfile << "#\n";
	myfile << "g Object001\n\n";

	//write out vertices 
	for(int a = 0; a < vCount; a++){
		myfile << "v " << verts[a * 3 + 0] << " " << verts[a * 3 + 1] << " " << verts[a * 3 + 2] << "\n";
	}
	myfile << "\n";
	
	//write out faces
	for(int a = 0; a < fCount; a++){
		myfile << "f " << faces[a * 4 + 0] << " " << faces[a * 4 + 1] << " " << faces[a * 4 + 2] << " " << faces[a * 4 + 3] << "\n";
	}
	
	//close the obj file
	myfile.close();
}

//writePosVel is an unfinished and unused function that is meant and another postprocessing and verification
//tool.  It's probably not necessary at this point, as the same information can be found in 
//the OBJ and VTK files
void writePosVel(dataStruct data){
	ofstream myfile;
	string filename = data.filename;
	filename += ".posvel";
	myfile.open (filename.c_str());
	
	
	myfile.close();
}

//runLagrit runs LaGriT from within the application if the user specifies it
void runLagrit(dataStruct data, parameterStruct param){
	cout << "Running LaGriT with " << data.size << " nodes ...\n";
	//get the base filename
	string filename(data.filename);
	//get the location for lagrit
	string commandLine(param.lagrit);

        //check for input files
        ifstream ckfile("tet5.inp");
        if (!ckfile) {
                cout << "tet5.inp not found. Mesh and FEHM files not created.\n";
                return;
        }
        ckfile.close();

	//build the command with logging output
	commandLine += " < " + filename + ".lgi >lagrit.log";
	//run lagrit
	system(commandLine.c_str());
}

//runWalkabout runs Walkabout from within the application if the user specifies it
void runWalkabout(dataStruct data, parameterStruct param){
	cout << "Running Walkabout...\n";
	//get the walkabout path
	string commandLine(param.walkabout);
	//append logging control
	commandLine += ">walkabout.log";
	//run walkabout
	system(commandLine.c_str());
}

//runAmanzi is an unfinished and unused function.  The intent is to run Amanzi, but the functionality
//has not been implemented yet.
void runAmanzi(dataStruct data, parameterStruct param){

}

//readConfig is a function to read in information from config.ini, which is used to specify locations of
// This file is required for lagrit and or walkabout and provides executable locations 
// Sample file content:
//   lagrit:/n/local_linux/lagrit
//   walkabout: ~/bin/walkabout

parameterStruct readConfig(){
	ifstream ifile("config.ini");

        //check to see if the file was actually opened
        if(!ifile.is_open()){
                cout << "Early Exit Error: config.ini not found.\n";
                 exit(-1); 
        }

	parameterStruct param;
	string line;
	//loop through each line in the file, and check if they contain application location info
	while(getline(ifile, line)){
		if(line.find("lagrit") != string::npos){
			int pos = line.find(":") + 1;
			param.lagrit = line.substr(pos);
		}
		if(line.find("walkabout") != string::npos){
			int pos = line.find(":") + 1;
			param.walkabout = line.substr(pos);
		}
		if(line.find("amanzi") != string::npos){
			int pos = line.find(":") + 1;
			param.amanzi = line.substr(pos);
		}
	}
	
	ifile.close();
	return param;
}

//readTets is a helper function that reads in the tets and the vertices that make them up from a LaGriT
//input/output file (inp extension)
bool readTets(dataStruct data, double** verts, int** faces, int &vCount, int& fCount){
	ifstream ifile;
	int dummy;
	string buf;
	string filename(data.filename);
	filename += ".inp";
	vCount = fCount = 0;
	ifile.open(filename.c_str());
	
	//check to see if the file was actually opened
	if(!ifile.is_open()){	
		return false;
	}
	//get the number of vertices and faces
	ifile >> vCount >> fCount >> dummy >> dummy >> dummy;
	
	//allocate memory for vert and face arrays
	*verts = new double[vCount * 3];
	*faces = new int[fCount * 4];
	
	//read in verts
	for(int a = 0; a < vCount; a++){
		ifile >> dummy >> (*verts)[a*3+0] >> (*verts)[a*3+1] >> (*verts)[a*3+2];
	}
	//read in faces
	for(int a = 0; a < fCount; a++){
		ifile >> buf >> buf >> buf >> (*faces)[a*4+0] >> (*faces)[a*4+1] >> (*faces)[a*4+2] >> (*faces)[a*4+3];
	}
	//close input file
	ifile.close();
	return true;
}

//setOperations is a function that reads in command line parameters to determine what operations the 
//user wants the application to perform.
// Syntax: particleDriver walkabout_h5_file [pre runLagrit runWalkabout post]

void setOperations(parameterStruct &param, int argc, char* argv[]){
	//assume by default that no operations will be performed
	param.pre = false;
	param.runLagrit = false;
	param.runWalkabout = false;
	param.post = false;
	//check if the program was provided an input file
	if(argc < 2){
		cout << "No h5 file specified!\n";
		exit(-1);
	}
	//then check the remaining parameters for control statements
	for(int a = 2; a < argc; a++){
		if(strcmp(argv[a],"pre") == 0){
			param.pre = true;
		}
		else if(strcmp(argv[a],"runLagrit") == 0){
			param.runLagrit = true;
		}
		else if(strcmp(argv[a],"runWalkabout") == 0){
			param.runWalkabout = true;
		}
		else if(strcmp(argv[a],"post") == 0){
			param.post = true;
		}
	}
	
}

//writeVTKs writes 3 vtk files as a post processing step after Walkabout has been run.  It creates
//one VTK for the path lines that particles take, 1 for particle origin points, and 1 for particle 
//destination points.  These files can be viewed using Paraview, and other VTK compatible viz tools.
void writeVTKs(dataStruct data){
	string buf;
	ifstream ifile;
	string filename("traj.out");
	ifile.open(filename.c_str());
	
	//verify that the file was opened
	if(!ifile.is_open()){
		cout << "traj.out not found.  VTK files not created.\n";
		return;
	}
	
	//skip the first two lines
	getline(ifile, buf);
	getline(ifile, buf);
	int particleCount;
	int tsCount;
	int pointCount = 0;
	
	//get the particle count
	ifile >> particleCount;
	
	//particles effectively becomes a 3D array of doubles. 
	//		dimension 1 is the particle,
	//		dimension 2 is the position for that particle at that timestep ([time, x, y, z]), and
	//		dimension 3 is the step in the path for the particle

	vector < vector < vector <double> > > particles;
	
	//populate the particle trajectories structure
	for(int a = 0; a < particleCount; a++){
		//get the number of timesteps for the particle
		ifile >> tsCount;
		//create vectors for per particle data
		vector <double> t,x,y,z;
		vector< vector <double> > particle;
		//populate per particle data
		for(int b = 0; b < tsCount; b++){
			double tt, tx, ty, tz;
			ifile >> tt >> tx >> ty >> tz;
			t.push_back(tt);
			x.push_back(tx);
			y.push_back(ty);
			z.push_back(tz);
			pointCount++;
		}
		//add per particle data to particles vector
		particle.push_back(t);
		particle.push_back(x);
		particle.push_back(y);
		particle.push_back(z);
		particles.push_back(particle);
	}
	//open files for writing
	ofstream pathlines;
	filename = data.filename + string(".vtk");
	pathlines.open(filename.c_str());
	ofstream pathorig;
	filename = data.filename + string(".orig.vtk");
	pathorig.open(filename.c_str());
	ofstream pathdest;
	filename = data.filename + string(".dest.vtk");
	pathdest.open(filename.c_str());
	
	//write full path vtk file

	pathlines << "# vtk DataFile Version 3.0\n";
	pathlines << "Particle advection lines\n";
	pathlines << "ASCII\n";
	pathlines << "DATASET POLYDATA\n";
	pathlines << "POINTS " << pointCount << " float\n";
	for(int a = 0; a < particles.size(); a++){
		for(int b = 0; b < particles[a][0].size(); b++){
			pathlines << particles[a][1][b] << " " << particles[a][2][b] << " " << particles[a][3][b] << "\n";
		}
	}
	pathlines << "LINES " << particles.size() << " " << (pointCount + particles.size()) << "\n";
	int tempPointCount = 0;
	for(int a = 0; a < particles.size(); a++){
		pathlines << particles[a][0].size() << " ";
		for(int b = 0; b < particles[a][0].size(); b++){
			pathlines << tempPointCount << " ";
			tempPointCount++;
		}
		pathlines << "\n";
	}
	
	//write origin vtk file
	
	pathorig << "# vtk DataFile Version 3.0\n";
	pathorig << "Particle advection origins\n";
	pathorig << "ASCII\n";
	pathorig << "DATASET POLYDATA\n";
	pathorig << "POINTS " << particles.size() << " float\n";
	for(int a = 0; a < particles.size(); a++){
		pathorig << particles[a][1][0] << " " << particles[a][2][0] << " " << particles[a][3][0] << "\n";
	}
	
	//write destination vtk file
	
	pathdest << "# vtk DataFile Version 3.0\n";
	pathdest << "Particle advection destinations\n";
	pathdest << "ASCII\n";
	pathdest << "DATASET POLYDATA\n";
	pathdest << "POINTS " << particles.size() << " float\n";
	for(int a = 0; a < particles.size(); a++){
		pathdest << particles[a][1][(particles[a][1].size() - 1)] << " " << particles[a][2][(particles[a][2].size() - 1)] << " " << particles[a][3][(particles[a][3].size() - 1)] << "\n";
	}
	
	//close output and input files
	pathlines.close();
	pathorig.close();
	pathdest.close();
	ifile.close();
}

