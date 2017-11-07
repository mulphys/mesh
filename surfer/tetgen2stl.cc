//
// Author: andrei.v.smirnov@gmail.com

#include <cstdlib>
#include <unistd.h>

#include <iostream>
#include <fstream>

#define INT int
#define REAL double


extern void write_stl(char*,INT,INT,REAL*,INT*);

using namespace std;

void doc() {
		cout 
		<< "Converts 3D Tetgen data to STL" <<endl
		<< "\t- Reads 3D points from file in .node format" <<endl
		<< "\t- Reads triangles data in .face format"         <<endl
		<< "\t- Outputs boundary in STL format"          <<endl;
}

void usage(char *prog) {
		cout << "Usage: " << prog 
		<< " [-hs] <input>.node <input>.face <output>.stl]"       <<endl 
		<< "Where:"                                        <<endl
		<< "\t<input>.node -input node coordinates in TETGEN format" <<endl
		<< "\t<input>.face -input face coordinates in TETGEN format" <<endl
		<< "\t<output>.stl - output file in STL format"    <<endl
		<< "\t-a - ASCII output"                           <<endl
		<< "\t-d - face data present in <input>.face file" <<endl
		<< "\t-s - silent mode: no output"                 <<endl;
}

int main ( int argc, char *argv[] )
{
	bool silent = false, facedata = false;
	int c;
	while ((c = getopt (argc, argv, "adhs")) != -1)
	switch(c)
	{
		break;
		case 'h':
		doc();
		usage(argv[0]);
		return 0;
		break;
		case 's':
		silent = true;
		break;
		case 'd':
		facedata = true;
		break;
	}

	int nfiles = argc - optind; 
	if (nfiles < 3)
	{	cerr << "ERROR: Wrong command line"<<endl;
		usage(argv[0]);
		return EXIT_FAILURE;
	}
	char 
		*node_filename = argv[optind],
		*face_filename = argv[optind+1],
		*output_filename = argv[optind+2];


	ifstream node_inp(node_filename, std::ifstream::in);
	if (!silent) cout <<"Reading file "<<node_filename<<endl;

	int nnodes, ndim, n;
	node_inp >> nnodes >> ndim >> n >> n;
	cout << "Number of nodes: "<<nnodes<<endl;
	REAL *X = new REAL[nnodes*ndim];
	for (int i=0; i< nnodes; i++) {
		node_inp >> n;
		for (int j=0; j<ndim; j++) {
			REAL x;
			node_inp >> x;
			X[(n-1)*ndim+j] = x;
		}
	}
	node_inp.close();


	ifstream face_inp(face_filename, std::ifstream::in);
	if (!silent) cout <<"Reading file "<<face_filename<<endl;

	int nfaces, nverts = 3;
	face_inp >> nfaces >> n;
	cout << "Number of faces: "<<nfaces<<endl;
	int *F = new int[nfaces*ndim];
	for (int i=0; i< nfaces; i++) {
		int data;
		face_inp >> n;
		for (int j=0; j<nverts; j++) {
			int k;
			face_inp >> k;
			F[(n-1)*ndim+j] = k-1;

		}
		if (facedata) face_inp >> data; 
	}
	face_inp.close();


	if (!silent) cout << "Writing STL data to "<<output_filename<<endl;
	write_stl(output_filename, nnodes, nfaces, X, F);

	delete F, X;

	return EXIT_SUCCESS;
}
