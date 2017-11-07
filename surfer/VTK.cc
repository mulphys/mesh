/*
	Output to VTK
	Author: andrei.v.smirnov@gmail.com
*/
#include <iostream>
#include <fstream>

/*
Examples:
http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/TriangleArea
http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataExtractNormalshttp://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolygonalSurfaceContourLineInterpolator
http://www.vtk.org/Wiki/VTK/Examples/Cxx/IO/ReadPlainTextTriangles
http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataGetPoint

*/
#include <vtkNew.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkLine.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkSTLReader.h>
#include <vtkSTLWriter.h>
#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkTriangle.h>
#include <vtkDataSetMapper.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkTriangleFilter.h>

#define INT int
#define REAL double

extern "C" {// FORTRAN WRAPPERS
	void write_stl_(char*,INT*,INT*,REAL*,INT*,INT);
	void write_vtp_(char*,INT*,REAL*,REAL*,INT);
	void write_vector_vtp_(char*,INT*,INT*,REAL*,REAL*,INT);
void output_boundary_ (
	char *name, 
	INT *points_num, 
	INT *fc_num, 
	REAL *points, 
	INT *fc, // fc array from dtris3
	INT *vm, // vm array from dtris3
	INT len
);
	void surfmesh_(char*,INT*,REAL*,INT);
} // END FORTRAN WRAPPERS

int read_stl(char*,INT*,INT*,REAL*,REAL*,INT*);
void write_stl(char*,INT,INT,REAL*,INT*);
void write_ply(char*,INT,INT,INT,REAL*,INT*);
void write_vtu(char*,INT,INT,INT,REAL*,INT*);
void write_vtp(char*,INT,REAL*,REAL*);
void write_vector_vtp(const char*,INT,INT,REAL*,REAL*);
void write_lines_vtp ( 
	char *filename,
	INT np, // number of points
	REAL *P, // line edges array: P[2*np*DIM]
	REAL *Q // scalar variable defined on lines
);
void output_boundary (
	char *name, 
	INT points_num, 
	INT fc_num, 
	REAL *points, 
	INT *fc, // fc array from dtris3
	INT *vm  // vm array from dtris3
);
void surfmesh(char*,INT,REAL*);

void show_stl ( 
	char *filename,
	INT nnodes,
	INT ntris,
	double *nodes,
//-	double *normals,
	double *triangles
) {
	if ( filename == NULL || strlen(filename) == 0 ) {
		cerr << "Required parameters: Filename" << endl; cerr.flush();
		exit( EXIT_FAILURE );
	}

#ifdef DISPLAY
	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(filename);
	reader->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(reader->GetOutputPort());

	// Visualize
	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->SetBackground(.3, .6, .3); // Background color green

	renderWindow->Render();
	renderWindowInteractor->Start();
#endif
//	exit( EXIT_SUCCESS );
}

//http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataExtractNormals
bool GetCellNormals(vtkPolyData* polydata)
{
	std::cout << "Looking for cell normals..." << std::endl; cout.flush();
 
	// Count points
	vtkIdType numCells = polydata->GetNumberOfCells();
	std::cout << "There are " << numCells << " cells." << std::endl;
 
	// Count triangles
	vtkIdType numPolys = polydata->GetNumberOfPolys();
	std::cout << "There are " << numPolys << " polys." << std::endl;
 
	////////////////////////////////////////////////////////////////
	// Double normals in an array
	vtkDoubleArray* normalDataDouble =
		vtkDoubleArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));
 
	if(normalDataDouble)
		{
		int nc = normalDataDouble->GetNumberOfTuples();
		std::cout << "There are " << nc
						<< " components in normalDataDouble" << std::endl;
		return true;
		}
 
	////////////////////////////////////////////////////////////////
	// Double normals in an array
	vtkFloatArray* normalDataFloat =
		vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));
 
	if(normalDataFloat)
		{
		int nc = normalDataFloat->GetNumberOfTuples();
		std::cout << "There are " << nc
						<< " components in normalDataFloat" << std::endl;
		return true;
		}
 
	////////////////////////////////////////////////////////////////
	// Point normals
	vtkDoubleArray* normalsDouble =
		vtkDoubleArray::SafeDownCast(polydata->GetCellData()->GetNormals());
 
	if(normalsDouble)
		{
		std::cout << "There are " << normalsDouble->GetNumberOfComponents()
							<< " components in normalsDouble" << std::endl;
		return true;
		}
 
	////////////////////////////////////////////////////////////////
	// Point normals
	vtkFloatArray* normalsFloat =
		vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetNormals());
 
	if(normalsFloat)
		{
		std::cout << "There are " << normalsFloat->GetNumberOfComponents()
							<< " components in normalsFloat" << std::endl;
		return true;
		}
 
	/////////////////////////////////////////////////////////////////////
	// Generic type point normals
	vtkDataArray* normalsGeneric = polydata->GetCellData()->GetNormals(); //works
	if(normalsGeneric)
		{
		std::cout << "There are " << normalsGeneric->GetNumberOfTuples()
							<< " normals in normalsGeneric" << std::endl;
 
		double testDouble[3];
		normalsGeneric->GetTuple(0, testDouble);
 
		std::cout << "Double: " << testDouble[0] << " "
							<< testDouble[1] << " " << testDouble[2] << std::endl;
 
		// Can't do this:
		/*
		float testFloat[3];
		normalsGeneric->GetTuple(0, testFloat);
 
		std::cout << "Float: " << testFloat[0] << " "
							<< testFloat[1] << " " << testFloat[2] << std::endl;
		*/
		return true;
		}
 
 
	// If the function has not yet quit, there were none of these types of normals
	std::cout << "Normals not found!" << std::endl;
	return false;
 
}
vtkFloatArray* CellNormalsFloat(vtkPolyData* polydata)
{
	std::cout << "Looking for cell normals..." << std::endl;
 
	// Count points
	vtkIdType numCells = polydata->GetNumberOfCells();
	std::cout << "There are " << numCells << " cells." << std::endl;
 
	// Count triangles
	vtkIdType numPolys = polydata->GetNumberOfPolys();
	std::cout << "There are " << numPolys << " polys." << std::endl;
 
	////////////////////////////////////////////////////////////////
	// Float normals in an array
	vtkFloatArray* normalDataFloat =
		vtkFloatArray::SafeDownCast(polydata->GetCellData()->GetArray("Normals"));
 
	if(normalDataFloat)
	{
		int nc = normalDataFloat->GetNumberOfTuples();
		std::cout << "There are " << nc
						<< " components in normalDataFloat" << std::endl;
		return normalDataFloat;
	}
 	return NULL;
}
 

int read_stl ( 
	char *filename,
	int *nnodes,
	int *ntris,
	double *nodes,
	double *normals,
	int *triangles
) {

	vtkSmartPointer<vtkSTLReader> reader =
		vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName(filename);
	reader->Update();

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkPolyData> polyData = reader->GetOutput();
//vtkSmartPointer<vtkPolyData>::New();
	int np = polyData->GetNumberOfPoints();
	if((nodes=(double*)malloc(3*np))==NULL) {
		fprintf(stderr,"ERROR: Can't allocate %d elements for points\n",3*np);
		return 1;
	};
	*nnodes = np;

	cout << "read_stl:Number of points: " << np << endl;//-
	vtkSmartPointer<vtkPoints> points	= polyData->GetPoints();
//http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataGetPoint
	for(vtkIdType i = 0; i < np; i++) {
		double p[3];
		points->GetPoint(i,p);
		//polyData->GetPoint(i,p);
		// This is identical to:
		// polyData->GetPoints()->GetPoint(i,p);
		std::cout << "Point " << i << " : (" << p[0] << " " << p[1] << " " << p[2] << ")" << std::endl;
		for(int j=0; j<3; j++) 
			nodes[3*i+j] = p[j];
	}
	
//http://www.vtk.org/Wiki/VTK/Examples/Cxx/IO/ReadPlainTextTriangles

//http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/TriangleArea

	int nc = polyData->GetNumberOfCells();
	cout << "read_stl:Ncells="<<nc<<endl;//-

	if((normals=(double*)malloc(3*nc))==NULL) {
		fprintf(stderr,"ERROR: Can't allocate %d elements for normals\n",3*nc);
		return 1;
	};
	if((triangles=(int*)malloc(3*nc))==NULL) {
		fprintf(stderr,"ERROR: Can't allocate %d elements for normals\n",9*nc);
		return 1;
	};
	*ntris = nc;

	for(vtkIdType i = 0; i < nc; i++) {
		vtkCell* cell = polyData->GetCell(i);
		vtkTriangle* triangle = dynamic_cast<vtkTriangle*>(cell);
		double p0[3];
		double p1[3];
		double p2[3];
		triangle->GetPoints()->GetPoint(0, p0);
		std::cout << "p0: " << p0[0] << " " << p0[1] << " " << p0[2] << std::endl;
		triangle->GetPoints()->GetPoint(1, p1);
		std::cout << "p1: " << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
		triangle->GetPoints()->GetPoint(2, p2);
		std::cout << "p2: " << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
 
		double area = vtkTriangle::TriangleArea(p0, p1, p2);
 
		std::cout << "area of triangle " << i << ": " << area << std::endl; 

		return 0;
	}

//http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/PolyDataExtractNormals

	vtkFloatArray* cellNormals = NULL;
	bool hasCellNormals = GetCellNormals(polyData);
	if (hasCellNormals) 
		cout << "Found cell normals\n";
	else {
		std::cout << "No cell normals were found. Computing normals..." << std::endl;
 
		// Generate normals
		vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
#if VTK_MAJOR_VERSION <= 5
		normalGenerator->SetInput(polyData);
#else
		normalGenerator->SetInputData(polyData);
#endif
		normalGenerator->ComputePointNormalsOff();
		normalGenerator->ComputeCellNormalsOn();
		normalGenerator->Update();
		/*
		// Optional settings
		normalGenerator->SetFeatureAngle(0.1);
		normalGenerator->SetSplitting(1);
		normalGenerator->SetConsistency(0);
		normalGenerator->SetAutoOrientNormals(0);
		normalGenerator->SetComputePointNormals(1);
		normalGenerator->SetComputeCellNormals(0);
		normalGenerator->SetFlipNormals(0);
		normalGenerator->SetNonManifoldTraversal(1);
		*/
 
		polyData = normalGenerator->GetOutput();
 
		// Try to read normals again
		hasCellNormals = GetCellNormals(polyData);
 
		std::cout << "On the second try, has cell normals? " << hasCellNormals << std::endl;

	}

	if (hasCellNormals) {
		cellNormals = CellNormalsFloat(polyData);
//		cout << "Number of face-normals: "<<cellNormals->size()<<endl;
		cout << "Number of face-normals: "<<cellNormals->GetNumberOfTuples()<<endl;
	}
	// Returns:
	*ntris = np;
}

void write_stl ( 
	char *filename, 
	INT np, // number of points
	INT nt, // number of triangles
	REAL *P, // points
//-	double *N, // triangle-normals
	INT *T // indexes to triangles
) { 
	//http://www.vtk.org/Wiki/VTK/Examples/Cxx/IO/ReadPlainTextTriangles
	if ( filename == NULL || strlen(filename) == 0 ) {
		cerr << "Required parameters: Filename" << endl; cerr.flush();
		exit( EXIT_FAILURE );
	}
	cout << "Writing "<<np<<" nodes and "<<nt<<" triangles to " << filename << endl; cout.flush();
	vtkIdType 
		number_of_points = np,
		number_of_triangles = nt;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(number_of_points);
	for (vtkIdType i = 0; i < number_of_points; i++) {
		int j=3*i;
		points->SetPoint(i, P[j], P[j+1], P[j+2]);
	}

	vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
	for (vtkIdType i = 0; i < number_of_triangles; i++)
	{	int j = 3*i;
		vtkIdType 
			a = T[j], 
			b = T[j+1], 
			c = T[j+2];

		polys->InsertNextCell(3);
		polys->InsertCellPoint(a);
		polys->InsertCellPoint(b);
		polys->InsertCellPoint(c);
	}
	vtkPolyData * polydata = vtkPolyData::New();
	polydata->SetPoints(points);
	polydata->SetPolys(polys);
//-	return polydata;

//-	vtkSmartPointer<vtkSphereSource> sphereSource =
//-		vtkSmartPointer<vtkSphereSource>::New();
//-	sphereSource->Update();

	vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
	stlWriter->SetFileName(filename);
//-	stlWriter->SetInputConnection(sphereSource->GetOutputPort());
#if VTK_MAJOR_VERSION <= 5
  stlWriter->SetInput(polydata);
#else
  stlWriter->SetInputData(polydata);
#endif
	stlWriter->Write();

	// Read and display for verification
//-	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
//-	reader->SetFileName(filename);
//-	reader->Update();

#ifdef DISPLAY
	
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
//-	mapper->SetInputConnection(sphereSource->GetOutputPort());
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(polydata);
#else
	mapper->SetInputData(polydata);
#endif

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->SetBackground(.3, .6, .3); // Background color green

	renderWindow->Render();
	renderWindowInteractor->Start();
#endif
//	exit( EXIT_SUCCESS );
}

void write_stl_ ( // FORTRAN WRAPPER 
	char *filename, 
	INT *np, // number of points
	INT *nt, // number of triangles
	REAL *P, // points
//-	double *N, // triangle-normals
	INT *T, // indexes to triangles
	INT len
) { 
	char name[510];
	strncpy(name, filename, len);
	name[len]='\0';
	write_stl(name, *np, *nt, P, T);
}

void write_ply ( 
	char *filename,
	INT nv, // number of vertices in a cell
	INT nn, // number of nodes in a mesh
	INT nc, // number of cells in a mesh
	REAL *nodes,
	INT *cells
) {
	if ( filename == NULL || strlen(filename) == 0 ) {
		cerr << "Required parameters: Filename" << endl; cerr.flush();
		exit( EXIT_FAILURE );
	}
	cout<<"write_ply: filename="<<filename<<", nv="<<nv<<", nn="<<nn<<", nc="<<nc<<endl; cout.flush();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(nn);
	for (INT i = 0; i < nn; i++) {
		INT j=3*i;
		points->SetPoint((vtkIdType)i, nodes[j], nodes[j+1], nodes[j+2]);
	}

	vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();
	for (vtkIdType i = 0; i < nc; i++) {
		polys->InsertNextCell(nv);
		for(INT j = 0; j<nv; j++) {
			INT k = nv*i+j;
			polys->InsertCellPoint((vtkIdType)cells[k]);
		}
	}
	vtkPolyData * polydata = vtkPolyData::New();
	polydata->SetPoints(points);
	polydata->SetPolys(polys);

	vtkSmartPointer<vtkPLYWriter> plyWriter =
		vtkSmartPointer<vtkPLYWriter>::New();
	plyWriter->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
  plyWriter->SetInput(polydata);
#else
  plyWriter->SetInputData(polydata);
#endif
	cout<<"write_ply: writing to "<<filename<<endl; cout.flush();
	plyWriter->Write();

	// Read and display for verification
//-	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
//-	reader->SetFileName(filename);
//-	reader->Update();

#ifdef DISPLAY
	
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInput(polydata);
#else
	mapper->SetInputData(polydata);
#endif

//+	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
//+	actor->SetMapper(mapper);
//+
//+	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
//+	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
//+	renderWindow->AddRenderer(renderer);
//+	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
//+	renderWindowInteractor->SetRenderWindow(renderWindow);
//+
//+	renderer->AddActor(actor);
//+	renderer->SetBackground(.3, .6, .3); // Background color green
//+
//+	renderWindow->Render();
//+	renderWindowInteractor->Start();
//+	
//+	exit( EXIT_SUCCESS );
//+
//+
//-	vtkSmartPointer<vtkSphereSource> sphereSource =
//-		vtkSmartPointer<vtkSphereSource>::New();
//-	sphereSource->Update();
//-	
//-	vtkSmartPointer<vtkPLYWriter> plyWriter =
//-		vtkSmartPointer<vtkPLYWriter>::New();
//-	plyWriter->SetFileName(filename);
//-	plyWriter->SetInputConnection(sphereSource->GetOutputPort());
//-	plyWriter->Write();
//-
//-	// Read and display for verification
//-	vtkSmartPointer<vtkPLYReader> reader =
//-		vtkSmartPointer<vtkPLYReader>::New();
//-	reader->SetFileName(filename);
//-	reader->Update();
//-	
//-	vtkSmartPointer<vtkPolyDataMapper> mapper =
//-		vtkSmartPointer<vtkPolyDataMapper>::New();
//-	mapper->SetInputConnection(reader->GetOutputPort());

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->SetBackground(.3, .6, .3); // Background color green

	renderWindow->Render();
	renderWindowInteractor->Start();
#endif
	
//	exit( EXIT_SUCCESS );
}


void write_vtu ( 
	char *filename,
	INT nv, // number of vertices in a cell
	INT nn, // number of nodes in a mesh
	INT nc, // number of cells in a mesh
	REAL *nodes,
	INT *cells
//	double *node_data, // data at node locations
//	double *cell_data // data at cell locations
) {
	if ( filename == NULL || strlen(filename) == 0 ) {
		cerr << "Required parameters: Filename" << endl; cerr.flush();
		exit( EXIT_FAILURE );
	}
	cout<<"write_vtu: filename="<<filename<<", nv="<<nv<<", nn="<<nn<<", nc="<<nc<<endl; cout.flush();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(nn);
	for (INT i = 0; i < nn; i++) {
		INT j=3*i;
		points->SetPoint((vtkIdType)i, nodes[j], nodes[j+1], nodes[j+2]);
	}

//-   vtkSmartPointer<vtkTetra> tetra =
//-    vtkSmartPointer<vtkTetra>::New();
//-
//-  tetra->GetPointIds()->SetId(0, 0);
//-  tetra->GetPointIds()->SetId(1, 1);
//-  tetra->GetPointIds()->SetId(2, 2);
//-  tetra->GetPointIds()->SetId(3, 3);

//-  vtkSmartPointer<vtkCellArray> cellArray =
//-    vtkSmartPointer<vtkCellArray>::New();
//-  cellArray->InsertNextCell(tetra);
//-

  vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
    vtkSmartPointer<vtkUnstructuredGrid>::New();

  unstructuredGrid->SetPoints(points);

	for (INT ic=0; ic<nc; ic++) {
		vtkIdType ptIds[nv];
		for (INT iv=0; iv<nv; iv++) {
			ptIds[iv] = cells[ic*nv+iv]-1;
		}
		unstructuredGrid->InsertNextCell( VTK_TETRA, nv, ptIds );
	}

  // Write file
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(unstructuredGrid);
#else
  writer->SetInputData(unstructuredGrid);
#endif
  writer->Write();

	// Read and display for verification

#ifdef DISPLAY	
	vtkSmartPointer<vtkDataSetMapper> mapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(unstructuredGrid->GetProducerPort());
#else
	mapper->SetInputData(unstructuredGrid);
#endif

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->SetBackground(.3, .6, .3); // Background color green

	renderWindow->Render();
	renderWindowInteractor->Start();
#endif
//	exit( EXIT_SUCCESS );
}
//
// Output Scalar in VTP format
//
void write_vtp ( 
	char *filename,
	INT np, // number of points
	REAL *P, // points array
	REAL *Q // scalar variable
//TODO: implement case if Q==NULL
) {
	if ( filename == NULL || strlen(filename) == 0 ) {
		cerr << "Required parameters: Filename" << endl; cerr.flush();
		exit( EXIT_FAILURE );
	}

	cout<<"write_vtp: filename="<<filename<<", np="<<np<<endl; cout.flush();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for ( INT i = 0; i < np; i++ ) {
		INT j = 3*i;
		points->InsertNextPoint ( P[j], P[j+1], P[j+2] );
	}

  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();

	for ( INT i = 0; i < np; i++ ) {
		INT j = 3*i;
    vtkIdType pid[1];

    //Add a point to the polydata and save its index, which we will use to create the vertex on that point.
		pid[0] = points->InsertNextPoint ( P[j], P[j+1], P[j+2] );

    //create a vertex cell on the point that was just added.
    vertices->InsertNextCell ( 1,pid );
	}

  // Create a polydata object and add the points to it.
//http://www.vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/ParallelFlowPaths/Testing/Cxx/TestPLagrangianParticleTracker.cxx

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
  polydata->SetVerts(vertices);

	vtkPointData *pointData = polydata->GetPointData();

	vtkSmartPointer<vtkDoubleArray> scalar = vtkSmartPointer<vtkDoubleArray>::New();
	scalar->SetNumberOfComponents(1);
	scalar->SetNumberOfTuples(np);
	scalar->SetName("Scalar");

	for (vtkIdType i=0; i<np; i++) {
		double q = Q[i];
		scalar->SetValue (i, q);

	}
	pointData->AddArray(scalar);

  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(polydata);
#else
  writer->SetInputData(polydata);
#endif

  // Optional - set the mode. The default is binary.
  //writer->SetDataModeToBinary();
  writer->SetDataModeToAscii();

  writer->Write();
}
//
// Output Lines in VTP format
//
void write_lines_vtp ( 
	const char *filename,
	INT nl, // number of points
	REAL *P, // line edges array: P[2*np*DIM]
	REAL *Q // scalar variable defined on lines
) {//TODO: did not work 
	if ( filename == NULL || strlen(filename) == 0 ) {
		cerr << "Required parameters: Filename" << endl; cerr.flush();
		exit( EXIT_FAILURE );
	}

	cout<<"write_vtp: filename="<<filename<<", nl="<<nl<<endl; cout.flush();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for ( INT i = 0; i < nl; i++ ) {
		for (int j=0; j<2; j++) {
			INT n = 3*(2*i + j);
			points->InsertNextPoint ( P[n], P[n+1], P[n+2] );
		}
	}

//?  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();
//?
//?	for ( INT i = 0; i < np; i++ ) {
//?		INT j = 3*i;
//?    vtkIdType pid[1];
//?
//?    //Add a point to the polydata and save its index, which we will use to create the vertex on that point.
//?		pid[0] = points->InsertNextPoint ( P[j], P[j+1], P[j+2] );
//?
//?    //create a vertex cell on the point that was just added.
//?    vertices->InsertNextCell ( 1,pid );
//?	}

  // Create a polydata object and add the points to it.
//http://www.vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/ParallelFlowPaths/Testing/Cxx/TestPLagrangianParticleTracker.cxx

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
//?  polydata->SetVerts(vertices);

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	for (int i=0; i<nl; i++) {
		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

		for (int j=0; j<2; j++) {
			line->GetPointIds()->SetId(i, 2*i+j);
		}
		lines->InsertNextCell(line);
	}

	polydata->SetLines(lines);

	vtkSmartPointer<vtkDoubleArray> scalar = vtkSmartPointer<vtkDoubleArray>::New();
	scalar->SetNumberOfComponents(1);
	scalar->SetNumberOfTuples(nl);
	scalar->SetName("Scalar");

	for (vtkIdType i=0; i<nl; i++) {
		double q = Q[i];
		scalar->SetValue (i, q);
	}

//-	vtkPointData *pointData = polydata->GetPointData();
//-	pointData->AddArray(scalar);
	polydata->GetCellData()->SetScalars(scalar);	

  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(polydata);
#else
  writer->SetInputData(polydata);
#endif

  // Optional - set the mode. The default is binary.
  //writer->SetDataModeToBinary();
  writer->SetDataModeToAscii();

  writer->Write();
}


// FORTRAN CALL WRAPPER
void write_vector_vtp_(
	char*filename,
	INT *np, 
	INT *nv,
	REAL *P,
	REAL *V,
	INT filename_length
) {
	char name[510];
	strncpy(name, filename, filename_length);
	name[filename_length]='\0';
	cout<<"Writing "<<*np<<" vectors to "<<name<<endl;
	cout.flush();
	write_vector_vtp((const char *)name, *np, *nv, P, V);
}

//
// Output Vector in VTP format
//
//https://lorensen.github.io/VTKExamples/site/Cxx/PolyData/PolyDataCellNormals/
void write_vector_vtp ( 
	const char *filename,
	INT np, // number of points
	INT nv, // number of components of a vector
	REAL *P, // points array
	REAL *V // vector variable
//TODO: implement case if Q==NULL
) {
	if ( filename == NULL || strlen(filename) == 0 ) {
		cerr << "Required parameters: Filename" << endl; cerr.flush();
		exit( EXIT_FAILURE );
	}

	cout<<"write_vtp: "<<np<<" data points to "<<filename<<endl; cout.flush();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkCellArray> vertices = vtkSmartPointer<vtkCellArray>::New();

	for ( INT i = 0; i < np; i++ ) {
		INT j = 3*i;
    vtkIdType pid[1];

    //Add a point to the polydata and save its index, which we will use to create the vertex on that point.
		pid[0] = points->InsertNextPoint ( P[j], P[j+1], P[j+2] );

    //create a vertex cell on the point that was just added.
    vertices->InsertNextCell ( 1,pid );
	}

  // Create a polydata object and add the points to it.
//http://www.vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/ParallelFlowPaths/Testing/Cxx/TestPLagrangianParticleTracker.cxx

//https://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/FieldData
//?  vtkSmartPointer<vtkDoubleArray> location = vtkSmartPointer<vtkDoubleArray>::New();

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

	polydata->SetPoints(points);
  polydata->SetVerts(vertices);

//?  location->SetNumberOfComponents(3);
//?  location->SetName("Coordinates");
//?
//?  // Create the data to store (here we just use (0,0,0))
//?	for ( INT i = 0; i < np; i++ ) {
//?  	double locationValue[3];
//?		for (int j=0; j<3; j++) {
//?			locationValue[j] = P[3*i+j];
//?		}
//?		location->InsertNextTuple (locationValue);
//?	} 
//?
//?  // The data is added to FIELD data (rather than POINT data as usual)
//?  polydata->GetFieldData()->AddArray(location);


//https://lorensen.github.io/VTKExamples/site/Cxx/PolyData/WarpVector/
//?	vtkPointData *pointData = polydata->GetPointData();


	vtkSmartPointer<vtkDoubleArray> vectors = vtkSmartPointer<vtkDoubleArray>::New();
	vectors->SetNumberOfComponents(3);
//?	vector->SetNumberOfTuples(np);
	vectors->SetNumberOfTuples(polydata->GetNumberOfPoints());
	vectors->SetName("Vector");

	for (vtkIdType ip=0; ip<np; ip++) {
//https://www.vtk.org/Wiki/VTK/Examples/Cxx/Utilities/VectorArrayUnknownLength
		double v[3];
		for(vtkIdType iv = 0; iv < nv; iv++) {
			v[iv] = V[ip*nv+iv];
		}
		vectors->SetTuple(ip, v);
	}
//?	pointData->AddArray(vector);
//?	polydata->GetFieldData()->AddArray(vector);
	polydata->GetCellData()->SetNormals(vectors);
//?	polydata->GetPointData()->SetNormals(vectors);

  // Write the file
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(filename);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(polydata);
#else
  writer->SetInputData(polydata);
#endif

  // Optional - set the mode. The default is binary.
  //writer->SetDataModeToBinary();
  writer->SetDataModeToAscii();

  writer->Write();

//?	delete tuple;
}

// FORTRAN CALL WRAPPER
void write_vtp_(
	char*filename,
	INT *np, 
	REAL *P,
	REAL *Q,
	INT filename_length
) {
	char name[510];
	strncpy(name, filename, filename_length);
	name[filename_length]='\0';
	cout<<"Writing "<<*np<<" points to "<<name<<endl;
	cout.flush();
	write_vtp(name, *np, P, Q);
}

void output_boundary (
	char *filename, 
	INT points_num, 
	INT fc_num, 
	REAL *points, 
	INT *fc, // fc array from dtris3
	INT *vm  // vm array from dtris3
) {
	INT 
		faces_num = 0,
		*faces=NULL;

	// Count faces:
	for (INT i=0; i< fc_num; i++) {
		if (fc[7*i+4] <= 0)
			faces_num++;
	}

	cout << "Found boundary of "<<faces_num<<" faces"<<endl; cout.flush();
	if (
		(faces = (INT*) calloc(3*faces_num, sizeof(*faces))) == NULL
	) {
		fprintf(stderr,"ERROR: Failed to allocate 3 x %d faces\n",faces_num);
		exit(1);
	}

	INT n = 0;	
	for (INT i=0; i<fc_num; i++) {
		INT k=7*i;
		if (fc[k+4] <= 0) {
			for (int j=0; j<3; j++) {
				faces[n++] = vm[fc[k+j]-1]-1;
			}
		}
	}

	write_stl (
		filename,
		points_num, 
		faces_num, 
		points,
		faces
	);
	free(faces);
}

void output_boundary_ ( //FORTRAN WRAPPER
	char *filename, 
	INT *points_num, 
	INT *fc_num, 
	REAL *points, 
	INT *fc, // fc array from dtris3
	INT *vm,  // vm array from dtris3
	INT len
) {

	char name[510];
	strncpy(name, filename, len);
	name[len]='\0';

		output_boundary (
			name, 
			*points_num, 
			*fc_num, 
			points, 
			fc,vm
		);

}

// 
// Runs Delaunay 3D
// Extracts surface mesh
//
void surfmesh_( //FORTRAN WRAPPER
	char*filename,
	INT *np, 
	REAL *P,
	INT filename_length
) {
	char name[510];
	strncpy(name, filename, filename_length);
	name[filename_length]='\0';
	cout<<"Writing "<<*np<<" points to "<<name<<endl;
	cout.flush();
	surfmesh(name, *np, P);
}

//
// Create Delaunay3D mesh
// and extract triangulated surface
//
void surfmesh ( 
	char *stl_filename,
	INT np, // number of vertices in a cell
	REAL *P // points array
) {
	if ( stl_filename == NULL || strlen(stl_filename) == 0 ) {
		cerr << "Required parameters: Filename" << endl; cerr.flush();
		exit( EXIT_FAILURE );
	}

	cout<<"surfmesh: filename="<<stl_filename<<", np="<<np<<endl; cout.flush();

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	for ( INT i = 0; i < np; i++ ) {
		INT j = 3*i;
		points->InsertNextPoint ( P[j], P[j+1], P[j+2] );
	}
  // Create a polydata object and add the points to it.
//http://www.vtk.org/gitweb?p=VTK.git;a=blob;f=Filters/ParallelFlowPaths/Testing/Cxx/TestPLagrangianParticleTracker.cxx

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);

	vtkPointData *pointData = polydata->GetPointData();


	bool ascii = false, silent = false;

	vtkSmartPointer<vtkDelaunay3D> mesh = vtkSmartPointer<vtkDelaunay3D>::New();
	mesh->SetInputData (polydata);

	if (!silent) cout << "Writing STL data to "<<stl_filename<<endl;

//https://public.kitware.com/pipermail/vtkusers/2009-October/054301.html
//https://www.vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/Medical/Cxx/GenerateCubesFromLabels.cxx

	mesh->AlphaLinesOn();
	mesh->AlphaTrisOn();
	mesh->AlphaTetsOn();
	mesh->AlphaVertsOn();
//	mesh->BoundingTriangulationOn();
	REAL 
		tolerance = mesh->GetTolerance(),
		alpha = mesh->GetAlpha(),
		offset = mesh->GetOffset();
	cout 
	<< "Tolerance: "<<tolerance<<endl//-
	<< "Alpha: "<<alpha<<endl//-
	<< "Offset: "<<offset<<endl;//-
//	tolerance = 0.0;
//	alpha = 0.0;
//	offset = 0.0;
	mesh->SetAlpha(alpha);
	mesh->SetTolerance(tolerance);
	mesh->SetOffset(tolerance);

	vtkGeometryFilter* geomFilter = vtkGeometryFilter::New();
	#if VTK_MAJOR_VERSION <= 5
	geomFilter->SetInput(mesh->GetOutput()); 
	#else
	geomFilter->SetInputConnection(mesh->GetOutputPort()); 
	#endif
	geomFilter->MergingOn();
	
	vtkTriangleFilter *triangles = vtkTriangleFilter::New();
	#if VTK_MAJOR_VERSION <= 5
	triangles->SetInput(geomFilter->GetOutput());
	#else
	triangles->SetInputConnection(geomFilter->GetOutputPort());
	#endif
	
	vtkSTLWriter *stl = vtkSTLWriter::New();
	stl->SetInputConnection(triangles->GetOutputPort());
	stl->SetFileName(stl_filename);
	if (ascii) stl->SetFileTypeToASCII();
	else stl->SetFileTypeToBinary();
	stl->Write();

}
