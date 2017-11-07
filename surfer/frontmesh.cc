/// \file 
/// \brief Implementation of the propagating front routines.
/// 
/// Author: Andrei Smirnov
/// andrei.v.smirnov@gmail.com
/// 

#include <cstdlib>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "frontmesh.h"

using namespace std;

bool 
	verbose = false,
	silent = false;

void doc() {
		cout 
		<< "Produces triangulated surface of point data" <<endl
		<< "\t- Reads 3D points and surface normals "    <<endl
		<< "\t- Performs boundary triangulation"         <<endl
		<< "\t- Outputs boundary "          <<endl;
}

void usage(char *prog) {
		cout << "Usage: " << prog 
		<< " [-s] <coordinates>.[node|dat] <normals>.[node|dat] [<output>.face]" <<endl 
		<< "Where:" <<endl
		<< "\t<coordinates>.node -input file of point coordinates in Tetgen node or Qhull dat format" <<endl
		<< "\t<normals>.norm -input file with surface normal vectors at points locations"    <<endl
		<< "\t<output>.face - output file in Tetgen face format" <<endl
		<< "\t-s - silent mode: no output" <<endl
		<< "\t-v - verbose mode: extra output" <<endl;
}


///
/// Main routine perfoming IO and calling the mesher.
///
int main ( int argc, char *argv[] )
{
	int frontmesh (int nPoints, REAL *Coords, REAL *Normals, int *Indices, int *&Faces); 
	silent = false;
	verbose = false;
	enum Format {dat_format, node_format} io_format;
	int c;
	while ((c = getopt (argc, argv, "hsv")) != -1)
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
		case 'v':
		verbose = true;
		break;
	}

	if (silent && verbose) {
		cerr << "CANNOT USE OPTIONS -s AND -v AT THE SAME TIME"<<endl;
		return 1;
	}

	int nfiles = argc - optind; 
	if (nfiles != 3)
	{	cerr << "ERROR: Wrong command line"<<endl;
		usage(argv[0]);
		return EXIT_FAILURE;
	}
	io_format = dat_format;
	char 
		*points_filename = argv[optind],
		*snorms_filename = argv[optind+1],
		*output_filename = argv[optind+2];

	if (
		   strstr(points_filename, ".dat")==nullptr 
		&& strstr(points_filename, ".node")==nullptr) {
		cerr << "File "<<points_filename<<" needs suffix '.node' or '.dat'" 
		<< endl;
		exit(1);
	}
	io_format = dat_format;
	if (strstr(points_filename, ".node")!=nullptr) io_format = node_format;

	if (strstr(output_filename, ".face")==nullptr) {
		cerr << "File "<<output_filename<<" needs suffix .face" << endl;
	}

	if (!silent) {
		cout <<"Reading "<<points_filename<<endl;
	}

	ifstream inp(points_filename);

	int nPoints, nDim, n;
	switch(io_format) {
		case dat_format:
		inp >> nDim >> nPoints; 
		break;
		case node_format:
		inp >> nPoints >> nDim >> n >> n;
		break;
	}
	if(!silent) {
	cout << "Number of Dimensions: "<< nDim   << endl
	     << "Number of nodes: "     << nPoints << endl;
	}

	REAL *Points = new REAL[nDim*nPoints];
	
	int 
		*Indices = new int[nPoints],
		*Faces = nullptr;

	for (int i=0; i<nPoints; i++) {
		int n = nDim*i, k;
		switch(io_format) {
			case dat_format:
			k = i+1;
			break;
			case node_format:
			inp >> k;
			break;
		}
		Indices[i] = k;
		for (int j=0; j<nDim; j++) {
			REAL x;
			inp >> x;
			Points[n+j] = x;
		}
	}

	inp.close();

	if (!silent) {
		cout <<"Reading "<<snorms_filename<<endl;
	}

	inp.open(snorms_filename);

	int nNorms, mDim;

	switch(io_format) {
		case dat_format:
		inp >> n >> nNorms; 
		break;
		case node_format:
		inp >> nNorms >> mDim >> n >> n;
		break;
	}
	if (mDim!=nDim) {
			cerr << "Dimensions dont agree: "<<n<<" != "<<nDim<<endl;
			exit(1);
	}
	if (nPoints!=nNorms) {
			cerr << "Number of points "<<nPoints<<" is not equal to number of surface normals: "<<nNorms<<endl;
			exit(1);
	}
	REAL 	*Snorms = new REAL[nDim*nPoints];
	
	for (int i=0; i<nPoints; i++) {
		int n = nDim*i, k;
		switch(io_format) {
			case dat_format:
			k = i+1;
			break;
			case node_format:
			inp >> k;
			break;
		}
		for (int j=0; j<nDim; j++) {
			REAL x;
			inp >> x;
			Snorms[n+j] = x;
		}
	}

	inp.close();

	int nFaces = frontmesh(nPoints, Points, Snorms, Indices, Faces);

	if (!silent) {
		cout <<"Writing "<<nFaces<<" faces to "<<output_filename<<endl;
	}
	ofstream out(output_filename);
	out << nFaces << '\t' << NV << endl;
	for (int i=0; i<nFaces; i++) {
		int n = i*NV;
		if (io_format == node_format) out << i+1;
		for (int j=0; j<NV; j++) {
			out << '\t' << Faces[n+j];
		}
		out << endl;
	}
	out.close();
	delete[] Points;
	delete[] Indices;
	delete[] Faces;
	delete[] Snorms;
	return 0; 
}

///
/// Checks if a segment connecting points
/// A and B intersects other segments in
/// the front. 
///
bool isIntersecting(
	Point *A, Point *B, // does segment A-B intersect a front?
	Vector &F, // surface-normal vector
	Node *node,  // do not consider this node
	Segment *seg, // do not consider this segment
	int level // level of recursion - how many 
	          // neighbor segments are considered
) {//TODO: fix this recursive implementation
	if (level <= 0) return false; 
	Collection<Segment> *segs = node->getSegs();
	if (segs->number() == 0) return false;
	segs->goFirst();
	do {
		Segment *next_seg = segs->getItem();
		if (next_seg != seg) {
			Node *next_node = next_seg->getOtherNode(node);
			if (
				intersecting(
					&(A->X), &(B->X),
					&(node->P->X), &(next_node->P->X),
					&F
				) ||
				isIntersecting(A, B, F, next_node, next_seg, level-1)
			) return true;
		}
		segs->goNext();
	} while(!segs->isFirst());
	return false;
}

/// 
/// Finds point with maximum bonding to a given point
///
REAL maxBond(// find closest point
	Point *P, // the given point
	Collection<Point> &points, // free points
	Segment *seg, // current segment
	Point *&Q // point with max bond to be returned
) {
	if (points.number() <= 0) return 0.0;

	Vector E;

	if (seg != nullptr) E.set(seg->getNorm());

	REAL max_score = 0.0;
	points.untag();
	points.goFirst();
	do {
		Point *A = points.getItem();

		if (dist(A->X, P->X) < P->range) {
			Vector D(A->X);
			D.sub(P->X);

			if (seg == nullptr || dot(D, E) > 0.0) { //exclude points
	                                      		// on the side of the mesh
				bool is_intersecting = false;
				REAL r = 1.0, p = 1.0;
				if (seg != nullptr) {
					// Check for intersection
					// between PA and neighboring segments
					Vector F;
					F.mid(P->F, A->F);
					for (int iv=0; iv<2; iv++) {
						Node *vert = seg->getNode(iv);
						if (is_intersecting = isIntersecting(P, A, F, vert, seg, 1)) break;
					}
					if (!is_intersecting) {
						Point
							*B = seg->getNodes()[0]->P,
							*C = seg->getNodes()[1]->P;
						REAL a = area(A, B, C);
						p = perimeter(A, B, C);
						r = a/p;
					}
				}//endif seg != nullptr
				if (!is_intersecting && r/p > MIN_ASPECT_RATIO) {//?
				REAL	
					b = P->bond(A),
					score = r*b;
					if (score > max_score) {
					max_score = score;
					points.setTag();
				}
				}
			}
		}
		points.goNext();
	} while(!points.isFirst());
	Q = points.getTaggedItem();
	return max_score;
}


/// 
/// Finds a front node with maximum bonding to a given point.
///
REAL maxBond(// find closest node on the front
	Point *P, // given point
	Collection<Node> &nodes, // front nodes 
	Segment *seg, // current segment
	Point *&Q // point with max bond on the front
) {

	Vector E;

	if (seg != nullptr) E.set(seg->getNorm());
 
	if (nodes.number() <= 0) return 0.0;
	Node **segnode = nullptr;
	if (seg!=nullptr) segnode= seg->getNodes();

	REAL max_score = 0.0;
	nodes.untag();
	nodes.goFirst();

	do { // for each node on the front:
		Node *node = nodes.getItem();
		REAL d = dist(P->X, node->P->X);
		if (d < P->range) {
			Vector D(node->P->X);
			D.sub(P->X);
			if (seg == nullptr || dot(D, E) > 0.0) { //exclude points
				                                   // behind the front
				int n=3;
				if (seg != nullptr) {
					for (n=0; n<3 && segnode[n] != node; n++);
				}
				if (n==3) { // exclude segment nodes 
					Point *A = node->P;
					REAL d = dist(A->X, P->X);
					if (d < P->range) {
						bool is_intersecting = false;
						if (seg != nullptr) {
							Vector F;
							F.mid(P->F, A->F);
							for (int iv=0; iv<2; iv++) {
								Node *vert = seg->getNode(iv);
								if (is_intersecting = isIntersecting(P, A, F, vert, seg, 1)) break;
							}
						}
						if (!is_intersecting) {
							REAL 
								b = A->bond(P),
								score = b;
							// Construct bond measure
							if (segnode != nullptr) {
								Point
									*B = segnode[0]->P,
									*C = segnode[1]->P;
								REAL 
									a = area(A, B, C),
									p = perimeter(A, B, C),
									r = a/p;
								if (r/p > MIN_ASPECT_RATIO) 
									score *= r;
								else
									score = 0.0;
							}
							if (score > max_score) {
								max_score = score;
								nodes.setTag();
							}
						}
					}
				}
			}
		}
		nodes.goNext();
	} while(!nodes.isFirst());
	Node *node = nodes.getTaggedItem();
	if (node == nullptr) return 0.0;
	Q = node->P;
	return max_score;
}

///
/// Moves node or segments from one collection to another.
///
void moveNode(Node *node, Collection<Node> &setA, Collection<Node> &setB) {
	if (node->numNeibs() > 0) return;
	setA.remove(node);
	setB.add(node);
}

/// Move current segment from set A to set B:
void moveSeg(Collection<Segment> &A, Collection<Segment> &B) {
	Segment *seg = A.getItem();
	A.remove();
	B.add(seg);
}

/// Move a given segment from A to B:
void moveSeg(Segment *seg, Collection<Segment> &A, Collection<Segment> &B) {
	A.remove(seg);
	B.add(seg);
}

#ifdef VTK
// 
// Output to VTK for debugging
// This will need VTK.cc file to link to
//
void write_vector (string filename, Collection<Segment> &segs) {

	extern void write_vector_vtp (const char *filename, INT nl, INT nv, REAL *P, REAL *V);

	int nl = segs.number();
	if (nl <= 0) return;

	if(!silent)cout << nl << " boundary segments"<<endl;

	REAL 
		*P = new REAL[nl*DIM], //coordinates of line vertices
		*V = new REAL[nl*DIM]; //line vectors

	int iseg = 0;
	segs.goFirst();
	do {
		Segment *seg = segs.getItem();
		Vector 
			A(seg->getNode(0)->P->X),
			B(seg->getNode(1)->P->X),
			C, AB;
		AB.sub(B,A);
		C.mid(A,B);
		for (int i=0; i<DIM; i++) {
			int n=DIM*iseg+i;
			P[n] = C.get(i);
			V[n] = AB.get(i);
		}
		iseg++;
		segs.goNext();
	} while (!segs.isFirst());

	write_vector_vtp (filename.c_str(), nl, DIM, P, V);

	delete[] P;
	delete[] V;
}
#endif

///
/// Generates Mesh Using Propagating Front Method.
///
int frontmesh 
(
	int nPoints, // number of points on the surface
	REAL *Coords, // points coordinates array
	REAL *Snorms, // surface normals at each point
	int *Indices, // points initial indices for output
	int *&Faces // point connectivity array to return
) {

	//Load points
	Collection<Point> points;

	for (int i=0; i<nPoints; i++) {
		int n = DIM*i,
			ind = Indices[i];
		REAL 
			*X = Coords + n,
			*F = Snorms + n;
		Point *P = new Point(ind, X, F);
		points.add(P);
	}

	if (!silent) cout << points.number() << " points loaded" << endl;

	Collection<Node> 
		fnodes, // front nodes
		mnodes; // mesh nodes - to avoid loose pointers
	Collection<Segment> 
		fsegs, // front segments
		bsegs, // boundary segments
		msegs; // mesh segments
	Collection<Triangle> mesh;

// 
// MAKE THE FIRST TRIANGLE
//
// Pick the first node
//
	points.goFirst();
	Point *pA = points.getItem();
	Node *nodeA = new Node(pA);
	fnodes.add(nodeA);
	points.remove();

	// Find closest point to A

	Point *pB = nullptr; 
	REAL bond = maxBond(pA, points, nullptr, pB);
	if (pB == nullptr) {
		cerr << "Cant locate second point - aborting"<<endl;
		exit(1);
	}
	Node *nodeB = new Node(pB);
	fnodes.add(nodeB);
	points.removeTagged();

	REAL range = RANGE*dist(nodeA->P->X, nodeB->P->X); 
	// adjust range if needed
	nodeA->P->setRange(range);
	nodeB->P->setRange(range);

	Vector 
		A(pA->X), 
		B(pB->X),
		AB(B), C, D, F, G; 
	AB.sub(A);
	C.mid(A, B); // middle between A and B

	F.mid(pA->F, pB->F); // surface normal vector
	Point pC(C, F);
	pC.setRange(range); 
	//Find closest point to C:
	Point *pE = nullptr; 
	bond = maxBond(&pC, points, nullptr, pE);
	if (pE == nullptr) {
		cerr << "Cannot find third point for the first triangle - aborting."<<endl;
		exit(1);
	}
	pE->setRange(range);
	Node *nodeE = new Node(pE);
	fnodes.add(nodeE); // add E to the front
	points.removeTagged(); // remove point pE from Points

	Segment // make new front segments:
		*segAB = new Segment(nodeA, nodeB, nodeE),
		*segAE = new Segment(nodeA, nodeE, nodeB),
		*segBE = new Segment(nodeB, nodeE, nodeA);

	nodeA->addSegs(segAB, segAE);
	nodeB->addSegs(segAB, segBE);
	nodeE->addSegs(segAE, segBE);
	// Add segments pA-pE, pB-pE to the front:
	fsegs.add(segAB);
	fsegs.add(segAE);// inserts behind current
	fsegs.add(segBE);
	// Add node pE and its neighbor segments to the front:
	// Make the first triangle:
	Triangle *triangle = new Triangle(pA, pB, pE);
	segAB->setFace0(triangle, 2);
	segAE->setFace0(triangle, 1);
	segBE->setFace0(triangle, 0);
	triangle->setSeg(2,segAB);
	triangle->setSeg(1,segAE);
	triangle->setSeg(0,segBE);

	mesh.add(triangle);
	range = triangle->perimeter();

	// Spread out range of influence uniformely
	points.goFirst();
	do {// TODO: Make range location-dependent
		points.getItem()->setRange(range);
		points.goNext();
	} while(!points.isFirst());

	//
	// MAIN LOOP OVER FRONT SEGMENTS
	//
	// Until there are still front segments left
	// loop over them:
	fsegs.goFirst();
	do {
		Segment *segAB = fsegs.getItem();
		nodeA = segAB->getNode(0),
		nodeB = segAB->getNode(1);
		pA = nodeA->P;
		pB = nodeB->P;
		A.set(pA->X);
		B.set(pB->X);
		C.mid(A, B);
		F.mid(pA->F, pB->F); // surface normal vector
		Point pC(C, F); 
		pC.setRange(0.5*(pA->range+pB->range));
		//Find closest point to C
		pE = nullptr;
		Point *pF = nullptr;
		REAL 
			bondE = maxBond(&pC, points, segAB, pE), // among remaining points
			bondF = maxBond(&pC, fnodes, segAB, pF); // among points on the front
			                           // excluding nodes on segment AB
		if (pE == nullptr && pF == nullptr) {
			// This segment has no neighbors = it's a boundary.
			// Remove if from the front and add to the boundary segments.
			bsegs.add(segAB);
			fsegs.remove();
			continue;
		}
		if (pE != nullptr && pF != nullptr) {// select strongest bond
			if (bondF > bondE) {
				pE = nullptr;
			} else {
				pF = nullptr;
			}
		}
		if (pF == nullptr) {// keep point E 
			nodeE = new Node(pE);
			segAE = new Segment(nodeA, nodeE, nodeB),
			segBE = new Segment(nodeB, nodeE, nodeA);
			nodeE->addSegs(segAE, segBE);
			// For nodeA: Replce current segment with segAE:
			nodeA->replace(segAB, segAE);
			// For nodeB: Replace current segment with segBE:
			nodeB->replace(segAB, segBE);
			fnodes.add(nodeE);
			Triangle *face = new Triangle(pA, pB, pE);
			segAE->setFace0(face, 1); // point pB is opposite to segAE 
			// and has index 1 inside the triangle (second argument above)
			segBE->setFace0(face, 0); // poiont pA is opposite to segBE
			// and has index 0 inside the triangle (first argument above)
			face->setSeg(2,segAB);
			face->setSeg(1,segAE);
			face->setSeg(0,segBE);
			mesh.add(face);
			segAB->setFace1(face, 2);
			segAB->setNadir(nodeE);
			moveSeg(fsegs, msegs);
//-			fsegs.remove(); // remove segment AB
			// Add segments pA-pE, pB-pE to the front:
			fsegs.insert(segAE);// inserts behind current
			fsegs.insert(segBE);
			points.removeTagged(); // removed point pE from Points
		} else {// keep point F (node on the front)
			Node *nodeF = fnodes.getTaggedItem();
			bool 
				neibA = nodeF->isNeib(nodeA), 
				neibB = nodeF->isNeib(nodeB);
			if (neibA && neibB) {
				// Isolated triangle 
				Segment 
					*segAF = nodeF->getSeg(nodeA),
					*segBF = nodeF->getSeg(nodeB);
				nodeA->remove(segAB);
				nodeA->remove(segAF);
				nodeB->remove(segAB);
				nodeB->remove(segBF);
				nodeF->remove(segAF);
				nodeF->remove(segBF);
				moveNode(nodeA, fnodes, mnodes);
				moveNode(nodeB, fnodes, mnodes);
				moveNode(nodeF, fnodes, mnodes);
				Triangle *face = new Triangle(pA, pB, pF);
				face->setSeg(2,segAB);
				face->setSeg(1,segAF);
				face->setSeg(0,segBF);
				mesh.add(face);
				segAB->setFace1(face, 2);
				segAB->setNadir(nodeF);
				moveSeg(fsegs, msegs);
//-				fsegs.erase();// erase segAB
				segAF->setFace1(face, 1);
				segAF->setNadir(nodeB);
				moveSeg(segAF, fsegs, msegs);
//-				fsegs.erase(segAF);
				segBF->setFace1(face, 0);
				segBF->setNadir(nodeA);
				moveSeg(segBF, fsegs, msegs);
		//-				fsegs.erase(segBF);
			} else
			if (neibA) {
				Segment 
					*segFA = nodeF->getSeg(nodeA),
					*segFB = new Segment(nodeF, nodeB, nodeA);
				nodeF->replace(segFA, segFB);
				nodeB->replace(segAB, segFB);
				nodeA->remove(segAB);
				nodeA->remove(segFA);
				moveNode(nodeA, fnodes, mnodes);
				// Add triangle FAB to mesh
				Triangle *face = new Triangle(pF, pA, pB);
				face->setSeg(0,segAB);
				face->setSeg(2,segFA);
				face->setSeg(1,segFB);
				segFB->setFace0(face, 1); // pA is opposite
				// to face FB and has index 1 (second argument)
				segAB->setFace1(face, 0);
				segAB->setNadir(nodeF);
				moveSeg(fsegs, msegs);
//-				fsegs.erase(); // erases current: segAB
				segFA->setFace1(face, 2);
				segFA->setNadir(nodeB);
				moveSeg(segFA, fsegs, msegs);
//-				fsegs.erase(segFA);
				fsegs.add(segFB);
				mesh.add(face);
			} else if (neibB) {// nodeF is a neighbor of B
				// For nodeF: replace segment FB with seg FA
				Segment 
					*segFB = nodeF->getSeg(nodeB),
					*segFA = new Segment(nodeF, nodeA, nodeB);
				nodeF->replace(segFB, segFA);
				// For nodeA: replace segment AB with FA
				nodeA->replace(segAB, segFA);
				nodeB->remove(segAB);
				nodeB->remove(segFB);
				moveNode(nodeB, fnodes, mnodes);
				// Add triangle FAB to mesh
				Triangle *face = new Triangle(pF, pA, pB);
				segFA->setFace0(face, 2); // pB is opposite to FA
				segAB->setFace1(face, 0);
				face->setSeg(0,segAB);
				face->setSeg(1,segFB);
				face->setSeg(2,segFA);
				segAB->setNadir(nodeF);
				moveSeg(fsegs, msegs);//move current seg which is segAB
//-				fsegs.erase(); // deletes current: segAB is current
				segFB->setFace1(face, 1);
				segFB->setNadir(nodeA);
				moveSeg(segFB, fsegs, msegs);
//-				fsegs.erase(segFB); // finds segFB and deletes is
				fsegs.add(segFA);
				mesh.add(face);
			} else {// nodeF is not neighbor of A or B
				Segment
					*segAF = new Segment(nodeA, nodeF, nodeB),
					*segBF = new Segment(nodeB, nodeF, nodeA);
				nodeF->addSegs(segAF, segBF);
				// Remove segAB from nodes A and B
				nodeA->remove(segAB);
				nodeB->remove(segAB);
				// Add segAF to nodeA
				nodeA->add(segAF);
				// Add segBF to nodeB
				nodeB->add(segBF);
				// Add triangle pA, pB, pF to mesh
				Triangle *face = new Triangle(pA, pB, pF);
				face->setSeg(2,segAB);
				face->setSeg(1,segAF);
				face->setSeg(0,segBF);
				segAF->setFace0(face, 1);
				segBF->setFace0(face, 0);
				mesh.add(face);
				// Erase segAB from fsegs
				segAB->setFace1(face, 2);
				segAB->setNadir(nodeF);
				moveSeg(fsegs, msegs);
//-				fsegs.erase(); // segAB is current seg
				fsegs.add(segAF);
				fsegs.add(segBF);
			}
		}
		fsegs.goNext();
	}	while(fsegs.number() != 0); 
	
	// 
	// MENDING HOLES
	// Zip holes between close boundary segments 
	//
/*
           * P
  A        D        B
  *--------*--------*
   \       |       /
    \      |      /
     \     |     /
      \    |    /
       \   |   /
        \  |  /
         \ | /
          \|/
           *
           C

Original triangle ABC is 
split into APC and PCB.

 */
	int npass=2, merr = 0, mwarn = 0; //Make nb passes to get them all
	for (int ipass=0; ipass<npass && bsegs.number() > 0; ipass++) {
		int nerr = 0, nwarn = 0;
		bsegs.goFirst();
		do {
			bool zipped = false;
			Segment *segAB = bsegs.getItem(); //segAB
			Triangle *triABC = segAB->getFace0();
			Node 
				*nodeA = segAB->getNode(0),
				*nodeB = segAB->getNode(1),
				*nodeC = segAB->getApex();
			Point 
				*pA = nodeA->P,
				*pB = nodeB->P,
				*pC = nodeC->P; // but we will get pC differently
			Vector
				*A = &(pA->X),
				*B = &(pB->X),
				AB, D;
			AB.sub(B, A);
			REAL ab = AB.len();
			AB.mul(1./ab); // normalize to 1
			D.mid(A, B);
			bsegs.setTag();
			bsegs.tagNext();
			do {
				Segment *that_seg = bsegs.getTaggedItem();
				for (int i=0; i<2; i++) {
					Node *nodeP = that_seg->getNode(i);
					Point *P = nodeP->P;
					if (dist(D, P->X)/ab < ZIP_PROXIMITY) {
						if (triABC == nullptr) {
							cerr << "ZIP HOLES WARNING: FAILED TO ZIP SOME HOLES IN PASS "<<ipass<<endl;
							nwarn++;
							break;
						}
						if (triABC->has(P)) continue;
						if (!triABC->replaceVert(pB, P)) {
							cerr << "ZIP HOLES ERROR: CANNOT REASSIGN TRIANGLE VERTEX IN PASS "<<ipass<<endl;
							nerr++;
							break;
						}
						segAB->setNode(1, nodeP); 
						//                             (0,  1,  2)
						Triangle *triPBC = new Triangle(P, pB, pC);
						Segment //            (    0,     1,     2)
							*segPB = new Segment(nodeP, nodeB, nodeC),
							*segPC = new Segment(nodeP, nodeC, nodeB); 
						segPB->setFace0(triPBC, 2); // 2->pC
						// segPB is a boundary seg => face1 is not set
						segPC->setFace0(triPBC, 1); // 1->pB
						int indA = triABC->getInd(pA);
//#ifdef ACERT
						if (indA < 0) {
							cerr << "ZIP HOLES ERROR: FAILED TO FIND VERTEX OF A TRIANGLE IN PASS "<<ipass<<endl;
							nerr++;
							break;
						}
//#endif
						segPC->setFace1(triABC, indA);
						Segment *segBC = triABC->getSeg(nodeB, nodeC);
						if (segBC == nullptr) {
							cerr << "ZIP HOLES WARNING: FAILED TO GET SEGMENT OF A TRIANGLE IN PASS"<<ipass<<endl;
							nwarn++;
							// TODO: fix this
//?							break;
						} else {
							segBC->replace(triABC, triPBC, 0);// 0->P
						}
						segBC->replaceCorner(nodeA, nodeP);
						segPC->setNadir(nodeA);
						triABC->replaceSeg(segBC, segPC);
						triPBC->setSeg(0,segBC); // 0->P
						triPBC->setSeg(1,segPC); // 1->pB
						triPBC->setSeg(2,segPB); // 2->pC
						bsegs.add(segPB);
						msegs.add(segPC);
						mesh.add(triPBC);
						zipped = true;
						break;
					}
				}
				if (zipped) break;
				bsegs.tagNext();
			} while (!bsegs.tagIsCurrent());
			bsegs.goNext();
		} while (!bsegs.isFirst());
		if (nwarn > 0 || nerr > 0) {
			cerr << "Pass "<< ipass << " of zip-holes produced "<<nerr<<" errors and "<<nwarn<<" warnings\n";
		}
		merr += nerr;
		mwarn += nwarn;
	}// end zip holes
//#define DEBUG
#ifdef DEBUG
	//
	// MESH OPTIMIZARION
	// Increase area/perimeter ratios
	//
	if (merr == 0 && mwarn == 0) {
		const REAL max_angle_cos = 0.9; //max cos(angle) between triangles
		int 
			niter = 0, 
			nswaps = 0, mswaps = 0;
			merr = 0, mwarn = 0;
		if (!silent)cout<<"Mesh optimization"<<endl;
		do {// niter-loop
			int nerr = 0, nwarn = 0;
			nswaps = 0;
			msegs.goFirst();
			do {// looping over msegs
				Segment *segBC = msegs.getItem();
				bool complete = true;
				for (int i=0; i<2; i++) {
					if (segBC->getFace(i) == nullptr) {
						complete = false;
						break;
					}
				}
				if (complete) {
					Triangle tris[2];
					int ivert[2]; //indexes of opposite vertices
					for (int i=0; i<2; i++) {
						tris[i].copy(segBC->getFace(i));
						ivert[i] = segBC->getVertInd(i);
					}
					Vector Norm[2]; // surface normal vectors
					REAL ap_min_old = BIG; // area-to-perimeter ratio
					for (int i=0; i<2; i++) {
						Triangle *tri = tris+i;
						tri->areaVec(Norm[i]); // now it's area
						REAL 
							area = Norm[i].len(),
							ap = area / tri->perimeter();
						if (ap < ap_min_old) {
							ap_min_old = ap;
						}
						Norm[i].mul(1./area); // normalize to unity
					}
					REAL angle = dot(Norm[0],Norm[1]);
					if (angle > max_angle_cos) {
				//
				// Swap segment
				//
	/*
	        B
	        .
	       / \
	      / ^ \
	     /  |  \
	   A. <-+-> .D
	     \  |  /
	      \ v /
	       \ /
	        .
	        C
	    BC -> AD
	*/
						int 
							iA = ivert[0],// A is opposite to segBC in triABC
							iB = (iA+1)%3,
							iC = (iA+2)%3,
							jD = ivert[1];// index of D in BCD (D is opposite to BC)
						Point
							*A = tris[0].getVert(iA),
							*B = tris[0].getVert(iB),
							*C = tris[0].getVert(iC),
							*D = tris[1].getVert(jD);
						int
							jB = tris[1].getInd(B),// index of B n triangle BCD
							jC = tris[1].getInd(C);// index of C in BCD
						if (jB < 0 || jC < 0) {
							cerr << "ERROR in mesh optimization: Failed to get vertex index\n";
							nerr++;
							break;
						}
						Vector
							AB, BC, AD, AC, BD,
							areaABC, areaBDC, // old triangles 
							areaABD, areaADC; // new triangles
						AB.sub(B->X, A->X);
						BC.sub(C->X, B->X);
						AD.sub(D->X, A->X);
						AC.sub(C->X, A->X);
						BD.sub(D->X, B->X);
						areaABC.cross(&AB, &AC);
						areaBDC.cross(&BD, &BC);
						areaABD.cross(&AB, &AD);
						areaADC.cross(&AD, &AC);
						REAL 
							ap_old = dot(areaABC, areaBDC)/BC.len(),
							ap_new = dot(areaABD, areaADC)/AD.len();

						if (ap_new > ap_old) {	
		
							// Triangle: ABC -> ABD
							tris[0].setVert(iC, D); // C->D
							ivert[0] = iB;  // since B is opposite to AD in ABD
							// Triangle: DCB -> DCA: 
							tris[1].setVert(jC, C);
							tris[1].setVert(jB, A);
							ivert[1] = jC;

//?						REAL ap_min_new = BIG;
//?						for (int i=0; i<2; i++) {
//?							Triangle *tri = tris+i;
//?							tri->areaVec(Norm[i]); // now it's area
//?							REAL 
//?								area = Norm[i].len(),
//?								ap = area / tri->perimeter();
//?							if (ap < ap_min_new) {
//?								ap_min_new = ap;
//?							}
//?							Norm[i].mul(1./area); // normalize to unity
//?						}
//?						REAL angle_new = dot(Norm[0],Norm[1]);
//?						if (angle_new >= angle_old && ap_min_new > ap_min_old) {// do swap
//-cout << "niter="<<niter<<", ido1="<<ido1<<", nswaps: "<<nswaps<<", ap: "<<ap_old<<" -> "<<ap_new<<endl;cout.flush();//-
							// Swap segment: BC -> AD
							// Get access to nodes A, B, C, D
							Node 
								*nodeB = segBC->getNode(B),
								*nodeC = segBC->getNode(C),
								*nodeA = segBC->getCorner(A),
								*nodeD = segBC->getCorner(D);
//-A->show("Point A: ");
//-B->show("Point B: ");
//-C->show("Point C: ");
//-D->show("Point D: ");
//-SegBC->show("segBC:");//-
//-Cout << "point A:"<<A<<", node A:"<<nodeA<<endl;
//-NodeA->show("nodeA.show:");//-
//-Cout << "node B:"<<nodeB<<" <-> C:"<<nodeC<<endl;//-
//-NodeB->show("nodeB.show:");//-
//-NodeC->show("nodeC.show:");//-
					
							// Access segments
							Segment 
								*segAB = segBC->getFace(0)->getSeg(nodeA, nodeB),
								*segAC = segBC->getFace(0)->getSeg(nodeA, nodeC),
								*segBD = segBC->getFace(1)->getSeg(nodeB, nodeD),
								*segCD = segBC->getFace(1)->getSeg(nodeC, nodeD);
	
							if (segAC == nullptr || segBD == nullptr) {
								cerr << "MESH OPTIMIZER WARNING "<<nwarn<<": FAILED TO GET SEGMENT IN PASS "<<niter<<endl;

//-cerr.flush();cout<<"segAB="<<segAB<<", segAC="<<segAC<<", segBD="<<segBD<<", segCD="<<segCD<<endl;cout.flush();//-
//-segBC->getFace(0)->show("segBC->face(0)");
//-segBC->getFace(1)->show("segBC->face(1)");
								nwarn++;
							} else {
								// In segment BC replace nodeB -> nodeA, nodeC -> nodeD, apex -> nodeB
								Segment *&segAD = segBC; // renaming for clarity
								segAD->replaceNode(nodeB, nodeA);
								segAD->replaceNode(nodeC, nodeD);
								segAD->setApex(nodeB);
								segAD->setNadir(nodeC);
		
								// Segment AC should repoint from ABC to ACD
								segAC->swap(segAD->getFace(0), segAD->getFace(1),jD);
								segAC->replaceCorner(nodeB, nodeD);
								// Segment BD should repoint from BCD to ABD
								segBD->swap(segAD->getFace(1), segAD->getFace(0),iA);
								segBD->replaceCorner(nodeC, nodeA);

								segAB->replaceCorner(nodeC, nodeD);
								segCD->replaceCorner(nodeB, nodeA);

								// AB,BC,Ac -> AB,BD,AD
								tris[0].setSeg(0, segAB);
								tris[0].setSeg(1, segBD);
								tris[0].setSeg(2, segAD);
								// BC,CD,BD -> AC,CD,AD
								tris[1].setSeg(0, segAC);
								tris[1].setSeg(1, segCD);
								tris[1].setSeg(2, segAD);
		
								for (int i=0; i<2; i++) {
									Triangle *tri = segBC->getFace(i);
									tri->copy(tris+i);
									segBC->setVertInd(i, ivert[i]);
								}
								nswaps++;
							}
						}
					}// endif max_angle_cos
				}// endif complete
				if (nerr>0) break;
				msegs.goNext();
			} while (!msegs.isFirst());
			if (nerr > 0 || nwarn > 0) {
				cerr << "Mesh optimization issued "<<nerr<<" errors and "<<nwarn<<" warnings in pass "<<niter<<endl;
				merr += nerr;
				mwarn += nwarn;
			}
			if (verbose && nswaps > 0) {
				cout << "Swapped "<<nswaps<<" segments in "<<niter<<" pass"<<endl;
			}
			niter++;
			mswaps += nswaps;
			if (nerr>0) break;
		} while(nswaps > 0); //end niter
		if (merr > 0 || mwarn > 0) {
			cout << "Mesh optimization issued "<<merr<<" errors and "<<mwarn<<" warnings"<<endl;
		}
		if (!silent && mswaps > 0) {
			cout << "Swapped "<<mswaps<<" segments in "<<niter<<" passes"<<endl;
		}
	}
#endif
	//
	// Pack all triangles into a single array for output
	//
	int nFaces = mesh.number();
	if(!silent) cout << "Constructed mesh of "<<nFaces<<" triangles" << endl;
	Faces = new int[NV*nFaces];
	mesh.goFirst();
	int n = 0;
	do {
		Triangle *face = mesh.getItem();
		for (int j=0; j<NV; j++) {
			int ind = face->getVert(j)->i;
			Faces[n+j] = ind;
		}
		n += NV;
		mesh.goNext();
	}	while(!mesh.isFirst());

	if(!silent) {
		cout 
		<< "Number of mesh triangles: "<<mesh.number()<<endl
		<< "Number of mesh edges: " << msegs.number()<<endl
		<< "Number of boundary edges: "<< bsegs.number()<<endl
		<< "Number of mesh nodes: "<< mnodes.number()<<endl;
	}

	fnodes.clean();
	fsegs.clean();
	mesh.wipe();
	bsegs.wipe();
	mnodes.wipe();
	points.wipe();
	msegs.clean();

	return nFaces;
}

