
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


//
// Main routine with IO
//
int main ( int argc, char *argv[] )
{
	int frontmesh (int nPoints, REAL *Coords, REAL *Normals, int *Indices, int *&Faces); 
	silent = false;
	verbose = false;
	enum Format {dat_format, node_format} io_format;
	int c;
	while ((c = getopt (argc, argv, "ahsv")) != -1)
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
		   strstr(points_filename, ".dat")==NULL 
		&& strstr(points_filename, ".node")==NULL) {
		cerr << "File "<<points_filename<<" needs suffix '.node' or '.dat'" 
		<< endl;
		exit(1);
	}
	io_format = dat_format;
	if (strstr(points_filename, ".node")!=NULL) io_format = node_format;

	if (strstr(output_filename, ".face")==NULL) {
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
		*Faces = NULL;

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
	delete Points;
	delete Indices;
	delete[] Faces;
	return 0; 
}

//
// Checking for intersections
//
bool isIntersecting(
	Point *A, Point *B, // does segment A-B intersect a front?
	Segment *seg0, // do not consider this segment
	Segment *seg, // do not consider this segment
	Node *node,  // do not consider this node
	int level // level of recursion - how many 
	          // neighbor segments are considered
) {
	if (level <= 0) return false; 
	for (int i=0; i<2; i++) {
		Node *vert = seg->getNode(i);
		if (vert == node) continue;
		Collection<Segment> *vertsegs = vert->getSegs();

		vertsegs->goFirst();
		do {
			Segment *segnext = vertsegs->getItem();
			if (segnext != seg && segnext != seg0) {
				Node *vertneib = segnext->getOtherNode(vert);
				if (
					intersecting(
						&(A->X), &(B->X),
						&(vert->P->X), &(vertneib->P->X),
						&(A->F)
					) 
				) return true;
				if (isIntersecting(A, B, seg, segnext, vertneib, --level)) return true;
			}
			vertsegs->goNext();
		} while(!vertsegs->isFirst());
	}
	return false;
}

// 
// Find point with maximum bonding to a given point
//
REAL maxBond(// find closest point
	Point *P, // the given point
	Collection<Point> &points, // free points
	Segment *seg, // current segment
	Point *&Q // point with max bond to be returned
) {
	if (points.number() <= 0) return 0.0;

	Vector E;

	if (seg != NULL) E.set(seg->getNorm());

	REAL max_score = 0.0;
	points.untag();
	points.goFirst();
	do {
		Point *A = points.getItem();

		if (dist(A->X, P->X) < P->range) {
			Vector D(A->X);
			D.sub(P->X);

			if (seg == NULL || dot(D, E) > 0.0) { //exclude points
	                                      		// on the side of the mesh
				bool is_intersecting = false;
				REAL r = 1.0, p = 1.0;
				if (seg != NULL) {
					// Check for intersection
					// between PA and neighboring segments
					is_intersecting = isIntersecting(P, A, seg, seg, NULL, 3);
					if (!is_intersecting) {
						Point
							*B = seg->getNodes()[0]->P,
							*C = seg->getNodes()[1]->P;
						REAL a = area(A, B, C);
						p = perimeter(A, B, C);
						r = a/p;
					}
				}//endif seg != NULL
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


// 
// Find front node with maximum bonding to a given point
//
REAL maxBond(// find closest node on the front
	Point *P, // given point
	Collection<Node> &nodes, // front nodes 
	Segment *seg, // current segment
	Point *&Q // point with max bond on the front
) {

	Vector E;

	if (seg != NULL) E.set(seg->getNorm());
 
	if (nodes.number() <= 0) return 0.0;
	Node **segnode = NULL;
	if (seg!=NULL) segnode= seg->getNodes();

	REAL max_score = 0.0;
	nodes.untag();
	nodes.goFirst();

	do { // for each node on the front:
		Node *node = nodes.getItem();
		REAL d = dist(P->X, node->P->X);
		if (d < P->range) {
			Vector D(node->P->X);
			D.sub(P->X);
			if (seg == NULL || dot(D, E) > 0.0) { //exclude points
				                                   // behind the front
				int n=3;
				if (seg != NULL) {
					for (n=0; n<3 && segnode[n] != node; n++);
				}
				if (n==3) { // exclude segment nodes 
					Point *A = node->P;
					REAL d = dist(A->X, P->X);
					if (d < P->range) {
						bool is_intersecting = false;
						if (seg != NULL) {
							is_intersecting = isIntersecting(P, A, seg, seg, NULL, 3);
						}
						if (!is_intersecting) {
							REAL 
								b = A->bond(P),
								score = b;
							// Construct bond measure
							if (segnode != NULL) {
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
	if (node == NULL) return 0.0;
	Q = node->P;
	return max_score;
}

//
// Removes node and puts it on the mesh
//
void moveNode(Node *node, Collection<Node> &setA, Collection<Node> &setB) {
	if (node->numNeibs() > 0) return;
	setA.remove(node);
	setB.add(node);
}

#ifdef VTK
// 
// Output to VTK for debugging
//
void write_lines (string filename, Collection<Segment> &segs) {

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

//
// Create Mesh Using Propagating Front
//
int frontmesh (int nPoints, REAL *Coords, REAL *Snorms, int *Indices, int *&Faces) {

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
		bsegs; // boundary segments
	Collection<Triangle> mesh;

	points.goFirst();
	Point *pA = points.getItem();
	Node *nodeA = new Node(pA);
	fnodes.add(nodeA);
	points.remove();

	// Find closest node

	Point *pB = NULL;
	REAL bond = maxBond(pA, points, NULL, pB);
	if (pB == NULL) {
		cerr << "Cant locate second point - aborting"<<endl;
		exit(1);
	}
	Node *nodeB = new Node(pB);
	fnodes.add(nodeB);
	points.removeTagged();

	REAL range = RANGE*dist(nodeA->P->X, nodeB->P->X); 
	//TODO: adjust range if needed
	nodeA->P->setRange(range);
	nodeB->P->setRange(range);

	Vector 
		A(pA->X), 
		B(pB->X),
		AB(B), C, D, F, G; 
	AB.sub(A);
	C.mid(A, B);

	F.mid(pA->F, pB->F); // surface normal vector
	Point pC(C, F);
	pC.setRange(range); 
	//Find closest point to C
	Point *pE = NULL;
	bond = maxBond(&pC, points, NULL, pE);
	if (pE == NULL) {
		cerr << "Cannot find third point for the first triangle - aborting."<<endl;
		exit(1);
	}
	pE->setRange(range);
	Node *nodeE = new Node(pE);
	fnodes.add(nodeE);
	points.removeTagged(); // removed point pE from Points

	Segment 
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
	// Add points A, B, E to mesh
	Triangle *triangle = new Triangle(pA, pB, pE);
	segAB->setTriangle(triangle, 2);
	segAE->setTriangle(triangle, 1);
	segBE->setTriangle(triangle, 0);
	mesh.add(triangle);
	range = triangle->perimeter();

	points.goFirst();
	do {
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
		pE = NULL;
		Point *pF = NULL;
		REAL 
			bondE = maxBond(&pC, points, segAB, pE), // among remaining points
			bondF = maxBond(&pC, fnodes, segAB, pF); // among points on the front
			                           // excluding nodes on segment AB
		if (pE == NULL && pF == NULL) {
			// This segment has no neighbors = it's a boundary.
			// Remove if from the front and add to the boundary segments.
			bsegs.add(segAB);
			fsegs.remove();
			continue;
		}
		if (pE != NULL && pF != NULL) {// select strongest bond
			if (bondF > bondE) {
				pE = NULL;
			} else {
				pF = NULL;
			}
		}
		if (pF == NULL) {// keep point E 
			nodeE = new Node(pE);
			segAE = new Segment(nodeA, nodeE, nodeB),
			segBE = new Segment(nodeB, nodeE, nodeA);
			nodeE->addSegs(segAE, segBE);
			// For nodeA: Replce current segment with segAE:
			nodeA->replace(segAB, segAE);
			// For nodeB: Replace current segment with segBE:
			nodeB->replace(segAB, segBE);
			// Remove current segment:
			fsegs.erase(); // erase segment AB
			// Add segments pA-pE, pB-pE to the front:
			fsegs.insert(segAE);// inserts behind current
			fsegs.insert(segBE);
			fnodes.add(nodeE);
			Triangle *face = new Triangle(pA, pB, pE);
			segAE->setTriangle(face, 1); // point pB is opposite to segAE 
			// and has index 1 inside the triangle (second argument above)
			segBE->setTriangle(face, 0); // poiont pA is opposite to segBE
			// and has index 0 inside the triangle (first argument above)
			mesh.add(face);
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
				mesh.add(face);
				fsegs.erase();// erase segAB
				fsegs.erase(segAF);
				fsegs.erase(segBF);
			} else
			if (neibA) {
				Segment 
					*segFA = nodeF->getSeg(nodeA),
					*segFB = new Segment(nodeF, nodeB, nodeA);
				nodeF->replace(segFA, segFB);
				nodeB->replace(segAB, segFB);
				nodeA->remove(segAB);
				nodeA->remove(segFA);
				// Erase segments FA and AB from fsegs
				fsegs.erase(); // erases current: segAB
				fsegs.erase(segFA);
				fsegs.add(segFB);
				moveNode(nodeA, fnodes, mnodes);
				// Add triangle FAB to mesh
				Triangle *face = new Triangle(pF, pA, pB);
				segFB->setTriangle(face, 1); // pA is opposite
				// to face FB and has index 1 (second argument)
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
				// Erase segments FB and AB from fsegs
				fsegs.erase(); // deletes current: segAB is current
				fsegs.erase(segFB); // finds segFB and deletes is
				fsegs.add(segFA);
				moveNode(nodeB, fnodes, mnodes);
				// Add triangle FAB to mesh
				Triangle *face = new Triangle(pF, pA, pB);
				segFA->setTriangle(face, 2); // pB is opposite to FA
				// inside triangle and has index 2 (third argument above)
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
				segAF->setTriangle(face, 1);
				segBF->setTriangle(face, 0);
				mesh.add(face);
				// Erase segAB from fsegs
				fsegs.erase(); // segAB is current seg
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
	int nb=2; //Ideally, ib=1 should be enough
	for (int ib=0; ib<nb && bsegs.number() > 0; ib++) {

	// Mend holes
	bsegs.goFirst();
	do {
		bool zipped = false;
		Segment *this_seg = bsegs.getItem();
		Node 
			*nodeA = this_seg->getNode(0),
			*nodeB = this_seg->getNode(1),
			*nodeC = this_seg->getApex();
		Point 
			*pA = nodeA->P,
			*pB = nodeB->P;
		Vector
			*A = &(pA->X),
			*B = &(pB->X),
			AB, C;
		AB.sub(B, A);
		REAL ab = AB.len();
		AB.mul(1./ab); // normalize to 1
		C.mid(A, B);
		bsegs.setTag();
		bsegs.tagNext();
		do {
			Segment *that_seg = bsegs.getTaggedItem();
			for (int i=0; i<2; i++) {
				Node *nodeP = that_seg->getNode(i);
				Point *P = nodeP->P;
				if (dist(C, P->X)/ab < ZIP_PROXIMITY) {
					// Get triangle of this_seg:
					Triangle *triangle = this_seg->getTriangle();
					if (triangle == NULL) {
						cerr << "FAILED TO ZIP SOME HOLES"<<endl;
						break;
					}
					// Re-assign triangle->pB to P
					if (!triangle->replaceVert(pB, P)) {
						cerr << "ERROR: CANNOT REASSIGN TRIANGLE VERTEX"<<endl;
						break;
					}
					this_seg->setNode(1, nodeP);
					// Form new Triangle: triangle->pA, pB, P
					Triangle *new_triangle =
						new Triangle(P, pB, this_seg->getVertex());
					// and add it to the mesh
					mesh.add(new_triangle);
					Segment *seg = new Segment(nodeA, nodeP, nodeC);
					seg->setTriangle(new_triangle, 2);
					bsegs.add(seg);
					zipped = true;
					break;
				}
			}
			if (zipped) {
//?			bsegs.erase(); //TODO: this should work
				break;
			}
			bsegs.tagNext();
		} while (!bsegs.tagIsCurrent());
		bsegs.goNext();
	} while (!bsegs.isFirst());
	}// end zip holes
#ifdef VTK
	write_lines("boundary.vtp", bsegs); 
#endif

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

	fnodes.clean();
	fsegs.clean();
	mesh.wipe();
	bsegs.wipe();
	mnodes.wipe();
	points.wipe();

	return nFaces;
}

