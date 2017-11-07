/** \file
 * \brief Header file with definitions of basic data structures.
 */
#define INT	int
#define REAL double
#define BIG	1.e30
#define SMALL	1.e-30
#define DIM 3
#define NV 3 

#define MIN_ASPECT_RATIO 0.005
#define ZIP_PROXIMITY 0.2
#define RANGE	1.2 // sensitive to mesh density

//#define VTK

using namespace std;

extern bool silent, verbose;

/**
 * Vector class.
 * Implements vector operations.
 */
class Vector {
	REAL V[DIM];

public:

	Vector() {
		for (int i=0; i<DIM; i++)
			V[i] = 0.0;
	}

	~Vector() {
	}

	Vector(REAL x, REAL y, REAL z) {
		V[0] = x;
		V[1] = y;
		V[2] = z;
	}

	Vector(REAL X[]) {
		for (int i=0; i<DIM; i++) {
			V[i] = X[i];
		}
	}

	Vector(Vector *X) {
		if (X==nullptr) return;
		for (int i=0; i<DIM; i++) {
			V[i] = X->get(i);
		}
	}

	Vector(Vector &X) {
		for (int i=0; i<DIM; i++) {
			V[i] = X.get(i);
		}
	}

	REAL *vec() {
		return V;
	}

	inline void set(int i, REAL x) {
		if (i>=0 && i<DIM)
			V[i] = x;
	}
	inline void set(REAL *X) {
		if (X==nullptr) return;
		for (int i=0; i<DIM; i++)
			V[i] = X[i];
	}

	inline REAL get(int i) {
		if (i>=0 && i<DIM)
			return V[i];
	}

	inline void set(Vector *X) {
		if (X==nullptr) return;
		for (int i=0; i<DIM; i++)
			V[i] = X->get(i);
	}

	inline void set(Vector &X) {
		for (int i=0; i<DIM; i++)
			V[i] = X.get(i);
	}

	inline void add(Vector &X) {
		for (int i=0; i<DIM; i++) {
			V[i] += X.get(i);
		}
	}

	inline void add(Vector *X) {
		for (int i=0; i<DIM; i++) {
			V[i] += X->get(i);
		}
	}

	inline void sub(Vector *X) {
		if (X==nullptr) return;
		for (int i=0; i<DIM; i++)
			V[i] -= X->get(i);
	}

	inline void sub(Vector &X) {
		for (int i=0; i<DIM; i++) {
			V[i] -= X.get(i);
		}
	}

	inline void sub(Vector &X, Vector &Y) {
		for (int i=0; i<DIM; i++) {
			V[i] = X.get(i) - Y.get(i);
		}
	}
	inline void sub(Vector &X, Vector *Y) {
		for (int i=0; i<DIM; i++) {
			V[i] = X.get(i) - Y->get(i);
		}
	}
	inline void sub(Vector *X, Vector *Y) {
		for (int i=0; i<DIM; i++) {
			V[i] = X->get(i) - Y->get(i);
		}
	}

	inline void mid(Vector &A, Vector &B) {
		for (int i=0; i<DIM; i++) {
			V[i] = 0.5*(A.get(i) + B.get(i));
		}
	}

	inline void mid(Vector *A, Vector *B) {
		for (int i=0; i<DIM; i++) {
			V[i] = 0.5*(A->get(i) + B->get(i));
		}
	}

	inline void mul(REAL c) {
		for (int i=0; i<DIM; i++) {
			V[i] *= c;
		}
	}

	inline void rot(REAL angle, Vector &axis) {
//TODO: complete if needed
	}

	inline REAL len() {
		REAL size = 0.0;
		for (int i=0; i<DIM; i++) {
			REAL v = V[i];
			size += v*v;
		}
		return sqrt(size);
	}
	
	inline void one() {// normalize to unity
		mul(1./len());
	}

	inline void cross(Vector *A, Vector *B) {
		REAL 
			*a = A->vec(),
			*b = B->vec();
		V[0] = a[1]*b[2] - a[2]*b[1];
		V[1] = a[2]*b[0] - a[0]*b[2];
		V[2] = a[0]*b[1] - a[1]*b[0];
	}

	void show(string name) {
		cout << name << " V:";
		for (int i=0;i<3;i++) cout << ' '<<V[i];
		cout << endl;
	}
};


inline REAL dist(Vector &A, Vector &B) {
	REAL d=0.0;
	for (int i=0; i<DIM; i++) {
		REAL r=A.get(i)-B.get(i);
		d+=r*r;
	}
	return sqrt(d);
}

inline REAL dist(Vector *A, Vector *B) {
	REAL d=0.0;
	for (int i=0; i<DIM; i++) {
		REAL r=A->get(i)-B->get(i);
		d+=r*r;
	}
	return sqrt(d);
}

inline REAL dot(Vector &A, Vector &B) {
	REAL c = 0.0;
	for (int i=0; i<DIM; i++) {
		c += A.get(i)*B.get(i);
	}
	return c;
}

inline REAL dot(Vector &A, Vector *B) {
	REAL c = 0.0;
	for (int i=0; i<DIM; i++) {
		c += A.get(i)*B->get(i);
	}
	return c;
}

inline REAL dot(Vector *A, Vector *B) {
	REAL c = 0.0;
	for (int i=0; i<DIM; i++) {
		c += A->get(i)*B->get(i);
	}
	return c;
}

inline void add(Vector &A, Vector &B, Vector &C) {
	REAL 
		*a = A.vec(),
		*b = B.vec(),
		*c = C.vec();
	for (int i=0; i<DIM; i++) 
		c[i] = a[i] + b[i];
}

inline void mid(Vector &A, Vector &B, Vector &C) {
	REAL 
		*a = A.vec(),
		*b = B.vec(),
		*c = C.vec();
	for (int i=0; i<DIM; i++) {
		c[i] = 0.5*(a[i] + b[i]);
	}
}

inline void sub(Vector &A, Vector &B, Vector &C) {
	REAL 
		*a = A.vec(),
		*b = B.vec(),
		*c = C.vec();
	for (int i=0; i<DIM; i++) 
		c[i] = a[i] - b[i];
}
inline void mul(REAL c, Vector &A) {
	REAL *a = A.vec();
	for (int i=0; i<DIM; i++) a[i] *= c;
}

inline void cross(Vector &A, Vector &B, Vector &C) {
	REAL 
		*a = A.vec(),
		*b = B.vec(),
		*c = C.vec();
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

bool clockwise(Vector *A, Vector *B, Vector *C, Vector *D) {
	Vector AB, BC, X;
	AB.sub(B,A);
	BC.sub(C,B);
	X.cross(&AB,&BC);
	if(dot(X,D)>=0.0) return true;
	return false;
} 

inline bool intersecting(
	Vector *A,
	Vector *B,
	Vector *C,
	Vector *D,
	Vector *E // direction normal to the surface
) { // Check intersectin between AB and CD
//http://jeffe.cs.illinois.edu/teaching/373/notes/x06-sweepline.pdf
//if (A==nullptr || B==nullptr || C == nullptr || D == nullptr || E == nullptr) return false;
	if (clockwise(A,C,D,E) == clockwise(B,C,D,E)) return false;
	if (clockwise(A,B,C,E) == clockwise(A,B,D,E)) return false;
	return true;
}

/** 
 * Point data strucure combines position 
 * vector with the surface normal vector.
 */
struct Point {
	REAL range = 0.0; // range of influence
	int i = 0; // index from the input
	Vector X, F; // coordinates and force
	Point(Vector &newX, Vector &newF) :
		range(0.0),
		i(0),
		X(),
		F()
	{
		range = BIG;
		i = -1;
		X.set(newX);
		F.set(newF);
	}
	Point(int ind, REAL *newX, REAL *newF) :
		X(), F()
	{
		range = BIG;
		i = ind;
		X.set(newX);
		F.set(newF);
	}
	void set(Vector &newX, Vector &newF) {
		X.set(newX);
		F.set(newF);
	}
	void setRange(REAL r) {
		range = r;
	}
	REAL bond (Point *Q) {
		// Bonding force to another point 
		REAL d = dist(X, Q->X);
		if (d < SMALL || d > range) return 0.0;
		return dot(F, Q->F) / d;
	}
	void show(string name) {
		cout << name << "ind="<<i<<", X:";
		for (int i=0; i<DIM; i++)
			cout << ' ' << X.get(i);
		cout << ", F:";
		for (int i=0; i<DIM; i++)
			cout << ' ' << F.get(i);
		cout << endl;
	}

};

inline REAL perimeter(Point *A, Point *B, Point *C) {
		return dist(A->X,B->X)+dist(B->X,C->X)+dist(C->X,A->X);
}

REAL area(Point *A, Point *B, Point *C) {
		Vector AB, AC, Area;
		if (A == nullptr || B == nullptr || C == nullptr) return 0.0;
		AB.set(B->X);
		AC.set(C->X);
		AB.sub(A->X);
		AC.sub(A->X);
		cross(AB, AC, Area);
		REAL a = 0.5*Area.len();
		return a;
}

/**
 * Collection class implements linked list
 * operations on other data strutures.
 */
template <class Item>
class Collection {
		struct List {
				Item *item;
				List *next, *prev;
				List() :
						item(nullptr),
						next(nullptr),
						prev(nullptr)
				{
						item = nullptr;
						next = prev = nullptr;
				}
				~List() {
						item = nullptr; 
						next = prev = nullptr;
				}
				void set(Item *item) {
						this.item = item;
				}
		};
		int n;
		List *current, *first, *tag;
		public:
		Collection() :
				n(0),
				current(nullptr),
				first(nullptr),
				tag(nullptr)
		{
		}
		~Collection() {
				if (n > 0 &&  first != nullptr) {
						current = first->next;
						while (current != first && current != nullptr) {
								current = current->next;
								delete current->prev;
						}
						if (current != nullptr) delete current;
				}
		}
		int number() { return n; }
		void add(Item *item) {
				if(n==0){
						first=new List();
						first->item=item;
						first->prev=first->next=first;
						current=first;
				} else {//Add before the first:
						List *newelement=new List();
						newelement->item=item;
						newelement->next=first;
						newelement->prev=first->prev;
						first->prev->next=newelement;
						first->prev=newelement;
				}
				n++;
		}
		void push(Item *item) {
				if(n==0){
						first=new List();
						first->item=item;
						first->prev=first->next=first;
						current=first;
				} else {//Add befpre the current:
						List *newelement=new List();
						newelement->item=item;
						newelement->prev=current;
						newelement->next=current->next;
						current->next->prev=newelement;
						current->next=newelement;
				}
				n++;
		}
		void insert(Item *item) {
				if(n==0){
						first=new List();
						first->item=item;
						first->prev=first->next=first;
						current=first;
				} else {//Add after the current:
						List *newelement=new List();
						newelement->item=item;
						newelement->next=current;
						newelement->prev=current->prev;
						current->prev->next=newelement;
						current->prev=newelement;
				}
				n++;
		}
		void append(Item *item) {
				if(n==0){
						first=new List();
						first->item=item;
						first->prev=first->next=first;
						current=first;
				} else {//Add after the first:
						List *newelement=new List();
						newelement->item=item;
						newelement->next=first;
						newelement->prev=first->prev;
						first->prev->next=newelement;
						first->prev=newelement;
				}
				n++;
		}
		void remove() {// unlinks without deleting the item
				if(n==0) return;
				if(current==nullptr) {
						cerr << "CANNOT DELETE VERTEX" << endl;
						exit(1);
				}
				current->item=nullptr;
				if(current==tag) {
						tag=nullptr;
				}
				if(current==first)first=first->next;
				current->prev->next=current->next;
				current->next->prev=current->prev;
				List *dead = current;
				current=current->next;
				delete dead;
				if(--n==0){
						current=first=tag=nullptr;
				}
		}
		void erase() {// unlinks and deletes the item
				if(n==0) return;
				if(current==nullptr) {
						cerr << "CANNOT DELETE VERTEX" << endl;
						exit(1);
				}
				delete current->item;
				current->item=nullptr;
				if(current==tag) {
						tag=nullptr;
				}
				if(current==first)first=first->next;
				current->prev->next=current->next;
				current->next->prev=current->prev;
				List *dead = current;
				current=current->next;
				delete dead;
				if(--n==0){
						current=first=tag=nullptr;
				}
		}
		void remove(Item *item) {// unlinks without deleting the item
				if(n==0) return;
				bool found=false;
				List *list=first;
				do {
						if(list->item==item) { found=true; break; }
						list=list->next;
				}	while (list!=first);
				if(!found) { 
						if (verbose) {
								cerr <<"WARNING: CANNOT REMOVE FROM LIST: ITEM NOT FOUND" << endl;
						}
						return;
				}
				list->item=nullptr;
				if(tag==list) {
						tag=nullptr;
				}
				if(list==current)current=current->next;
				if(list==first)first=first->next;
				list->prev->next=list->next;
				list->next->prev=list->prev;
				delete list;
				list=nullptr;
				if(--n==0){
						current=first=tag=nullptr;
				}
		}
		void erase(Item *item) {// unlinks and deletes the item
				if(n==0) return;
				bool found=false;
				List *list=first;
				do {
						if(list->item==item) { found=true; break; }
						list=list->next;
				}	while (list!=first);
				if(!found) { 
						if (verbose) {
								cerr <<"CANNOT ERASE FROM LIST: ITEM NOT FOUND" << endl;
						}
						return;
				}
				delete list->item;
				list->item=nullptr;
				if(tag==list) {
						tag=nullptr;
				}
				if(list==current)current=current->next;
				if(list==first)first=first->next;
				list->prev->next=list->next;
				list->next->prev=list->prev;
				delete list;
				list=nullptr;
				if(--n==0){
						current=first=tag=nullptr;
				}
		}
		void goFirst() {
				current=first;
		} 
		void goLast() {
				if (first!=nullptr) current=first->prev;
		} 
		void goNext() {
				if(current!=nullptr) current=current->next;
		} 
		void goPrev() {
				if(current!=nullptr) current=current->prev;
		} 
		bool isFirst() {
				return current==first?true:false;
		}
		bool isLast() {
				if (first!=nullptr) return current==first->prev?true:false;
				return nullptr;
		}
		void goTag() {
				if (tag != nullptr) current=tag;
		}
		void setTag() {
				tag=current;
		}
		void untag() {
				tag=nullptr;
		}
		bool isTagged() {
				return tag == nullptr ? false : true;
		}
		void removeTagged() {
				goTag();
				remove();
		}
		void tagFirst() {
				tag = first;
		}
		void tagNext() {
				if(tag!=nullptr) tag=tag->next;
		} 
		void tagPrev() {
				if(tag!=nullptr) tag=tag->prev;
		} 
		bool tagIsFirst() {
				return tag==first?true:false;
		}
		bool tagIsCurrent() {
				return tag==current?true:false;
		}
		Item *getItem() {
				if(current==nullptr)return nullptr;
				return current->item;
		}
		void replaceItem(Item *newitem) {
				if(current==nullptr)return;
				current->item=newitem;
		}
		Item *getTaggedItem() {
				if(tag==nullptr)return nullptr;
				return tag->item;
		}
		void clean() {
				if (n == 0) return;
				current = first->next;
				while(current != first) {
						current = current->next;
						delete current->prev;
				}
				delete first;
				current=first=tag=nullptr;
				n=0;
		}
		void wipe() {
				if (first == nullptr) return;
				current = first;
				do {
						delete current->item;
						current = current->next;
				} while(current != first);
				clean();
		}
		void show(string name) {
				cout << "Collection "<<name<<" of "<<n<<" items" << endl;
				current = first;
				do {
						current->item->show(name);
						current = current->next;
				} while (current != first);
		}
};

struct Node;
class Triangle;

/**
 * Segment cosists of two nodes and a vector
 * in the direction perpendicular to the line
 * connecting the nodes. The direction points
 * inside the already generated mesh and is
 * used to avoid points already in the mesh
 * when selecting a new point for a surface
 * triangle.
 */
class Segment {

		Node *nodes[2],*apex,*nadir; 
		// nodes[0:1] end vertices
		// appex: opposite of the front segment
		// nadir: opposite of the apex (used in mesh optimization)
		Vector E; // normal to edge node[0] - node[1]
	  // pointing out of the mesh
		// node[0], node[1] are the end nodes of the segment
		// node[2] is the opposite node of the mesh triangle:
		// apex node
		Triangle *faces[2]; // points to two adjacent triangular faces
		// used in hole-patching and mesh optimization routines

		int vertex[2]; // index of triangle vertex opposite to this segment 

		public:
		Segment() : E(), apex(nullptr), nadir(nullptr) {
				for (int i=0;i<2;i++) {
						nodes[i] = nullptr;
						vertex[i] = -1;
						faces[i] = nullptr;
				}
		}
		Segment(Node *n0, Node *n1, Node *n2) : 
				E(), nadir(nullptr) {
						for (int i=0;i<2;i++) {
								faces[i]=nullptr;
								vertex[i] = -1;
						}
						nodes[0] = n0;
						nodes[1] = n1;
						apex = n2;
						setNorm();
				}
		~Segment() {
		}
		void clean();
		void setFace0(Triangle *tri, int vert) {
				vertex[0] = vert;
				*faces = tri;
		}
		void setFace1(Triangle *tri, int vert) {
				vertex[1] = vert;
				faces[1] = tri;
		}
		void setFace(int i, int vert, Triangle *tri) {
#ifdef ACERT
				if (i<0 || i>1) return;
#endif
				faces[i] = tri;
				vertex[i] = vert;
		}
		Triangle *getFace0() {
				return *faces;
		}
		Triangle **getFaces() {
				return faces;
		}
		Triangle *getFace(int i) {
				return faces[i];
		}
		void setEndNodes(Node *n0, Node *n1) {
				nodes[0]=n0;
				nodes[1]=n1;
		}
		Vector *getNorm() {
				return &E;
		}

		void setNode(int i, Node *n) {
#ifdef ACERT
				if(i<0||i>1) {
						cerr <<"#! CANNOT GET Node "<<i<<" of a segment"<<endl;
						return nullptr;
				}
#endif
				nodes[i] = n;
		}

		Node *getNode(int i) {
#ifdef ACERT
				if(i<0||i>1) {
						cerr <<"#! CANNOT GET Node "<<i<<" of a segment"<<endl;
						return nullptr;
				}
#endif
				return nodes[i];
		}

		Node *getNode(Point *A);

		Node **getNodes() {
				return nodes;
		}
		Node *getOtherNode(Node *n) {
				for (int i=0; i<2; i++) {
						if (nodes[i] != n) return nodes[i];
				}
		}
		Node *getApex() {
				return apex;
		}
		void setApex(Node *nd) {
				apex = nd;
		}
		Node *getNadir() {
				return nadir;
		}
		void setNadir(Node *nd) {
				nadir = nd;
		}
		REAL length();

		int getVertInd(int i) {
				return vertex[i];
		}

		int setVertInd(int i, int ind) {
				vertex[i] = ind;
		}

		void replace(Triangle *A, Triangle *B, int ind) {
				for (int i=0; i<2; i++) {
						if (faces[i] == A) {
								faces[i] = B;
								vertex[i] = ind;
								return;
						}
				}
				cerr << "ERROR: Failed to replace face pointer\n";
		}

		void swap(Triangle *A, Triangle *B, int ivert) {
				for (int i=0; i<2; i++) {
						if (faces[i] == A) {
								faces[i] = B;
								vertex[i] = ivert;
								return;
						}
						if (faces[i] == B) {
								faces[i] = A;
								vertex[i] = ivert;
								return;
						}
				}
				cerr << "ERROR: Failed to swap face pointer\n";
		}


		void replaceNode(Node *A, Node *B) {
				for (int i=0; i<2; i++) {
						if (nodes[i] == A) {
								nodes[i] = B;
								return;
						}
				}
				cerr << "ERROR: Failed to replace end node in segment\n";
		}

		void replaceCorner(Node *A, Node *B) {
				if (A == apex) {apex = B; return;}
				if (A == nadir) {nadir = B; return;}
				cerr << "ERROR: Failed to replace corner node in segment\n";
		}

		void setNorm();
		Point *getVertex(int i);	
		Node *getCorner(Point *P);
		void show(string name);
};

/**
 * Node contains a point and a list of 
 * segments which are connected to it
 */
struct Node {
		int index;//for file output
		Collection<Segment> segs;
		Point *P = nullptr;
		Node() : P(nullptr), segs() {
				P = nullptr;
		}
		~Node() {
				if (P != nullptr) delete P;
				P = nullptr;
		}
		Node(Point *Q): P(nullptr), segs() {
				P = Q;
		}
		Node(Point *Q, Segment *A, Segment *B) :
				P(nullptr),
				segs()
		{
				P = Q;
				segs.add(A);
				segs.add(B);
		}

		REAL coordinate(int i) { return P->X.get(i); }
		REAL *coordinates() { return P->X.vec();}
		REAL distance(REAL y[]) {
				REAL d=0.0F, *x = P->X.vec();
				for(int i=0;i<DIM;i++) {
						REAL r=y[i]-x[i];
						d+=r*r;
				}
				return (REAL)sqrt(d);
		}
		REAL distance(Node *node) {
				REAL d=0.0F, *x=P->X.vec();
				for(int i=0;i<DIM;i++) {
						REAL r=node->coordinate(i)-x[i];
						d+=r*r;
				}
				return (REAL)sqrt(d);
		}
		REAL distance(Node *node, REAL dist[]) {
				REAL d=0.0F, *x=P->X.vec();
				for(int i=0;i<DIM;i++) {
						REAL r=(REAL)(node->coordinate(i)-x[i]);
						d+=r*r;
						dist[i]=(REAL)r;
				}
				return (REAL)sqrt(d);
		}
		void addSegs(Segment *A, Segment *B) {
				segs.add(A);
				segs.add(B);
		}
		void add(Segment *seg) {
				segs.add(seg);
		}
		int numNeibs() {
				return segs.number();
		}
		void replace(Segment *oldseg, Segment *newseg) {
				if(segs.number()==0) return;
				segs.goFirst();
				do {
						Segment *seg=segs.getItem();
						if(seg==oldseg) {
								segs.replaceItem(newseg);
								return;
						}
						segs.goNext();
				}	while(!segs.isFirst());
				cerr<<"#! CANNOT REPLACE "<<oldseg<<" WITH "<<newseg<<endl;cerr.flush();
				exit(1);
		}
		bool isNeib(Node *node) {
				if(node==this) {
						cerr<<"#! THIS NODE "<<node<<" CANNOT BE NEIGHBOR OF ITSELF"<<endl;cerr.flush();
						if(!silent)cout<<"ERROR: THIS NODE "<<node<<" CANNOT BE NEIGHBOR OF ITSELF"<<endl;
						return false;
				}
				if(segs.number()==0) return false;
				segs.goFirst();
				do {
						Node **neibs=segs.getItem()->getNodes();
						for (int i=0; i<2; i++) {
								if(neibs[i]==node) return true;
						}
						segs.goNext();
				}	while(!segs.isFirst());
				return false;
		}
		void show(string name) {
				REAL *x = P->X.vec();
				cout<<"#"<<name<<" P->i="<<P->i<<", X:";cout.flush();
				for (int i=0; i<DIM; i++)
						cout <<' '<<x[i];
				cout << endl;cout.flush();
		}
		void setSegs(Segment *s0, Segment *s1) {
				segs.add(s0);
				segs.add(s1);
		}
		Collection<Segment> *getSegs() {
				return &segs;
		}
		Segment *getSeg(Node *node) {
				//Return segment containing node 
				if(node==this) {
						cerr<<"#! CANNOT RETURN SEGMENT FOR this NODE"<<endl;
						return nullptr;
				}
				if(segs.number()>0) {
						segs.goFirst();
						do {
								Segment *seg=segs.getItem();
								for(int in=0;in<2;in++) {
										Node *nb=seg->getNode(in);
										if(node==nb) return seg;
								}
								segs.goNext();
						}	while(!segs.isFirst());
				}
				return nullptr;
		}
		Segment *getOtherSeg(Segment *seg) {
				//Return other segment but seg
				if(segs.number()>0) {
						segs.goFirst();
						do {
								Segment *s=segs.getItem();
								if (s != seg) return s;
								segs.goNext();
						}	while(!segs.isFirst());
				}
				return nullptr;
		}
		void remove(Segment *seg) {
				if(segs.number()==0) return;
				segs.goFirst();
				do {
						if(seg==segs.getItem()) {
								segs.remove();
								return;
						}
						segs.goNext();
				}	while(!segs.isFirst());
				cerr<<"#! FAILED TO DELETE NODE SEGMENT: "<<seg<<endl;
		}
};

REAL Segment::length() {
		REAL d=0.0;
		for(int i=0;i<2;i++) {
				REAL r=nodes[1]->coordinate(i) - nodes[0]->coordinate(i);
				d+=r*r;
		}
		return (REAL)sqrt(d);
}

void Segment::setNorm() {
		Point 
				*pA = nodes[0]->P,
				*pB = nodes[1]->P,
				*P = apex->P;
		Vector
				A(pA->X),
				B(pB->X),
				C, D, AB;
		AB.sub(B,A);
		AB.one();
		C.mid(A, B);
		D.sub(C, P->X);
		REAL d = dot(D,AB);
		AB.mul(d);
		E.sub(D,AB); // normal to AB pointing out of meshed triangle
		E.one();
}

Node *Segment::getNode(Point *A) {
		for (int i=0; i<2; i++) {
				Node *node = nodes[i];
				if (A == node->P) return node;
		}
}

Node *Segment::getCorner(Point *A) {
		if (apex->P == A) return apex;
		if (nadir->P == A) return nadir;
		return nullptr;
}

void Segment::show(string name) {
		cout<<"#! Segment "<<name<<endl;
		for (int i=0; i<2; i++) {
				cout << "\tNode "<<i<<": "<<nodes[i]<<":"<<nodes[i]->P;
				nodes[i]->show(" show:");
		}
		cout << " Apex: "<<apex<<":"<<apex->P;
		apex->show(" show:");
		cout << "Nadir: "<<nadir<<":"<<nadir->P;
		nadir->show(" show:");
}

/**
 * Triangle holds pointers to three points
 * which form its vertices
 */
class Triangle {
		const static int nvert=3;
		Point *verts[nvert];
		Segment *segs[nvert]; // used in mesh optimization
	public:

	Triangle() {
		for (int i=0; i<<NV; i++) {
			verts[i] = nullptr;
			segs[i] = nullptr;
		}
	}
	Triangle(Point *A, Point *B, Point *C) {
		verts[0]=A;
		verts[1]=B;
		verts[2]=C;
		for (int i=0; i<<NV; i++) {
			segs[i] = nullptr;
		}
	}

	void setVert(int i, Point *P) {
#ifdef ACERT
		if (i>=0&&i<NV) 
#endif
			verts[i] = P;
	}
	Point *getVert(int i) { return verts[i]; }

	void copy(Triangle *tri) {
		for (int i=0; i<nvert; i++) {
			verts[i] = tri->getVert(i);
			segs[i] = tri->getSeg(i);
		}
	}

	bool replaceVert(Point *P, Point *Q) {
		for (int i=0; i<nvert; i++) {
			if (P == verts[i]) {
				verts[i] = Q;
				return true;
			}
		}
		return false;
	}

	int getInd(Point *P) {
		for (int i=0; i<nvert; i++) {
			if (P == verts[i]) return i;
		}
		return -1;
	}

	void setSeg(int i, Segment *seg) {
#ifdef ACERT
		if (i>=0&&i<NV) 
#endif
			segs[i] = seg;
	}

	Segment *getSeg(int i) { return segs[i]; }

	Segment *getSeg(Node *A, Node *B) {
		for (int i=0; i<nvert; i++) {
			Segment *seg = segs[i];
			int match = 0;
			for (int j=0; j<2; j++) {
				Node *node = seg->getNode(j);
				if (node == A || node == B) match++;
			}
			if (match == 2) return seg;
		}
		return nullptr;
	}

	bool replaceSeg(Segment *P, Segment *Q) {
		for (int i=0; i<nvert; i++) {
			if (P == segs[i]) {
				segs[i] = Q;
				return true;
			}
		}
		return false;
	}

	bool has(Point *point) {
		for (int i=0; i<NV; i++) {
			if (verts[i] == point) return true;
		}
		return false;
	}

	REAL perimeter() {
		REAL p = 0.0;
		for (int i=0; i<3; i++) {
			Vector 
				a = verts[i]->X,
				b = verts[(i+1)%3]->X;
			p += dist(a, b);
		}
		return p;
	}

	void areaVec(Vector &Area) {
#ifdef ACERT
		for (int i=0;i<NV;i++) {
			if (verts[i] == nullptr) return 0.0;
		}
#endif
		Vector
			&A = verts[0]->X,
			&B = verts[1]->X,
			&C = verts[2]->X, 
			AB, AC;
		sub(B, A, AB);
		sub(C, A, AC);
		cross(AB, AC, Area);
		Area.mul(0.5);
	}

	REAL area() {
		Vector Area;
		areaVec(Area);
		REAL a = Area.len();
		return a;
	}

	void surfaceNormalVec(Vector &Norm) {
		areaVec(Norm);
		Norm.one();
	}	

	void show(string name){
		cout<<name<<":"<<endl;
		for (int i=0;i<3;i++) {
			char msg[16];
			sprintf(msg,"Vert %d: ",i);
			if (verts[i]!=nullptr) verts[i]->show(msg);
		}
		for (int i=0;i<3;i++) {
			char msg[16];
			sprintf(msg,"Seg %d:%0x",i,segs[i]);
			if (segs[i]!=nullptr) segs[i]->show(msg);
		}
		cout << endl;
	}
};

Point *Segment::getVertex(int i) {
#ifdef ACERT
	if (i<0 || i>1) return;
#endif
	if (faces[i] != nullptr) {
		return faces[i]->getVert(vertex[i]);
	}
}


