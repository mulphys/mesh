#define INT	int
#define REAL double
#define BIG	1.e30
#define SMALL	1.e-30
#define DIM 3
#define NV 3 

#define MIN_ASPECT_RATIO 0.005
#define ZIP_PROXIMITY 0.2
#define RANGE	2.5

//TODO: move to class

using namespace std;

extern bool silent, verbose;


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
		if (X==NULL) return;
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
		if (X==NULL) return;
		for (int i=0; i<DIM; i++)
			V[i] = X[i];
	}

	inline REAL get(int i) {
		if (i>=0 && i<DIM)
			return V[i];
	}

	inline void set(Vector *X) {
		if (X==NULL) return;
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
		if (X==NULL) return;
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
		V[3] = a[0]*b[1] - a[1]*b[0];
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
	c[3] = a[0]*b[1] - a[1]*b[0];
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
//if (A==NULL || B==NULL || C == NULL || D == NULL || E == NULL) return false;
if (clockwise(A,C,D,E) == clockwise(B,C,D,E)) return false;
if (clockwise(A,B,C,E) == clockwise(A,B,D,E)) return false;
return true;
}


struct Point {
	REAL range = 0.0;
	int i = 0;
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
		return (dot(F, Q->F)+0.25) / d;
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
	if (A == NULL || B == NULL || C == NULL) return 0.0;
	AB.set(B->X);
	AC.set(C->X);
	AB.sub(A->X);
	AC.sub(A->X);
	cross(AB, AC, Area);
	REAL a = 0.5*Area.len();
	return a;
}


template <class Item>
class Collection {
	struct List {
		Item *item = NULL;
		List *next = NULL, *prev = NULL;
		List() :
			item(NULL),
			next(NULL),
			prev(NULL)
		{
			item = NULL;
			next = prev = NULL;
		}
		~List() {
			item = NULL; 
			next = prev = NULL;
		}
		void set(Item *item) {
			this.item = item;
		}
	};
	int n;
	List *current=NULL, *first=NULL, *tag=NULL;
	public:
	Collection() :
		n(0),
		current(NULL),
		first(NULL),
		tag(NULL)
	{
		n=0;
		first=current=tag=NULL;
	}
	~Collection() {
		if (first != NULL) {
			current = first->next;
			if (current != first)
			do {
				tag = current;
				current = current->next;
				delete tag;
			} while (current != first);
			delete current;
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
		if(current==NULL) {
			cerr << "CANNOT DELETE VERTEX" << endl;
			exit(1);
		}
		current->item=NULL;
		if(current==tag) {
			tag=NULL;
		}
		if(current==first)first=first->next;
		current->prev->next=current->next;
		current->next->prev=current->prev;
		List *dead = current;
		current=current->next;
		delete dead;
		if(--n==0){
			current=first=tag=NULL;
		}
	}
	void erase() {// unlinks and deletes the item
		if(n==0) return;
		if(current==NULL) {
			cerr << "CANNOT DELETE VERTEX" << endl;
			exit(1);
		}
		delete current->item;
		current->item=NULL;
		if(current==tag) {
			tag=NULL;
		}
		if(current==first)first=first->next;
		current->prev->next=current->next;
		current->next->prev=current->prev;
		List *dead = current;
		current=current->next;
		delete dead;
		if(--n==0){
			current=first=tag=NULL;
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
		list->item=NULL;
		if(tag==list) {
			tag=NULL;
		}
		if(list==current)current=current->next;
		if(list==first)first=first->next;
		list->prev->next=list->next;
		list->next->prev=list->prev;
		delete list;
		list=NULL;
		if(--n==0){
//?			delete first;
			current=first=tag=NULL;
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
		list->item=NULL;
		if(tag==list) {
			tag=NULL;
		}
		if(list==current)current=current->next;
		if(list==first)first=first->next;
		list->prev->next=list->next;
		list->next->prev=list->prev;
		delete list;
		list=NULL;
		if(--n==0){
//?			delete first;
			current=first=tag=NULL;
		}
	}
	void goFirst() {
		current=first;
	} 
	void goLast() {
		if (first!=NULL) current=first->prev;
	} 
	void goNext() {
		if(current!=NULL) current=current->next;
	} 
	void goPrev() {
		if(current!=NULL) current=current->prev;
	} 
	bool isFirst() {
		return current==first?true:false;
	}
	bool isLast() {
		if (first!=NULL) return current==first->prev?true:false;
		return NULL;
	}
	void goTag() {
		if (tag != NULL) current=tag;
	}
	void setTag() {
		tag=current;
	}
	void untag() {
		tag=NULL;
	}
	bool isTagged() {
		return tag == NULL ? false : true;
	}
	void removeTagged() {
		goTag();
		remove();
	}
	void tagFirst() {
		tag = first;
	}
	void tagNext() {
		if(tag!=NULL) tag=tag->next;
	} 
	void tagPrev() {
		if(tag!=NULL) tag=tag->prev;
	} 
	bool tagIsFirst() {
		return tag==first?true:false;
	}
	bool tagIsCurrent() {
		return tag==current?true:false;
	}
	Item *getItem() {
		if(current==NULL)return NULL;
		return current->item;
	}
	void replaceItem(Item *newitem) {
		if(current==NULL)return;
		current->item=newitem;
	}
	Item *getTaggedItem() {
		if(tag==NULL)return NULL;
		return tag->item;
	}
	void clean() {
		while (number()>0) remove();
		current=first=tag=NULL;
		n=0;
	}
	void wipe() {
		if (first == NULL) return;
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

class Segment {

	Node *node[3]; 
	Vector E; // normal to edge node[0] - node[1]
	          // pointing out of the mesh
	// node[0], node[1] are the end nodes of the segment
	// node[2] is the opposite node of the mesh triangle:
	// apex node
//-	bool boundary; // is this a boundary segment?
	Triangle *triangle; // points to the triangle the segment belongs to
	int vertex; // index of triangle vertex opposite to this segment 

	public:
	Segment() : E(), triangle(NULL), vertex(-1) {
		for (int i=0;i<DIM;i++) node[i]=NULL;
		triangle = NULL;
//-		boundary = false;
		vertex = -1;
	}
	Segment(Node *n0, Node *n1, Node *n2) : E(), vertex(-1), triangle(NULL) {
		triangle = NULL;
//-		boundary = false;
		vertex = -1;
		node[0]=n0;
		node[1]=n1;
		node[2]=n2;
		setNorm();
	}
	~Segment() {
	}
//-	void setBoundary(bool is_boundary) {
//-		boundary = is_boundary;
//-	}
//-	bool isBoundary() {
//-		return boundary;
//-	}
	void setTriangle(Triangle *tri, int vert) {
		triangle = tri;
		vertex = vert;
	}
	Triangle *getTriangle() {
		return triangle;
	}

	void setEndNodes(Node *n0, Node *n1) {
		node[0]=n0;
		node[1]=n1;
	}
	Vector *getNorm() {
		return &E;
	}

	void setNode(int i, Node *n) {
#ifdef ACERT
		if(i<0||i>2) {
			cerr <<"#! CANNOT GET Node "<<i<<" of a segment"<<endl;
			return NULL;
		}
#endif
		node[i] = n;
	}

	Node *getNode(int i) {
#ifdef ACERT
		if(i<0||i>2) {
			cerr <<"#! CANNOT GET Node "<<i<<" of a segment"<<endl;
			return NULL;
		}
#endif
		return node[i];
	}
	Node **getNodes() {
		return node;
	}
	Node *getOtherNode(Node *n) {
		for (int i=0; i<2; i++) {
			if (node[i] != n) return node[i];
		}
	}
	Node *getApex() {
		return node[2];
	}
	void setApex(Node *nd) {
		node[2] = nd;
	}
	REAL length();
	void show(string name) {
		cout<<"#! Segment "<<name<<": "<<node[0]<<", "<<node[1]<<", "<<node[2]<<endl;
	}

	void setNorm();

	Point *getVertex();	
};


struct Node {
	int index;//for file output
	Collection<Segment> segs;
	Point *P = NULL;
	Node() : P(NULL),segs() {
		P = NULL;
	}
	~Node() {
	}
	Node(Point *Q):segs() {
		P = Q;
	}
	Node(Point *Q, Segment *A, Segment *B) {
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
			return NULL;
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
		return NULL;
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
		return NULL;
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
		REAL r=node[1]->coordinate(i) - node[0]->coordinate(i);
		d+=r*r;
	}
	return (REAL)sqrt(d);
}

void Segment::setNorm() {
		Point 
			*pA = node[0]->P,
			*pB = node[1]->P,
			*P = node[2]->P;
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

class Triangle {
	const static int nvert=3;
	Point *verts[nvert];
	public:
	Triangle(Point *A, Point *B, Point *C) {
		verts[0]=A;
		verts[1]=B;
		verts[2]=C;
	}

	void setVert(int i, Point *P) {
#ifdef ACERT
		if (i>=0&&i<NV) 
#endif
			verts[i] = P;
	}
	Point *getVert(int i) { return verts[i]; }

	bool replaceVert(Point *P, Point *Q) {
		for (int i=0; i<NV; i++) {
			if (P == verts[i]) {
				verts[i] = Q;
				return true;
			}
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
	}

	REAL area() {
		for (int i=0;i<NV;i++) {
			if (verts[i] == NULL) return 0.0;
		}
		Vector
			&A = verts[0]->X,
			&B = verts[1]->X,
			&C = verts[2]->X, 
			AB, AC, Area;
		sub(B, A, AB);
		sub(C, A, AC);
		cross(AB, AC, Area);
		REAL a = 0.5*Area.len();
		return a;
	}
	void show(string name){
		cout<<name<<":";
		for (int i=0;i<3;i++) {
			cout << ' ' << verts[i];
			if (verts[i]!=NULL) cout << ':' << verts[i]->i;
		}
		cout << endl;
	}
};

Point *Segment::getVertex() {
	return triangle->getVert(vertex);
}


