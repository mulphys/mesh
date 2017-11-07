import java.io.*;

class Vector2D {
	float[] x;
	Vector2D() {
		x=new float[2];
		x[0]=0.0F;
		x[1]=0.0F;
	}
	Vector2D(Vector2D a) {
		x=new float[2];
		x[0]=a.component(0);
		x[1]=a.component(1);
	}
	Vector2D(float x0, float x1) {
		x=new float[2];
		x[0]=x0;
		x[1]=x1;
	}
	public void equal(Vector2D A) {
		x[0]=A.component(0);
		x[1]=A.component(1);
	}	
	public void set(float x0, float x1) {
		x[0]=x0;
		x[1]=x1;
	}
	public void set(Node node) {
		for(int i=0;i<2;i++) 
			x[i]=node.coordinate(i);
	}
	public float component(int i) {
		if(i>1||i<0) {
			System.out.println("Vector component "+i+" out of bounds");
			return 0.0F;
		}
		return x[i];
	}
	public void rotate90() {
		float y=x[0];
		x[0]=x[1];
		x[1]=-y;
	} 
	public void scale(float c) {
		for(int i=0;i<2;i++)x[i]*=c;
	}
	public void add(Vector2D A) {
		for(int i=0;i<2;i++) x[i]+=A.component(i);
	}
	public void subtract(Vector2D A) {
		for(int i=0;i<2;i++) x[i]-=A.component(i);
	}
	public float length() {
		return (float)Math.sqrt(x[0]*x[0]+x[1]*x[1]);
	}
	public float distance(Vector2D A) {
		float 
			x0=A.component(0)-x[0],
			x1=A.component(1)-x[1];
		return (float)Math.sqrt(x0*x0+x1*x1);
	}
	public void midpoint(Vector2D A, Vector2D B) {
		equal(A);add(B);scale(0.5F);
	}
	public float scalprod(Vector2D A) {
		return x[0]*A.component(0)+x[1]*A.component(1);
	}
	public float vecprod(Vector2D A) {
		return x[0]*A.component(1)-x[1]*A.component(0);
	}
}
class BoundaryVertex2D {
	float[] x;
	BoundaryVertex2D () {
		x = new float[2];
	}
}
class BoundarySegment2D {
	int[] v;
	BoundarySegment2D () {
		v = new int[2];
	}
}
class Node {
	static final int dim=2;//2D IMPLEMENTATION
	int index,//for file output
		type;
	float x[],y[],area;
	Collection<Node> neibs,dead;
	Collection<Segment> segs;
	Node() {
		type=0;
		x=new float[dim];
		y=new float[dim];
		area=0.0F;
		neibs=new Collection<Node>();//front neighbors
		dead=new Collection<Node>();//previous front neighbors
		segs=new Collection<Segment>();//front segments
	}
	Node(float x0, float x1) {
		this();
		x[0]=x0;
		x[1]=x1;
		for(int i=0;i<dim;i++) y[i]=x[i];
	}
	Node(Vector2D V) {
		this(V.component(0),V.component(1));
	}
	public void type(int type) { this.type=type; }
	public int type() { return type; }
	public float coordinate(int i) { return x[i]; }
	public void coordinate(int i, float value) { x[i]=value; }
	public float[] coordinates() { return x; }
	public void setCoordinates(float[] x) { for(int i=0;i<dim;i++) this.x[i]=x[i]; }
	public void area(float value) { area=value; }
	public void incArea(float value) { area+=value; }
	public float area() { return area; }
	public void newCoordinate(int i, float value) { y[i]=value; }
	public void incNewCoordinate(int i, float value) { y[i]+=value; }
	public float newCoordinate(int i) { return y[i]; }
	public float distance(float y[]) {
		double d=0.0F;
		for(int i=0;i<dim;i++) {
			double r=y[i]-x[i];
			d+=r*r;
		}
		return (float)Math.sqrt(d);
	}
	public float distance(Node node) {
		double d=0.0F;
		for(int i=0;i<dim;i++) {
			double r=node.coordinate(i)-x[i];
			d+=r*r;
		}
		return (float)Math.sqrt(d);
	}
	public float distance(Node node, float[] dist) {
		double d=0.0F;
		for(int i=0;i<dim;i++) {
			double r=(double)(node.coordinate(i)-x[i]);
			d+=r*r;
			dist[i]=(float)r;
		}
		return (float)Math.sqrt(d);
	}
	public void setNeibs(Node n0, Node n1) {
		neibs.add(n0);
		neibs.add(n1);
	}
	public void add(Node node) {
		neibs.add(node);
	}
	public void add(Segment seg) {
		segs.add(seg);
	}
	public int numberOfNeibs() {
		return neibs.number();
	}
	public void replace(Node oldneib,Node newneib) {
//-		for(int i=0;i<2;i++) if(neib[i]==oldneib) neib[i]=newneib;
		if(neibs.number()==0) return;
		neibs.goFirst();
		do {
			Node node=neibs.getItem();
			if(node==oldneib) {
				neibs.replaceItem(newneib);
				return;
			}
			neibs.goNext();
		}	while(!neibs.isFirst());
		System.out.println("#! CAN'T REPLACE "+oldneib+" WITH "+newneib);
		System.exit(1);
	}
	public void replace(Segment oldseg,Segment newseg) {
//-		for(int i=0;i<2;i++) if(segs[i]==oldseg) segs[i]=newseg;
		if(segs.number()==0) return;
		segs.goFirst();
		do {
			Segment seg=segs.getItem();
			if(seg==oldseg) {
				segs.replaceItem(newseg);
				return;
			}
			segs.goNext();
		}	while(!segs.isFirst());
		System.out.println("#! CAN'T REPLACE "+oldseg+" WITH "+newseg);
		System.exit(1);
	}
	public boolean isNeib(Node node) {
//-		for(int i=0;i<2;i++) if(neib[i]==node) return true;
		if(node==this) {
			System.out.println("#! THIS NODE "+node+" CAN'T BE NEIGHBOR OF ITSELF");
			return false;
		}
		if(neibs.number()==0) return false;
		neibs.goFirst();
		do {
			Node neib=neibs.getItem();
			if(neib==node) return true;
			neibs.goNext();
		}	while(!neibs.isFirst());
		return false;
	}
	public boolean wasNeib(Node node) {
		if(this==node) {
			System.out.println("#! NODE "+node+" CAN'T BE COMPARED WITH ITSELF");
			return false;
		}
		if(dead.number()==0) return false;
		dead.goFirst();
		do {
			Node oldneib=dead.getItem();
			if(oldneib==node) {
				return true;
			}
			dead.goNext();
		}	while(!dead.isFirst());
		return false;
	}
	public Node getNeib(Node node) {
		//Get the other neighbor node of this node,
		//presuming that 'node' is a neighbor
		if(node==this) {
			System.out.println("#! THIS NODE "+node+" CAN'T BE NEIGHBOR OF ITSELF");
			return null;
		}
		if(neibs.number()==0) return null;
		neibs.goFirst();
		do {
			Node neib=neibs.getItem();
			if(neib!=node) {return neib;}
			neibs.goNext();
		}	while(!neibs.isFirst());
		System.out.println("WARNING in getNeib: Opposite neighbor node of "+this+" to node="+node+" is not found");
		return null;
	}
	public void show(String tag) {
		System.out.println("#"+tag+": "+this+": x:"+x[0]+" "+x[1]+"; neibs:"+neibs.number()+":");
		if(neibs.number()>0) {
			neibs.goFirst();
			do {
				Node neib=neibs.getItem();
				System.out.println("#\t"+neib);
				neibs.goNext();
			}	while(!neibs.isFirst());
		}
//+		System.out.println("#\t segs: "+segs.number()+":");
//+		if(segs.number()>0) {
//+			segs.goFirst();
//+			do {
//+				Segment seg=segs.getItem();
//+				System.out.println("#\t"+seg);
//+				segs.goNext();
//+			}	while(!segs.isFirst());
//+		}
	}
	public void setSegs(Segment s0, Segment s1) {
//-		segs[0]=s0; segs[1]=s1;
		segs.add(s0);segs.add(s1);
	}
	public Segment getSeg(Node node) {
		//Return segment containing node a
		if(node==this) {
			System.out.println("#! CAN'T RETURN SEGMENT FOR this NODE");
			return null;
		}
//-		for(int is=0;is<2;is++) {
//-			Segment seg=segs[is];
//-			for(int in=0;in<2;in++) {
//-				Node nb=seg.getNode(in);
//-				if(node==nb) return seg;
//-			}
//-		}
		if(segs.number()>0) {
			segs.goFirst();
			do {
				Segment seg=segs.getItem();
				for(int in=0;in<2;in++) {
					Node nb=seg.getNode(in);
					if(node==nb) return seg;
				}
				segs.goNext();
			}	while(!segs.isFirst());
		}
		return null;
	}
	public void delete(Node node) {
		if(neibs.number()==0) return;
		neibs.goFirst();
		do {
			if(node==neibs.getItem()) {
				dead.add(node);
				neibs.delete();
				return;
			}
			neibs.goNext();
		}	while(!neibs.isFirst());
		System.out.println("#! FAILED TO DELETE NODE NEIGHBOR: "+node);
	}
	public void delete(Segment seg) {
		if(segs.number()==0) return;
		segs.goFirst();
		do {
			if(seg==segs.getItem()) {
				segs.delete();
				return;
			}
			segs.goNext();
		}	while(!segs.isFirst());
		System.out.println("#! FAILED TO DELETE NODE SEGMENT: "+seg);
	}
	public void move(float dx, float dy) {
		x[0]+=dx; 
		x[1]+=dy;
	}
	public void move(float[] dx) {
		for(int i=0;i<dim;i++) x[i]+=dx[i]; 
	}
}
class Segment {
	Node[] node;
	Segment() {
		node=new Node[2];
	}
	Segment(Node n0, Node n1) {
		this();
		node[0]=n0;
		node[1]=n1;
	}
	public void setNodes(Node n0, Node n1) {
		node[0]=n0;
		node[1]=n1;
	}
	public Node getNode(int i) {
		if(i<0||i>1) {
			System.out.println("#! CAN'T GET Node "+i+" of a segment");
			return null;
		}
		return node[i];
	}
	public float length() {
		double d=0.0F;
		for(int i=0;i<2;i++) {
			double r=node[1].coordinate(i) - node[0].coordinate(i);
			d+=r*r;
		}
		return (float)Math.sqrt(d);
	}
}
class Triangle {
	final static int nvert=3;
	Node[] verts;
	Triangle(Node n0, Node n1, Node n2) {
		verts=new Node[nvert];
		verts[0]=n0;
		verts[1]=n1;
		verts[2]=n2;
	}
	public Node getVert(int i) { return verts[i]; }
	public float area() {//2D IMPLEMENTATION
		final int dim=2;
//-		Vector2D 
//-			A=new Vector2D(Node[0]),
//-			B=new Vector2D(Node[1]);
//-		return a=0.5F*A.vprod(B);
		float[] ab=new float[dim], ac=new float[dim];
		for(int i=0;i<dim;i++) {
			ab[i]=verts[1].coordinate(i)-verts[0].coordinate(i);
			ac[i]=verts[2].coordinate(i)-verts[0].coordinate(i);
		}
		return 0.5F*Math.abs(ab[0]*ac[1]-ab[1]*ac[0]);
	}
}
public class PFMesh {
	final static int 
		dim=2,
		type_boundary_node=0,
		type_internal_node=1,
		type_control_node =2;
	int 
		mbv, //maximum number of boundary vertexes
		nbv, //actual number of boundary vertexes
		mbs, //maximum number of boundary segments
		nbs; //number of boundary segments
	float cellsize,domainsize;
	public Contour contour;
	public BoundaryVertex2D[] b2Dverts;
	public BoundarySegment2D[] b2Dsegs;
	Collection<Node> 
		points,//DEBUG
		bnodes,//boundary vertex
		mnodes,//mesh nodes
		fnodes;//front nodes
	Collection<Segment> 
		msegs,//mesh boundary segements
		bsegs,//boundary segments
		fsegs;//front segments
	Collection<Triangle> triangles;
	public PFMesh() {
		mbv=0;
		nbv=0;
		nbs=0;
		cellsize=1.0F;
		b2Dverts=null;
		b2Dsegs=null;
		points=bnodes=mnodes=fnodes=null;
		msegs=fsegs=null;
		contour=null;
	}
	public void reset() {
		mbv=0;
		nbv=0;
		nbs=0;
		cellsize=1.0F;
		b2Dverts=null;
		b2Dsegs=null;
		points=mnodes=fnodes=null;
		msegs=fsegs=null;
		contour=null;
	}
	public Collection<Node> Nodes() { return mnodes; }
	public Collection<Segment> Segments() { return msegs; }
	public Collection<Triangle> Triangles() { return triangles; }
	void ReadBoundary(InputStream is) throws IOException {
		StreamTokenizer st = new StreamTokenizer(
			new BufferedReader(new InputStreamReader(is, "UTF-8")));
		st.eolIsSignificant(true);
		st.commentChar('#');
		scan:
		while (true) {
			switch (st.nextToken()) {
				default:
				break;
//-			break scan;
//-			  case StreamTokenizer.TT_EOL:
//-			break;
			  case StreamTokenizer.TT_EOF:
			break scan;
			case StreamTokenizer.TT_WORD:
				if ("nv".equals(st.sval)) {//number of boundary vertexes
					int nbv=0;
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						nbv=(int)st.nval;
						AllocBVerts(nbv);
					}
				} else
				if ("v".equals(st.sval)) {//boundary vertex
					double x = 0, y = 0, z = 0;
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						x = readDouble(st);
						if (st.nextToken() == StreamTokenizer.TT_NUMBER) { 
							y = readDouble(st);
//								if (st.nextToken() == StreamTokenizer.TT_NUMBER)
//									z = readDouble(st);
							addBVert((float) x, (float) y);
						}
					}
				} else  
				if ("nbs".equals(st.sval)) {//initialize boundary segments
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						int nbs = (int) st.nval;
						AllocBSegs(nbs);
					}
				} else
				if ("bs".equals(st.sval)) {//boundary segment
					int v1=0,v2=0;
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						v1 = (int)st.nval;
						if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
							v2 = (int)st.nval;
							addBSeg(v1-1, v2-1);
						}
					}
				}
			}
		}
		is.close();
		if (st.ttype != StreamTokenizer.TT_EOF){
			System.out.println("File Format Exception in "+st.toString());
			System.exit(1);
		}
	}
	void ReadControlPoints(InputStream is) throws IOException {
		StreamTokenizer st = new StreamTokenizer(
			new BufferedReader(new InputStreamReader(is, "UTF-8")));
		st.eolIsSignificant(true);
		st.commentChar('#');
		scan:
		while (true) {
			switch (st.nextToken()) {
				default:
				break;
			  case StreamTokenizer.TT_EOF:
			break scan;
			case StreamTokenizer.TT_WORD:
				if ("np".equals(st.sval)) {//number of boundary vertexes
					int np=0;
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						np=(int)st.nval;
						contour=new Contour(np);
					}
				} else
				if ("nv".equals(st.sval)) {//number of boundary vertexes
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						mbv=(int)st.nval;
						AllocBVerts(mbv);
					}
				} else
				if ("p".equals(st.sval)) {//boundary vertex
					double x = 0, y = 0, z = 0;
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						x = readDouble(st);
						if (st.nextToken() == StreamTokenizer.TT_NUMBER) { 
							y = readDouble(st);
//-							if (st.nextToken() == StreamTokenizer.TT_NUMBER)
//-								z = readDouble(st);
							contour.addPoint((float) x, (float) y);
						}
					}
				}
			}
		}
		is.close();
		if (st.ttype != StreamTokenizer.TT_EOF){
			System.out.println("File Format Exception in "+st.toString());
			System.exit(1);
		}
	}
	public void setControlPoints(InputStream is, float x, float y) throws IOException {
		final float eps=0.15F;
		contour=new Contour();//Cubic bezier
		int np=contour.maxNumberPoints();
		
		contour.addPoint(0.0F,0.0F);
		contour.addPoint(0.0F,1.0F);
		contour.addPoint(x,y);
		StreamTokenizer st = new StreamTokenizer(
			new BufferedReader(new InputStreamReader(is, "UTF-8")));
		st.eolIsSignificant(true);
		st.commentChar('#');
		scan:
		while (true) {
			switch (st.nextToken()) {
				default:
				break;
			  case StreamTokenizer.TT_EOF:
			break scan;
			case StreamTokenizer.TT_WORD:
				if ("nv".equals(st.sval)) {//number of boundary vertexes
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						mbv=(int)st.nval;
						AllocBVerts(mbv);
					}
				}
			}
		}
		is.close();
		if (st.ttype != StreamTokenizer.TT_EOF){
			System.out.println("File Format Exception in "+st.toString());
			System.exit(1);
		}
	}
	public void setControlPointsRandom(InputStream is) throws IOException {
		final float eps=0.15F;
		contour=new Contour();//Cubic bezier
		int np=contour.maxNumberPoints();
		
		contour.addPoint(0.0F,0.0F);
		contour.addPoint(0.0F,1.0F);
		contour.addPoint(1.0F+eps*(float)(1.0-2.0*Math.random()),1.0F+eps*(float)(1.0-2.0*Math.random()));
		StreamTokenizer st = new StreamTokenizer(
			new BufferedReader(new InputStreamReader(is, "UTF-8")));
		st.eolIsSignificant(true);
		st.commentChar('#');
		scan:
		while (true) {
			switch (st.nextToken()) {
				default:
				break;
			  case StreamTokenizer.TT_EOF:
			break scan;
			case StreamTokenizer.TT_WORD:
				if ("nv".equals(st.sval)) {//number of boundary vertexes
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						mbv=(int)st.nval;
						AllocBVerts(mbv);
					}
				}
			}
		}
		is.close();
		if (st.ttype != StreamTokenizer.TT_EOF){
			System.out.println("File Format Exception in "+st.toString());
			System.exit(1);
		}
	}
	public void GenerateBoundary() {
		final int dim=2;
		if(mbv==0){
			System.out.println("Invalid number of boundary vertexes "+nbv);
			return;
		}
		final float dt=1.0F/(float)mbv;
		for(int iv=0;iv<mbv;iv++) {
			float x[]=contour.vertex((float)iv*dt);
			addBVert(x[0], x[1]);
		}
		mbs=mbv;
		AllocBSegs(mbs);
		for(int is=0;is<mbs-1;is++) {
			addBSeg(is, is+1);
		}
		addBSeg(mbs-1,0);
	}
	public void Init() {
		initFront();
		initMesh();
	}
	public void initFront() {
		final float powonesixthoftwo=(float)Math.pow(2.0,1.0/6.0);
		Node[] nodes=new Node[nbv];//create integer index for nodes
		bnodes=new Collection<Node>();
		fnodes=new Collection<Node>();
		for(int inode=0;inode<nbv;inode++) {
			BoundaryVertex2D bv=b2Dverts[inode];
			Node node=new Node(bv.x[0],bv.x[1]);
			node.type(type_boundary_node);
			fnodes.add(node);
			bnodes.add(node);
			nodes[inode]=node;
		}
		bsegs=new Collection<Segment>();
		fsegs=new Collection<Segment>();
		cellsize=0.0F;
		for(int is=0;is<nbs;is++) {
			BoundarySegment2D bs=b2Dsegs[is];
			Node n0=nodes[bs.v[0]],n1=nodes[bs.v[1]];
			Segment seg=new Segment(n0,n1);
			float seglen=seg.length();
			bsegs.add(seg);
			fsegs.add(seg);
			n0.add(seg);
			n0.add(n1);
			n1.add(seg);
			n1.add(n0);
			cellsize+=seglen;
		}
		cellsize/=(float)nbs;
		nodes=null;//delete
	}
	public void initMesh() {
		//copy all front nodes and segments into mesh:
		mnodes=new Collection<Node>();
		fnodes.goFirst();
		do {
			mnodes.add(fnodes.getItem());
			fnodes.goNext();
		}	while(!fnodes.isFirst());
		mnodes.goLast(); mnodes.tag();// remember the last boundary node
		msegs=new Collection<Segment>();
		fsegs.goFirst();
		do {
			msegs.add(fsegs.getItem());
			fsegs.goNext();
		}	while(!fsegs.isFirst());
		triangles=new Collection<Triangle>();
	}
	public Node addNode(Vector2D D, Node a, Node b, Segment s) {
		Node	c=new Node(D);
		c.type(type_internal_node);
		Segment ac=new Segment(a,c);
		a.replace(b,c);
		a.replace(s,ac);
		Segment cb=new Segment(c,b);
		b.replace(a,c);
		b.replace(s,cb);
		c.add(a);
		c.add(b);
		c.add(ac);
		c.add(cb);	
		fnodes.add(c);
		mnodes.add(c);
		fsegs.add(ac);
		msegs.add(ac);///to draw the mesh 
		fsegs.add(cb);
		msegs.add(cb);///to draw the mesh
		return c;
	}
	public void Step() {
		Steps(1);
	}
	public void Triangulate() {
		Steps(Integer.MAX_VALUE);
	}
	public void Steps(int nsteps) {
		//should follow Init()
		final double pi=4.0*(float)Math.atan(1.0);
		final float smallangle=(float)Math.cos(pi/6.0);
/*
 * Triangulation using heuristic propagating front algorithm
 *
            D
   G*       *
           /|\
          / | \
         / F*  \
        /  E*   \
       /    |    \
      /     |     \
     *------*------*
     A      C      B
D is a new candidate node such that ABD is a right triangle: AB=BD=AD
C is the middle point between A and B
E is the center of the cucumscribed circle for triangle ABD, 
  with the radius r=ED
F is the center of the circle with the radius r>ED
G is the node on the front closest to F
ab=length(AB)=length(A.separation(B))
ALGORITHM:
- Find C=(A+B)/2
- Find F:
-- Rotate vector AB by PI/2: AB90=rotate90(AB)
-- Normalize: CD1=AB90/length(AB90); //unit vector in direction CD
-- Find point F:
--- Scale: cd=ab*3^(1/2)/2; ce=cd/3; cf=0.7*cd=ce*scale; (heuristic)
--- F=C+CD1*cf
-- Radius: af=length(AF)-eps (eps is needed to avoid nodes A and B)
- Find all the front nodes inside the circle: F,af
- if (such nodes exist):
-- Out of these nodes, find node G closest to C.
-- Adopt G as the third vertex for the new triangle ABG, that is, D=G 
- else (!exist):
-- Adopt point D as the third vertex for the new triangle ABG

Smooth mesh by applying forces along struts
*/
//-points=new Collection<Node>();///DDD
		Vector2D //initialize all vectors before the loop for better efficiency:
			A=new Vector2D(),
			B=new Vector2D(),
			AB=new Vector2D(),
			AG=new Vector2D(),
			C=new Vector2D(),
			CD=new Vector2D(),
			CD1=new Vector2D(),//unit vector in direction CD
			D=new Vector2D(),
			//-	E=new Vector2D(),
			F=new Vector2D(),
			CF=new Vector2D(),
			G=new Vector2D();
//+		float maxlength=cellsize;//can be made position dependent
//int iseg=0;///DDD
		if(fsegs!=null&&fsegs.number()>0) {
			fsegs.goFirst();
//-			while(true) {
			int istep=0;
			while(istep<nsteps && fsegs.number()>0) {
				Segment seg=fsegs.getItem();
				//Retrieve two vertexes of a segment:
				Node 
					a=seg.getNode(0),//first vertex of a segment
					b=seg.getNode(1);//second vertex
				A.set(a);
				B.set(b);
//points.clean();///DDD
//points.add(a);///DDD
//points.add(b);///DDD
//points.add(newnode);///DDD
//if(++iseg==90)break;///DDD
				AB.equal(B); AB.subtract(A);
				float ab=AB.length(),
					maxlength=1.5F*ab; //heuristic
			// C=(A+B)/2
				C.midpoint(A,B);
			// Rotate vector AB by PI/2: AB90=rotate90(AB)
				CD1.equal(AB); CD1.scale(1.0F/ab); CD1.rotate90();
			// Normalize to a unit vector in direction CD
				float cd=0.5F*(float)Math.sqrt(3.0)*ab,
			// Find point F:
				cf=1.1F*cd;//heuristic: does not have to be 0.7 ?
				//F=C+CF=C+CD1*cf:
				CF.equal(CD1); CF.scale(cf);
				F.equal(C); F.add(CF);
				//Find front node 'neib' (point G) inside circle(F,r=af-eps) 
				// closest to point C:
				float af=A.distance(F),
					r=0.95F*af;//heuristic: to exclude nodes A,B
				Node newnode=null;fnodes.untag();
				if(fnodes.number()>0) {
					fnodes.goFirst();
					float gcmin=Float.MAX_VALUE;
					do {
						newnode=fnodes.getItem();
						G.set(newnode);
						if(G.distance(F)<r) {
							if(!newnode.wasNeib(a)&&!newnode.wasNeib(b)) {
								float gc=G.distance(C);
								if(gc<gcmin) {
									gcmin=gc;fnodes.tag();
								}
							}
						}
						fnodes.goNext();
					}	while(!fnodes.isFirst());
				}//end if(fnodes.number()>0)
				if((newnode=fnodes.getTaggedItem())!=null) {//Assign node to the new triangle vertex:
					//Check if a and node are neighbors
					Segment delsa=null,delsb=null;
					boolean dela=false,delb=false,addnode=false,addseg=false;
					if(a.isNeib(newnode)) {
						delsa=a.getSeg(newnode);
						newnode.delete(delsa);
						newnode.delete(a);
						//delete node a from the front
						dela=true;
					} else {
						//Compute angle B-A-newnode and compare it to BAC:
						C.set(newnode); C.subtract(A); 
						float an=C.length(),
							angleBAN=AB.scalprod(C)/(ab*an);
						C.set(a.getNeib(b)); C.subtract(A);
						float ac=C.length(),
							angleBAC=AB.scalprod(C)/(ab*ac);
						if(angleBAC+smallangle>=angleBAN) {
							fsegs.goNext();
							continue;
						}
						if(!b.isNeib(newnode)) {
						//Compute angle A-B-newnode and compare it to ABC:
							C.set(newnode); C.subtract(B); 
							float bn=C.length(),
								angleABN=-AB.scalprod(C)/(ab*bn);
							C.set(b.getNeib(a));
							C.subtract(B);
							float bc=C.length(),
							//angle ABC:
								angleABC=-AB.scalprod(C)/(ab*bc);
							if(angleABC+smallangle>=angleABN) {
								fsegs.goNext();
								continue;
							}
						}
						float newseglength=a.distance(newnode);
						if(newseglength>maxlength) {//&&newseglength>b.distance(newnode)&&newseglength>maxlength) {//SPLIT newseg 
							G.set(newnode);
							D.midpoint(A,G);
							addnode=true;
							newnode=addNode(D,a,b,seg);
						} else {// ADD newseg
							Segment newseg=new Segment(a,newnode);
							fsegs.add(newseg);
							msegs.add(newseg);// for mesh output
							a.replace(b,newnode);
							a.replace(seg,newseg);
							newnode.add(a);
							newnode.add(newseg);
							addseg=true;
						}
					}
					if(!addnode) {
//?					if(addseg) {
						if(b.isNeib(newnode)) {
							delsb=b.getSeg(newnode);
							newnode.delete(delsb);
							newnode.delete(b);
							delb=true;
							//delete newnode b from the front
						} else {
							float newseglength=b.distance(newnode);
							if(addseg=false&&newseglength>maxlength) {//&&newseglength>b.distance(newnode)&&newseglength>maxlength) {//SPLIT newseg 
								G.set(newnode);
								D.midpoint(B,G);
								addnode=true;
								newnode=addNode(D,a,b,seg);
							} else {
								Segment newseg=new Segment(newnode,b);
								fsegs.add(newseg);
								msegs.add(newseg);// for mesh output
								b.replace(a,newnode);
								b.replace(seg,newseg);
								newnode.add(b);
								newnode.add(newseg);
							}
						}
						int ndeleted=0;
						if(dela) {
							a.delete(delsa);
							a.delete(newnode);
							a.delete(b);
	//						b.delete(a);
							fsegs.delete(delsa);
							if(a.numberOfNeibs()==0)fnodes.delete(a);
							ndeleted++;
						}
						if(delb) {
							b.delete(delsb);
							b.delete(newnode);
							b.delete(a);
	//						a.delete(b);
							fsegs.delete(delsb);
							if(b.numberOfNeibs()==0)fnodes.delete(b);
							ndeleted++;
						}
						if(ndeleted==2) {//delete newnode from the front
							if(newnode.numberOfNeibs()==0)fnodes.delete(newnode);
						}
//?					} //if addseg
					} //end if(!addnode)
				} else {//Create new triangle vertex: D
					CD.equal(CD1);CD.scale(cd);
					D.equal(C);D.add(CD);
					//Compute angle A-B-newnode and compare it to ABC:
					C.equal(D); C.subtract(A); 
					float an=C.length(),
						angleBAN=AB.scalprod(C)/(ab*an);
					C.set(a.getNeib(b));
					C.subtract(A);
					float ac=C.length(),
						angleBAC=AB.scalprod(C)/(ab*ac);
					if(angleBAC+smallangle>=angleBAN) {
						fsegs.goNext(); 
						continue;
					}
					C.equal(D); C.subtract(B); 
					float bn=C.length(),
						angleABN=-AB.scalprod(C)/(ab*bn);
					C.set(b.getNeib(a));
					C.subtract(B);
					float bc=C.length(),
						angleABC=-AB.scalprod(C)/(ab*bc);
					if(angleABC+smallangle>=angleABN) {
						fsegs.goNext(); 
						continue;
					}
					newnode=addNode(D,a,b,seg);
				}
				triangles.add(new Triangle(a,b,newnode));
				fsegs.delete();
//-				if(fsegs.number()==0) break;
//-				fsegs.goNext();///DDD
//+				fsegs.goLast();///DDD
				istep++;
			}// end while
		}//end if(fsegs.number()>0)
	}
	public void Relax() {
		final int niter=10; //number of iterations
		if(triangles.number()<=0||mnodes.number()<=0) return;
		for(int iter=0;iter<niter;iter++) {
		mnodes.goFirst();
		do {
			Node node=mnodes.getItem();
			for(int i=0;i<dim;i++)node.newCoordinate(i,0.0F);
			node.area(0.0F);
			mnodes.goNext();
		}	while(!mnodes.isFirst());
		int nvert=Triangle.nvert;
		triangles.goFirst();
		do {
			Triangle tri=triangles.getItem();
			float area=tri.area(),
				xc[]=new float[dim];
			for(int i=0;i<dim;i++) {
				float sum=0.0F;
				for(int j=0;j<nvert;j++) {
					sum+=tri.getVert(j).coordinate(i);
				}
				xc[i]=sum/(float)nvert;
			}
			for(int iv=0;iv<nvert;iv++) {
				Node vert=tri.getVert(iv);
				for(int i=0;i<dim;i++)vert.incNewCoordinate(i,xc[i]*area);
				vert.incArea(area);
			}
			triangles.goNext();
		}	while(!triangles.isFirst());
		mnodes.goFirst();
		do {
			Node node=mnodes.getItem();
			if(node.type()==type_internal_node) {
				float areai=1.0F/node.area();
				for(int i=0;i<dim;i++) {
					node.coordinate(i,node.newCoordinate(i)*areai);
				}
			}
			mnodes.goNext();
		}	while(!mnodes.isFirst());
		}
	}
	public void Output(OutputStream outstream) throws IOException {
		boolean output_front=false;
		Writer out=null;
		try {// Set up output stream
			out = new BufferedWriter(new OutputStreamWriter(outstream, "UTF8"));
		} catch (Throwable t) {
			t.printStackTrace();
		}
/*DDD
		out.write("nbv\t"+nbv+"\n");
		for(int ibv=0;ibv<nbv;ibv++) {
			BoundaryVertex2D bv=b2Dverts[ibv];
//+			out.write("bv\t"+bv.x[0]+"\t"+bv.x[1]+"\n");
			out.write("bv\t"+bv.x[0]+"\t"+bv.x[1]+"\t0.0\n");///DDD
		}
		out.write("nbs\t"+nbs+"\n");
		for(int ibs=0;ibs<nbs;ibs++) {
			BoundarySegment2D bs=b2Dsegs[ibs];
			out.write("bs\t"+bs.v[0]+"\t"+bs.v[1]+"\n");
		}
*///DDD
		out.write("#Debugging Points: np=number of points\n");
		out.write("np\t"+points.number()+"\n");
		if(points.number()>0) {
			points.goFirst();
			do {
				Node node=points.getItem();
				out.write("p\t"+node.coordinate(0)+"\t"+node.coordinate(1)+"\t0.0\n");
				///node.index=inode;//set index for further file output
				points.goNext();
			}	while(!points.isFirst());
		}
		if(output_front) {// FOR DEBUGGING
			out.write("#Front nodes: fn=number of front nodes\n");
			out.write("fn\t"+fnodes.number()+"\n");
			if(fnodes.number()>0) {
				int inode=0;
				fnodes.goFirst();
				do {
					Node node=fnodes.getItem();
					out.write("n\t"+node.coordinate(0)+"\t"+node.coordinate(1)+"\t0.0\n");
					node.index=inode;//set index for further file output
					fnodes.goNext();inode++;
				}	while(!fnodes.isFirst());
			}
			out.write("#Front segments: nfs=number of front segments\n");
			out.write("nfs\t"+fsegs.number()+"\n");
			if(fsegs.number()>0) {
				fsegs.goFirst();
				int iseg=0;
				do {
					Segment seg=fsegs.getItem();
					int is1=seg.getNode(0).index+1,is2=seg.getNode(1).index+1;
					out.write("s\t"+is1+"\t"+is2+"\n");
					fsegs.goNext();iseg++;
				}	while(!fsegs.isFirst());
			}
		}// END DEBUGGING
		out.write("#Mesh nodes: mn=number of mesh nodes\n");
		out.write("mn\t"+mnodes.number()+"\n");
		if(mnodes.number()>0) {
			int inode=0;
			mnodes.goFirst();
			do {
				Node node=mnodes.getItem();
				out.write("n\t"+node.coordinate(0)+"\t"+node.coordinate(1)+"\t0.0\n");
				node.index=inode;//set index for further file output
				mnodes.goNext();inode++;
			}	while(!mnodes.isFirst());
		}
		out.write("#Mesh segments: ns=number of mesh segments\n");
		out.write("ns\t"+msegs.number()+"\n");
		if(msegs.number()>0) {
			msegs.goFirst();
			int iseg=0;
			do {
				Segment seg=msegs.getItem();
				int is1=seg.getNode(0).index+1,is2=seg.getNode(1).index+1;
				out.write("s\t"+is1+"\t"+is2+"\n");
				msegs.goNext();iseg++;
			}	while(!msegs.isFirst());
		}
		out.write("#Mesh triangles: nt=number of mesh triangles\n");
		out.write("nt\t"+triangles.number()+"\n");
		if(triangles.number()>0) {
			triangles.goFirst();
			do {
				Triangle t=triangles.getItem();
				out.write("t");
				for(int iv=0;iv<Triangle.nvert;iv++) {
					int index=t.getVert(iv).index+1;
					out.write("\t"+index);
				}
				out.write("\n");
				triangles.goNext();
			}	while(!triangles.isFirst());
		}
		out.close();
	}
	public void AllocBVerts (int nbv) {
		this.mbv=nbv;
		b2Dverts=new BoundaryVertex2D[mbv];
	}
	public void AllocBSegs (int mbs) {
		this.mbs=mbs;
		b2Dsegs=new BoundarySegment2D[mbs];
	}
	public void addBVert(float x, float y) {
		if(nbv==mbv){
			System.out.println("BND-VERTEX NO. "+nbv+" EXCEEDS MAX "+mbv);
			return;
		}
		b2Dverts[nbv]=new BoundaryVertex2D();
		b2Dverts[nbv].x[0]=x;
		b2Dverts[nbv].x[1]=y;
		nbv++;
	}
	public void addBSeg(int v1, int v2) {
		if(nbs==mbs){
			System.out.println("BND-SEGMENT NO. "+nbs+" EXCEEDS MAX. "+mbs);
			return;
		}
		b2Dsegs[nbs]=new BoundarySegment2D();
		b2Dsegs[nbs].v[0]=v1;
		b2Dsegs[nbs].v[1]=v2;
		nbs++;
	}
	static double readDouble(StreamTokenizer st) throws IOException {
		//Reading scientific notation: X.XeX
		double num = st.nval;
		int exp = 0;
		st.ordinaryChars('\0', ' ');
		st.nextToken();
		st.whitespaceChars('\0', ' ');
		if (st.ttype == StreamTokenizer.TT_WORD &&
			Character.toUpperCase(st.sval.charAt(0)) == 'E') {
				try {
					exp = Integer.parseInt(st.sval.substring(1));
				} catch (NumberFormatException e) {
				st.pushBack();
			}
		} else if (!('\0' <= st.ttype && st.ttype <= ' ')) {
			st.pushBack();
		}
		return num * Math.pow(10, exp);
	}
}
