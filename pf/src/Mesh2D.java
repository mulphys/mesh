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
class Potential {
	static float maxpotential=Float.MAX_VALUE,///10.0F,
		rmin=0.5F, rmax=1.5F;
	float strength, radius;
	Potential() {
		strength=radius=0.0F;
	}
	Potential(float strength, float radius) {
		this.strength=strength;
		this.radius=radius;
	}
	public float potential(float distance) {//Lennart-Jones
		if(distance==0) {
			System.out.println("Potential at distance="+distance+" is undefined");
			System.exit(1);
		} else 
		if(strength<Float.MIN_VALUE||radius<Float.MIN_VALUE) {
			System.out.println("Potential strength="+strength+" or radius="+radius+" are not set");
			System.exit(1);
		}
		double r=(double)radius/(double)distance;
		if(r<rmin) return 0.0F;//heuristic
		if(r>rmax) r=rmax;//heuristic
		float p=strength*(float)(Math.pow(r,12.0) - Math.pow(r,6));
		if(Math.abs(p)>maxpotential)p=p<0.0F?-maxpotential:maxpotential;
		return p;
	}
}
class Node {
	static final int dim=2;//2D IMPLEMENTATION
	int index,//for file output
		type;
	float x[],y[],v[],area,energy;
	Potential potential;
	Collection<Node> neibs,dead;
	Collection<Segment> segs;
	Node() {
		type=0;
		x=new float[dim];
		y=new float[dim];
		v=new float[dim];
		area=energy=0.0F;
		potential=new Potential();
		neibs=new Collection<Node>();//front neighbors
		dead=new Collection<Node>();//front neighbors
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
	public float energy() { return energy; }
	public void energy(float energy) { this.energy=energy; }
	public float radius() { return potential.radius; }
	public void radius(float r) { potential.radius=r; }
	public float strength() { return potential.strength; }
	public void strength(float value) { potential.strength=value; } 
	public float potential(float dist) { return potential.potential(dist); }
	public float velocity(int i) {
		if(i>dim||i<0) {
			System.out.println("Velocity index "+i+" out of bounds");
			System.exit(1);
		}
		return v[i];
	}
	public void velocity(int i, float value) {
		if(i>dim||i<0) {
			System.out.println("Velocity index "+i+" out of bounds");
			System.exit(1);
		}
		v[i]=value;
	}
	public void setZeroVelocity() { 
		for(int i=0;i<dim;i++) v[i]=0.0F;
	}
	public void setVelocity(float[] velocity) {
		for(int i=0;i<dim;i++) v[i]=velocity[i]; 
	}
	public void incVelocity(float[] dv) {
		for(int i=0;i<dim;i++) v[i]+=dv[i]; 
	}
	public void getVelocity(float[] velocity) {
		for(int i=0;i<dim;i++) velocity[i]=v[i]; 
	}
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
		System.out.println("#\t segs: "+segs.number()+":");
		if(segs.number()>0) {
			segs.goFirst();
			do {
				Segment seg=segs.getItem();
				System.out.println("#\t"+seg);
				segs.goNext();
			}	while(!segs.isFirst());
		}
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
public class Mesh2D {
	final static int 
		dim=2,
		model_monte_carlo=0,
		model_propagating_front=1,
		model_force_field=2,
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
	ForceField forcefield;
	MonteCarlo montecarlo;
	public Mesh2D() {
		mbv=0;
		nbv=0;
		nbs=0;
		cellsize=1.0F;
		b2Dverts=null;
		b2Dsegs=null;
		points=bnodes=mnodes=fnodes=null;
		msegs=fsegs=null;
		contour=null;
		forcefield=null;
		montecarlo=null;
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
	public void setControlPoints(InputStream is) throws IOException {
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
			node.radius(0.0F);
			node.strength(1.0F);// heuristic
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
			n0.radius(0.5F*(n0.radius()+0.5F*seglen/powonesixthoftwo));
			n1.add(seg);
			n1.add(n0);
			n1.radius(0.5F*(n1.radius()+0.5F*seglen/powonesixthoftwo));
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
	public void init(int imodel) {
		switch(imodel) {
		case model_force_field:
			if(forcefield==null) {forcefield=new ForceField(); } 
			else forcefield.init();
			break;
		case model_monte_carlo:
			if(montecarlo==null) {montecarlo=new MonteCarlo(); } 
			else montecarlo.init();
			break;
		}
	}
	public void insert(int imodel) {
		switch(imodel) {
		case model_force_field:
			if(forcefield==null) {forcefield=new ForceField(); forcefield.insert(); } 
			else if(forcefield.done) 
				System.out.println("Force-Field Mesh generated");
			else
				forcefield.insert();
			break;
		case model_monte_carlo:
			if(montecarlo==null) {montecarlo=new MonteCarlo(); montecarlo.insert(); } 
			else if(montecarlo.done) 
				System.out.println("Monte-Carlo Mesh generated");
			else
				montecarlo.insert();
			break;
		}
	}
	public void Triangulate(int imodel) {
		switch(imodel) {
			case model_force_field: 
				if(forcefield==null) {forcefield=new ForceField(); forcefield.triangulate(); } 
				else if(forcefield.done) 
					System.out.println("Force-Field Mesh generated");
				else
					forcefield.triangulate();
				break;
			case model_propagating_front: TriangulatePF(); break;
			case model_monte_carlo: 
				if(forcefield==null) {montecarlo=new MonteCarlo(); montecarlo.triangulate(); } 
				else if(montecarlo.done) 
					System.out.println("Monte-Carlo Mesh generated");
				else
					montecarlo.triangulate();
				break;
		}
	}
	class MonteCarlo {
		boolean done;
		float distav,distmax,time,totenergy,energy,minenergy,
			xmin[], xmax[];
		MonteCarlo() {	
			init();
//-			Potential.rmin=0.3F;
//-			Potential.rmax=3.0F;
//-			Potential.maxpotential=0.1F*(Float.MAX_VALUE-100.0F);//heuristic
			Potential reference_potential=new Potential(1.0F,1.0F);
			minenergy=reference_potential.potential(0.6F);//heuristic
			if(fnodes==null||fnodes.number()<=0||fsegs==null||fsegs.number()==0) {
				System.out.println("Boundary is empty: meshing aborted");
				return;
			}
			done=false;
			xmin=new float[dim];
			xmax=new float[dim];
			float[] xcenter=new float[dim];
			for(int i=0;i<dim;i++) {
				xcenter[i]=0.0F;
				xmin[i]=Float.MAX_VALUE;
				xmax[i]=-Float.MAX_VALUE;
			}
			fnodes.goFirst();
			do {
				Node node=fnodes.getItem();
				node.radius(0.0F); node.strength(1.0F);
				node.setZeroVelocity();
				for(int i=0;i<dim;i++) {
					float x=node.coordinate(i);
					xcenter[i]+=x;
					if(x<xmin[i]) xmin[i]=x;
					if(x>xmax[i]) xmax[i]=x;
				}
				fnodes.goNext();
			}	while(!fnodes.isFirst());
			domainsize=0.0F;
			int nnodes=fnodes.number();//initial number of nodes
			for(int i=0;i<dim;i++) {
				float x=xmax[i]-xmin[i];
				xcenter[i]/=(float)nnodes;
				domainsize=x*x;
			}
			domainsize=(float)Math.sqrt(domainsize);
			for(int i=0;i<dim;i++) {
				xmin[i]-=0.2*domainsize;
				xmax[i]+=0.2*domainsize;
			}
			//Compute the average inter-node distance along the boundary:
			distav=0.0F;
			distmax=0.0F;
			fsegs.goFirst();
			do {
				Segment seg=fsegs.getItem();
				float seglen=seg.length(),
					r=0.5F*seglen;
				distav+=seglen;
				Node 
					a=seg.getNode(0),//first vertex of a segment
					b=seg.getNode(1);//second vertex
				float 
					ra=a.radius(),
					rb=b.radius();
				if(seglen>distmax)distmax=seglen;
				a.radius(1.0F*(ra+r));///0.5F*
				b.radius(1.0F*(rb+r));///0.5F*
				fsegs.goNext();
			}	while(!fsegs.isFirst());
			distav/=(float)fsegs.number();
//-			dx=0.01F*distav;//heuristics
			energy=0.0F;//average equilibrium interaction energy per node
		}
		void init() {
			Potential.rmin=0.3F;
			Potential.rmax=3.0F;
			Potential.maxpotential=0.1F*(Float.MAX_VALUE-100.0F);//heuristic
		}
		void insert() {// Insert randomly:
System.out.println("insert nodes="+mnodes.number()+", energy="+energy);System.out.flush();///DDD
//+			bnodes.goFirst();
//+			do {
				Node
					bnode=bnodes.getItem(), 
					node=new Node();
				node.type(type_internal_node);
				float radius=bnode.radius();
				for(int i=0;i<dim;i++) {
					float x=bnode.coordinate(i)+radius*(float)(1.0F-2.0F*Math.random());
					node.coordinate(i,x);
				}
//-				node.radius(distav); 
				node.radius(radius);
				float energy=energy(node.coordinates());
				node.energy(energy);
				node.strength(1.0F);//heuristic
				totenergy+=energy;
//-			mnodes.goTag(); mnodes.goNext(); mnodes.insert(node);
				mnodes.add(node);
				bnodes.goNext();
//+			}	while(!bnodes.isFirst());
		}
		void relax(int maxiter) {
			if(mnodes==null||mnodes.number()<=0) {
				System.out.println("Relax: node set is empty - operation aborted");
				return;
			}
			//Reach equilibrium
int nnodes=0;///DDD
			float 
				x[]=new float[dim];
			for(int iter=0;iter<maxiter;iter++) {
nnodes=0;
				mnodes.goTag(); mnodes.goNext();//first internal node
				do {
					Node node=mnodes.getItem();
					float dx=0.9F*node.radius();//heuristic
					float oldenergy=node.energy();
					do {//test new coordinates
						for(int i=0;i<dim;i++) 
							x[i]=node.coordinate(i)+dx*(float)(1.0F-2.0F*Math.random());
					}	while(!inside(x));
					float newenergy=energy(x);
					if(newenergy<oldenergy) {
						node.setCoordinates(x);
						node.energy(newenergy);
						totenergy+=newenergy-oldenergy;
					}
nnodes++;
					mnodes.goNext();
				}	while(!mnodes.isFirst());
			}//END for(iter=0..maxiter)
			if(nnodes==0) {
				System.out.println("Relax: No internal nodes found");
				return;
			}
			energy=totenergy/(float)mnodes.number();
			if(energy>minenergy) done=true; else done=false;
System.out.println("Realx: maxiter="+maxiter+", nnodes="+nnodes+", minenergy="+minenergy+", energy="+energy+", totenergy="+totenergy);///DDD
		}
		boolean inside(float[] x) {
			for(int i=0;i<dim;i++) 
				if(x[i]<xmin[i]||x[i]>xmax[i]) return false;
			return true;
		}
		float energy(float x[]) {
			if(mnodes==null||mnodes.number()==0) {
				System.out.println("setVelocity: nodes set empty");
				System.exit(1);
			}
			float energy=0.0F;
			mnodes.tag(1);
			do {
				Node node=mnodes.getItem();
				//-neibradius=mnodes.getItem().radius(),
				energy+=node.potential(node.distance(x));
				mnodes.goNext(); 
			}	while(!mnodes.isTag(1));//end for
			return energy;
		}
		float energy() {
			if(mnodes==null||mnodes.number()==0) {
				System.out.println("setVelocity: nodes set empty");
				System.exit(1);
			}
			mnodes.tag(1);
			Node node=mnodes.getItem();
			float	noderadius=node.radius(),
				neibradius=distav,
				neibdist=domainsize,
				energy=0.0F,
				distance=0.0F;
			int 
				neibtype=-1,//undefined
				nnodes=0;
			for (
				mnodes.goNext(); 
				!mnodes.isTag(1);
				mnodes.goNext()
			) {
				Node next=mnodes.getItem();
				//-neibradius=mnodes.getItem().radius(),
				int nextype=next.type();
				float	nextradius=next.radius();
				while((distance=next.distance(node))<0.1*nextradius) {//heuristic
					//Randomize vector:
					for(int i=0;i<dim;i++) 
						node.coordinate(i,0.2F*noderadius*(float)(1.0-2.0*Math.random()));//heuristic
				}
				energy+=next.potential(distance);
				if(distance<neibdist) {//find the closes neighbor
					neibdist=distance;
					neibtype=nextype;
					neibradius=nextradius;
				}
			} //end for
			if(neibtype==type_boundary_node)//adjust the interaction radius
				node.radius(neibradius);
			else
				node.radius(0.5F*(node.radius()+neibradius));
			if(nnodes==0) {
				System.out.println("setVelocity: nodes set empty");
				System.exit(1);
			}
			return energy;
		}
		void triangulate() {
				//ALGORITHM: Keep inserting new nodes as long as 
				// the average equilibrium interaction energy per node is less than minenergy	
			while (energy<minenergy) {
				insert();
				relax(mnodes.number()-bnodes.number());
if(mnodes.number()==100) break;///DDD
			}//END while(energy<minenergy);
		}
	}
	class ForceField {
		boolean done;
		float distav,distmax,time,dt,energy,minenergy,
			xmin[], xmax[];
		ForceField() {	
			init();
		//-	Potential.maxpotential=10.F;//heuristic
			Potential reference_potential=new Potential(1.0F,1.0F);
			minenergy=reference_potential.potential(0.9F);//heuristic
			if(fnodes==null||fnodes.number()<=0||fsegs==null||fsegs.number()==0) {
				System.out.println("Boundary is empty: meshing aborted");
				return;
			}
			done=false;
			time=0.0F;
			xmin=new float[dim];
			xmax=new float[dim];
			float[] xcenter=new float[dim];
			for(int i=0;i<dim;i++) {
				xcenter[i]=0.0F;
				xmin[i]=Float.MAX_VALUE;
				xmax[i]=-Float.MAX_VALUE;
			}
			fnodes.goFirst();
			do {
				Node node=fnodes.getItem();
				node.radius(0.0F); node.strength(1.0F);
				node.setZeroVelocity();
				for(int i=0;i<dim;i++) {
					float x=node.coordinate(i);
					xcenter[i]+=x;
					if(x<xmin[i]) xmin[i]=x;
					if(x>xmax[i]) xmax[i]=x;
				}
				fnodes.goNext();
			}	while(!fnodes.isFirst());
			domainsize=0.0F;
			int nnodes=fnodes.number();//initial number of nodes
			for(int i=0;i<dim;i++) {
				float x=xmax[i]-xmin[i];
				xcenter[i]/=(float)nnodes;
				domainsize=x*x;
			}
			domainsize=(float)Math.sqrt(domainsize);
			for(int i=0;i<dim;i++) {
				xmin[i]-=0.2*domainsize;
				xmax[i]+=0.2*domainsize;
			}
//-			Vector2D C=new Vector2D(xcenter[0],xcenter[1]);
			//Compute the average inter-node distance along the boundary:
			distav=0.0F;
			distmax=0.0F;
			fsegs.goFirst();
			do {
				Segment seg=fsegs.getItem();
				float seglen=seg.length(),
					r=0.5F*seglen;
				distav+=seglen;
				Node 
					a=seg.getNode(0),//first vertex of a segment
					b=seg.getNode(1);//second vertex
				float 
					ra=a.radius(),
					rb=b.radius();
				if(seglen>distmax)distmax=seglen;
				a.radius(1.0F*(ra+r));///0.5F*
				b.radius(1.0F*(rb+r));///0.5F*
				fsegs.goNext();
			}	while(!fsegs.isFirst());
			distav/=(float)fsegs.number();
			dt=0.01F*distav;//heuristics
			energy=0.0F;//average equilibrium interaction energy per node
		}
		void init() {
			Potential.maxpotential=10.F;//heuristic
			Potential.rmin=0.5F;
			Potential.rmax=1.5F;
		}
		void triangulate() {
				//ALGORITHM: Keep inserting new nodes as long as 
				// the average equilibrium interaction energy per node is less than minenergy	
			while (energy<minenergy) {
				insert();
				relax(mnodes.number()-bnodes.number());
if(mnodes.number()==100) break;///DDD
			}//END while(energy<minenergy);
		}
		void insert() {// Insert randomly:
System.out.println("insert nodes="+mnodes.number()+", energy="+energy);System.out.flush();///DDD
			Node node=new Node();
			node.type(type_internal_node);
			for(int i=0;i<dim;i++) {
				float x=xmin[i]+(xmax[i]-xmin[i])*(float)Math.random();
				node.coordinate(i,x);
			}
//+			node.radius(distav); 
node.radius(distmax);///DDD 
			node.strength(1.0F);//heuristic
//-			mnodes.goTag(); mnodes.goNext(); mnodes.insert(node);
			mnodes.add(node);
//-node.show("insert.node");///DDD
System.out.println("Inserting node "+node);///DDD
//+			relax(mnodes.number()-bnodes.number());
//-node.show("noderelaxed");///DDD
		}
		void relax(int maxiter) {
			if(mnodes==null||mnodes.number()<=0) {
				System.out.println("Relax: node set is empty - operation aborted");
				return;
			}
			//Reach equilibrium
			float eng=0.0F;
			int nnodes=0;
			float 
				v[]=new float[dim],
				x[]=new float[dim];
maxiter=1;///DDD
			for(int iter=0;iter<maxiter;iter++) {
//-System.out.println("\tinter="+iter);///DDD
				mnodes.goTag(); mnodes.goNext();//first internal node
				do {
					Node node=mnodes.getItem();
					node.setZeroVelocity();
					mnodes.goNext();
				}	while(!mnodes.isFirst());
				eng=0.0F; nnodes=0;
				mnodes.goTag(); mnodes.goNext();//first internal node
//-mnodes.goFirst();///DDD
				do {
//-System.out.println("\t\tnode "+nnodes+", energy="+eng);System.out.flush();///DDD
					Node node=mnodes.getItem();
//-if(node.type()==type_internal_node) {///DDD
					eng+=setVelocity();
//-					node.getVelocity(v);
//-System.out.println("v="+v[0]+"\t"+v[1]);///DDD
					for(int i=0;i<dim;i++) 
						x[i]=node.coordinate(i)+node.velocity(i)*dt;
					if(inside(x)) {node.setCoordinates(x);time+=dt;}
//-					handleBoundary();
					nnodes++;
//-}///DDD
					mnodes.goNext();
				}	while(!mnodes.isFirst());
			}//END for(iter=0..maxiter)
			if(nnodes==0) {
				System.out.println("Relax: No internal nodes found");
				return;
			}
			energy=eng/(float)nnodes;
			if(energy>minenergy) done=true;
System.out.println("Realx: energy="+energy+"; minenergy="+minenergy);///DDD
		}
		boolean inside(float[] x) {
			for(int i=0;i<dim;i++) 
				if(x[i]<xmin[i]||x[i]>xmax[i]) return false;
			return true;
		}
		float setVelocity() {
			float energy=0.0F,
				vector[]=new float[dim]; //velocity to be returned
			if(mnodes==null||mnodes.number()==0) {
				System.out.println("setVelocity: nodes set empty");
				System.exit(1);
			}
			mnodes.tag(1);
			Node node=mnodes.getItem();
			float	noderadius=node.radius(),
				neibradius=distav,
				neibdist=domainsize,
				distance=0.0F;
			int 
				neibtype=-1,//undefined
				nnodes=0;
			for (
				mnodes.goNext(); 
				!mnodes.isTag(1);
				mnodes.goNext()
			) {
				Node next=mnodes.getItem();
				//-neibradius=mnodes.getItem().radius(),
				int nextype=next.type();
				float	nextradius=next.radius();
				
	//next.show("nnodes="+nnodes+": next="+next);///DDD
				while((distance=next.distance(node,vector))<0.1*nextradius) {//heuristic
	//-System.out.println("*** ZERO distance="+distance);///DDD
					//Randomize vector:
					for(int i=0;i<dim;i++) 
						vector[i]=0.2F*noderadius*(float)(1.0-2.0*Math.random());//heuristic
					node.move(vector);
				}
				float
					potential=next.potential(1.2F*distance),///DDD//heuristic
					scale=potential/distance;
	//-System.out.println("node "+nnodes+": distance="+distance+", nextradius="+nextradius+", scale="+scale);///DDD
				for(int i=0;i<dim;i++) vector[i]*=scale;
				node.incVelocity(vector);
				energy+=potential;nnodes++;
				if(distance<neibdist) {//find the closes neighbor
					neibdist=distance;
					neibtype=nextype;
					neibradius=nextradius;
				}
			} //end for
			if(neibtype==type_boundary_node)//adjust the interaction radius
				node.radius(neibradius);
			else
				node.radius(0.5F*(node.radius()+neibradius));
			if(nnodes==0) {
				System.out.println("setVelocity: nodes set empty");
				System.exit(1);
			}
	//-System.out.println("setVeolcity: energy="+energy+", nnodes="+nnodes);///DDD
			return energy/(float)nnodes;
		}
		void handleBoundary() {
			//Find the closest boundary node to current mesh node (mnodes) 
			//and check if the current mesh node crossed over 
			//any of the segments attached to that boundary node (bnodes):
			bnodes.goFirst();
			do {
	
				bnodes.goNext();
			}	while(!bnodes.isFirst());
		}
	}
	void TriangulatePF() {
		//should follow Init()
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
points=new Collection<Node>();///DDD
		Vector2D //initialize all vectors before the loop for better efficiency:
			A=new Vector2D(),
			B=new Vector2D(),
			AB=new Vector2D(),
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
			while(true) {
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
					maxlength=1.5F*ab;
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
					if(!addnode) {
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
					} 
					}
//-					if(addnode) {
//-						CD.equal(CD1);CD.scale(cd);
//-						D.equal(C);D.add(CD);
//-						newnode=addNode(D,a,b,seg);
//-					}
				} else {//Create new triangle vertex: D
					CD.equal(CD1);CD.scale(cd);
					D.equal(C);D.add(CD);
					newnode=addNode(D,a,b,seg);
				}
				fsegs.delete();
				triangles.add(new Triangle(a,b,newnode));
				if(fsegs.number()==0) break;
//-				fsegs.goNext();///DDD
//+				fsegs.goLast();///DDD
			}
			smooth();
		}//end if(fsegs.number()>0)
	}
	public void smooth() {
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
