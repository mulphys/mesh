class Point {
	static final int dim=2;
	public float[] coordinates;
	public Point() {
		coordinates=new float[dim];
	}
	public Point(float x, float y) {
		this();
		coordinates[0]=x;
		coordinates[1]=y;
	}
	public float coordinate(int i) { 
		if(i>dim||i<0) {
			System.out.println("Point index "+i+" exceeds dimension "+dim);
			System.exit(1);
		}
		return coordinates[i]; 
	}
	public void coordinate(int i, float value) { 
		if(i>dim||i<0) {
			System.out.println("Point index "+i+" exceeds dimension "+dim);
			System.exit(1);
		}
		coordinates[i]=value; 
	}
	public void coordinates(float x0, float x1) { 
		coordinates[0]=x0; 
		coordinates[1]=x1; 
	}
}
class Contour {
	int 
		mp,np; //allocated and actual numbers of reference points
//-		mv,nv; //allocated and actual number of vertexes
	Point[] points;//reference points
//-		verts;//contour vertexes
	Contour() {
		mp=3;//Cubic
		np=0;
//-nv=0;
		points=new Point[mp];
//-vertexes=null;
	}
	Contour(int n) {
		mp=n;
		np=0;
		points=new Point[n];
	}
	public int maxNumberPoints() { return mp; }
	public int numberPoints() { return np; }
//-	public void allocateVertexes(int n) {
//-		mv=n;
//-		verts=new Point[n];
//-	}
	public void addPoint(float x, float y) {
		if(np==mp){
			System.out.println("CONTROL POINT NO. "+np+" EXCEEDS MAX "+mp);
			return;
		}
		points[np++]=new Point(x, y);
	}
	public void setPoint(int i, float x, float y) {
		if(i>=np){
			System.out.println("CONTROL POINT INDEX "+i+" EXCEEDS MAX "+np);
			return;
		}
		points[i].coordinates(x,y);
	}
	public Point point(int i) {
		if(i>np) {
			System.out.println("CAN'T ACCESS POINT "+i+" BEYOND MAX "+np);
			System.exit(1);
		}
		return points[i];
	}
	public float[] vertex(
		float t //time-parameter: position on the contour line
	) {
		final int dim=2;
		float[]   
			a=new float[dim], 
			b=new float[dim], 
			c=new float[dim],
			x=new float[dim];
		float   tSquared, tCubed;
		/* calculate the polynomial coefficients */
		for(int i=0;i<dim;i++) {
			c[i] = 3.0F * (points[1].coordinate(i) - points[0].coordinate(i));
			b[i] = 3.0F * (points[2].coordinate(i) - points[1].coordinate(i)) - c[i];
			//+ a[i] = points[3].coordinate(i) - points[0].coordinate(i) - c[i] - b[i];
			a[i] = - (c[i] + b[i]);// makes it a contour
		}
		    
		/* calculate the curve point at parameter value t */
		    
		tSquared = t * t;
		tCubed = tSquared * t;

		for(int i=0;i<dim;i++) 
			x[i] = a[i]*tCubed + b[i]*tSquared + c[i]*t + points[0].coordinate(i);	    
		return x;
	}
}

