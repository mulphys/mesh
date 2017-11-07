/*
 */

/*
 * @(#)Demo.java	1.18 06/02/22
 */

/* A set of classes to parse, represent and display 3D wireframe models
   represented in Wavefront .obj format. */

import java.awt.*;
import java.util.*;
import java.awt.event.*;
import java.applet.Applet;
//import java.awt.Graphics;
//import java.awt.Color;
import java.awt.Event;
import java.io.*;
import java.net.URL;
import javax.swing.JButton;

class FileFormatException extends Exception {
	public FileFormatException(String s) {
		super(s);
	}
}
class Conf {
	public URL url;
	public String mdname;
	public int itransparent;
	public int model;
	public float zoom,xrot,yrot,zrot;
	public Conf() {
		url=null;
		mdname = null;
		itransparent=-1;
		zoom=1.0F;
		xrot=20.0F;
		yrot=20.0F;
		zrot=0.0F;
	}
}
/** The representation of a 3D model */
class Model3D {
	Conf conf;
	float vert[];//real coordinates
	int tvert[];//screen coordinates
	int vtype[];//vertex types: 0-reference points, 1-vertexes
	int nvert, maxvert;
	int con[];
	int ncon, maxcon;
	boolean transformed;
	Matrix3D mat;

	float xmin, xmax, ymin, ymax, zmin, zmax;

	Model3D () {
		mat = new Matrix3D ();
	}
	Model3D(Conf conf) throws IOException {
	/** Create a 3D model from a mesh */
		this();
		this.conf=conf;
		mat.xrot(conf.xrot);
		mat.yrot(conf.yrot);
	}
	/** Create a 3D model by parsing an input stream */
	Model3D (InputStream is, Conf conf) throws IOException, FileFormatException {
		this();
		this.conf=conf;
		mat.xrot(conf.xrot);
		mat.yrot(conf.yrot);
		StreamTokenizer st = new StreamTokenizer(
			new BufferedReader(new InputStreamReader(is, "UTF-8")));
		st.eolIsSignificant(true);
		st.commentChar('#');
//	st.parseNumbers();
		int ntype0=0;//not to add segments to reference points (type=0)
		scan:
		while (true) {
			switch (st.nextToken()) {
			  default:
			break scan;
			  case StreamTokenizer.TT_EOL:
			break;
			  case StreamTokenizer.TT_WORD:
			if (
				"p".equals(st.sval)
				||"v".equals(st.sval)
				||"n".equals(st.sval)
			) {
				int type=0;
				if("v".equals(st.sval))type=1;
				else if("n".equals(st.sval))type=2;
				double x = 0, y = 0, z = 0;
				if(type==0) ntype0++;
				if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
					x = readDouble(st);
					if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
						y = readDouble(st);
						if (st.nextToken() == StreamTokenizer.TT_NUMBER)
							z = readDouble(st);
					}
				}
				//addVert((float) x, (float) y, (float) z);
				addVert(type, (float) x, (float) y, (float) z);
				while (st.ttype != StreamTokenizer.TT_EOL &&
					st.ttype != StreamTokenizer.TT_EOF)
				st.nextToken();
			} else 
			if ("bs".equals(st.sval) || "s".equals(st.sval)) {
				int start = -1;
				int prev = -1;
				int n = -1;
				while (true)
				if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
					n = ntype0+(int) st.nval;//dont link type0 nodes into segments
					if (prev >= 0) add(prev - 1, n - 1);
					if (start < 0) start = n;
					prev = n;
				} else if (st.ttype == '/')
					st.nextToken();
					else
						break;
					if (start >= 0) add(start - 1, prev - 1);
					if (st.ttype != StreamTokenizer.TT_EOL)
					break scan;
				} else {
					while (st.nextToken() != StreamTokenizer.TT_EOL
						&& st.ttype != StreamTokenizer.TT_EOF);
				}
			}
		}
		is.close();
		if (st.ttype != StreamTokenizer.TT_EOF)
			throw new FileFormatException(st.toString());
	}
	public void reset() {
		vert=null;//real coordinates
		tvert=null;//screen coordinates
		vtype=null;//vertex types: 0-reference points, 1-vertexes
		nvert=maxvert=0;
		con=null;
		ncon=maxcon=0;
		transformed=false;
	}
	public void renderMesh(Mesh2D mesh) {
		//? this.mesh=mesh;
		Collection<Node> nodes=mesh.Nodes();
		if(nodes!=null&&nodes!=null&&nodes.number()>0) {
			Node[] pointer=new Node[nodes.number()];
			int inode=0;
			nodes.goFirst();
			do {
				Node node=nodes.getItem();
				int type=node.type();
				float 
					x=node.coordinate(0),
					y=node.coordinate(1),
//+					z=node.coordinate(2);
					z=0.0F;
				addVert(type, x, y, z);
				pointer[inode]=node;
				node.index=inode++;
				nodes.goNext();
			}	while(!nodes.isFirst());
		}
		Collection<Segment> segs=mesh.Segments();
		if(segs!=null&&segs!=null&&segs.number()>0) {
			segs.goFirst();
			int[] inode=new int[2];
			do {
				Segment seg=segs.getItem();
				for(int i=0;i<2;i++) {
					inode[i]=seg.getNode(i).index;
				}
				add(inode[0],inode[1]);
				segs.goNext();
			}	while(!segs.isFirst());
		}
	}
	/** Add a vertex to this model */
	int addVert(float x, float y, float z) {
		int i = nvert;
		if (i >= maxvert)
			if (vert == null) {
				maxvert = 100;
				vert = new float[maxvert * 3];
			} else {
				maxvert *= 2;
				float nv[] = new float[maxvert * 3];
				System.arraycopy(vert, 0, nv, 0, vert.length);
				vert = nv;
			}
		i *= 3;
		vert[i] = x;
		vert[i + 1] = y;
		vert[i + 2] = z;
		return nvert++;
	}
	/** Add a vertex and its type to this model */
	int addVert(int type, float x, float y, float z) {
		int i = nvert;
		if (i >= maxvert)
			if (vert == null) {
				maxvert = 100;
				vert = new float[maxvert * 3];
				vtype = new int[maxvert];
			} else {
				maxvert *= 2;
				float nv[] = new float[maxvert * 3];
				System.arraycopy(vert, 0, nv, 0, vert.length);
				vert = nv;
				int t[] = new int[maxvert];
				System.arraycopy(vtype, 0, t, 0, vtype.length);
				vtype = t;
			}
		vtype[nvert]=type;
		i *= 3;
		vert[i] = x;
		vert[i + 1] = y;
		vert[i + 2] = z;
		return nvert++;
	}
	void setVert(int i, float x, float y, float z) {
		if(i >= nvert) {
			System.out.println("VERTEX INDEX "+i+" EXCEEDS "+nvert);
			System.exit(1);
		}
		i *= 3;
		vert[i] = x;
		vert[i + 1] = y;
		vert[i + 2] = z;
	}
	/** Add a line from vertex p1 to vertex p2 */
	void add(int p1, int p2) {
		int i = ncon;
		if (p1 >= nvert || p2 >= nvert) return;
		if (i >= maxcon)
			if (con == null) {
				maxcon = 100;
				con = new int[maxcon];
			} else {
				maxcon *= 2;
				int nv[] = new int[maxcon];
				System.arraycopy(con, 0, nv, 0, con.length);
				con = nv;
			}
		if (p1 > p2) {
			int t = p1;
			p1 = p2;
			p2 = t;
		}
		con[i] = (p1 << 16) | p2;
		ncon = i + 1;
	}
	/** Transform all the points in this model */
	void transform() {
		if (transformed || nvert <= 0)
			return;
		if (tvert == null || tvert.length < nvert * 3)
			tvert = new int[nvert*3];
		mat.transform(vert, tvert, nvert);
		transformed = true;
	}

	/* Quick Sort implementation
	*/
	private void quickSort(int a[], int left, int right)
	{
		int leftIndex = left;
		int rightIndex = right;
		int partionElement;
		if ( right > left)
	{
		 /* Arbitrarily establishing partition element as the midpoint of
		  * the array.
		  */
		 partionElement = a[ ( left + right ) / 2 ];

		 // loop through the array until indices cross
		 while( leftIndex <= rightIndex )
		 {
			/* find the first element that is greater than or equal to
			 * the partionElement starting from the leftIndex.
			 */
			while( ( leftIndex < right ) && ( a[leftIndex] < partionElement ) )
				++leftIndex;

			/* find an element that is smaller than or equal to
			 * the partionElement starting from the rightIndex.
			 */
			while( ( rightIndex > left ) &&
					( a[rightIndex] > partionElement ) )
				--rightIndex;

			// if the indexes have not crossed, swap
			if( leftIndex <= rightIndex )
			{
				swap(a, leftIndex, rightIndex);
				++leftIndex;
				--rightIndex;
			}
		}

		 /* If the right index has not reached the left side of array
		  * must now sort the left partition.
		  */
		 if( left < rightIndex )
			quickSort( a, left, rightIndex );

		 /* If the left index has not reached the right side of array
		  * must now sort the right partition.
		  */
		if( leftIndex < right )
			quickSort( a, leftIndex, right );

		}
	}

	private void swap(int a[], int i, int j)
	{
		int T;
		T = a[i];
		a[i] = a[j];
		a[j] = T;
	}
	/** eliminate duplicate lines */
	void compress() {
	int limit = ncon;
	int c[] = con;
	quickSort(con, 0, ncon - 1);
	int d = 0;
	int pp1 = -1;
	for (int i = 0; i < limit; i++) {
		int p1 = c[i];
		if (pp1 != p1) {
			c[d] = p1;
			d++;
		}
		pp1 = p1;
	}
	ncon = d;
	}

	static int ncolors=4,iblack=3,ired=0,iblue=1,igreen=2;
	static Color color[][];

	/** Paint this model to a graphics context.  It uses the matrix associated
	with this model to map from model space to screen space.
	The next version of the browser should have double buffering,
	which will make this *much* nicer */
	void paint(Graphics g) {
	if (vert == null || nvert <= 0)
		return;
	transform();
	if (color == null) {
		color = new Color[ncolors][16];
		for (int i = 0; i < 16; i++) {
			int shade = (int) (170*(1-Math.pow(i/15.0, 2.3)));
			color[iblack][i] = new Color(shade, shade, shade);
			color[ired][i] = new Color(255, shade, shade);
			color[igreen][i] = new Color(shade, 255, shade);
			color[iblue][i] = new Color(shade, shade, 255);
		}
	}
//-	int lg = 0;
	float V[] = vert;
	int v[] = tvert;
	if (nvert <= 0) return;
	for (int i = 0; i < nvert; i++) {
		int 
			icolor=vtype[i],
			ip=3*i,
			shade = v[ip + 2],
			size=shade;
		if(icolor==conf.itransparent) continue;
		if(icolor>=ncolors)icolor=ncolors-1;
		if (shade < 0) shade = 0;
		if (shade > 15) shade = 15;
//-		if (shade != lg) {
//-			lg = shade;
			g.setColor(color[icolor][shade]);
//-		}
		size=shade;
		int radius=size/2;
//-		if(icolor==0)radius*=1.5;///DDD
		g.fillOval( v[ip]-radius, v[ip + 1]-radius, size,size);
	}
	//if (lim <= 0 || nvert <= 0) return;
	int lg = 0;
	int lim = ncon;
	int c[] = con;
	if (lim <= 0) return;
	for (int i = 0; i < lim; i++) {
		int T = c[i];
		int p1 = ((T >> 16) & 0xFFFF) * 3;
		int p2 = (T & 0xFFFF) * 3;
		int shade = v[p1 + 2] + v[p2 + 2];
		if (shade < 0) shade = 0;
		if (shade > 15) shade = 15;
		if (shade != lg) {
			lg = shade;
			g.setColor(color[iblack][shade]);
		}
		g.drawLine(v[p1], v[p1 + 1], v[p2], v[p2 + 1]);
	}
	}
	/** Find the bounding box of this model */
	void findBB() {
	if (nvert <= 0)
		return;
	float v[] = vert;
	float xmin = v[0], xmax = xmin;
	float ymin = v[1], ymax = ymin;
	float zmin = v[2], zmax = zmin;
	for (int i = nvert * 3; (i -= 3) > 0;) {
		float x = v[i];
		if (x < xmin)
		xmin = x;
		if (x > xmax)
		xmax = x;
		float y = v[i + 1];
		if (y < ymin)
		ymin = y;
		if (y > ymax)
		ymax = y;
		float z = v[i + 2];
		if (z < zmin)
		zmin = z;
		if (z > zmax)
		zmax = z;
	}
	this.xmax = xmax;
	this.xmin = xmin;
	this.ymax = ymax;
	this.ymin = ymin;
	this.zmax = zmax;
	this.zmin = zmin;
	}
	double readDouble(StreamTokenizer st) throws IOException, FileFormatException {
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
//				System.out.println("Num " + num * Math.pow(10, exp));
		return num * Math.pow(10, exp);
	}
}
class ControlPanel extends Panel { 
	ActionListener al;
	ItemListener il;
	GraphicsPanel graphicsPanel;
	Choice selectModel;
	JButton buttonRun;
	ControlPanel(EventListener listener, GraphicsPanel graphicsPanel) {
		al = (ActionListener)listener;
		il = (ItemListener)listener;
		this.graphicsPanel=graphicsPanel;
		Button b = new Button("Contour");
		b.addActionListener(al);
		add(b);
		buttonRun=new JButton("Mesh");
		buttonRun.addActionListener(al);
		add(buttonRun);
		b=new Button("Relax");
		b.addActionListener(al);
		add(b);
		selectModel = new Choice();
		add(selectModel);
		selectModel.addItem("Monte-Carlo");
		selectModel.addItem("Propagating Front");
		selectModel.addItem("Force-Field");
		selectModel.select("Monte-Carlo");
//+		b = new Button("previous");
//+		b.addActionListener(al);
//+		add(b);
//+
//+		add(new Label("go to:", Label.RIGHT));
//+
//+		Choice c = new Choice();
//+		c.addItemListener(il);
//+		add(c);
//+
//+		c.addItem("Arc");
//+		c.addItem("Oval");
//+		c.addItem("Polygon");
//+		c.addItem("Rect");
//+		c.addItem("RoundRect");
	}
}

class GraphicsPanel extends Panel
implements Runnable, MouseListener, MouseMotionListener
{ 
	Conf conf;
//-	ActionListener al;
//-	ItemListener il;
	Mesh2D mesh;
	Model3D md;
	boolean painted = true;
	float xfac;
	int prevx, prevy;
	float xtheta, ytheta;
	float scalefudge = 1;
	Matrix3D amat = new Matrix3D(), tmat = new Matrix3D();
	String mdname = null;
	String mdfile = null;
	String message = null;

	float deltaX = 0.0f, deltaY = 0.0f, scale = 1.0f;//new
	boolean left = false, center = false, right = false;//new

	GraphicsPanel(Conf conf) {
		this.conf=conf;
		mdname=conf.mdname;
		scalefudge = conf.zoom;//Float.valueOf(getParameter("zoom")).floatValue();
		loadBoundary();
		amat.xrot(1.0);///DDD: to prevent bug with when xrot=yrot=zrot=0
		amat.xrot(conf.xrot);
		amat.yrot(conf.yrot);
		amat.zrot(conf.zrot);
		setSize(500, 500);
		init();
	}
	public void loadBoundary() {
		mesh = new Mesh2D();
		InputStream is = null;
		try {
			is = new URL(conf.url, mdname).openStream();
			//-mesh.ReadBoundary(is); //replaces the line below.
			mesh.ReadControlPoints(is); mesh.GenerateBoundary(); //replaces mesh.ReadBoundary(is);
		}	catch(Exception e) {
			System.out.println(e.toString());
		}
		try {
			if (is != null)
			is.close();
		} catch(Exception e) {
		}
	}
	public void setBoundary() {
		mesh = new Mesh2D();
		InputStream is = null;
		try {
			is = new URL(conf.url, mdname).openStream();
			mesh.setControlPoints(is);
		}	catch(Exception e) {
			System.out.println(e.toString());
		}
		try {
			if (is != null)
			is.close();
		} catch(Exception e) {
		}
		mesh.GenerateBoundary();
	}
	public void init() {
		addMouseListener(this);
		addMouseMotionListener(this);
		try {
			Thread.currentThread().setPriority(Thread.MIN_PRIORITY);
			mesh.Init();
			Model3D m = new Model3D (conf);
			m.renderMesh(mesh);
			md = m;
			m.findBB();
			m.compress();
			float xw = m.xmax - m.xmin;
			float yw = m.ymax - m.ymin;
			float zw = m.zmax - m.zmin;
			if (yw > xw)
			xw = yw;
			if (zw > xw)
			xw = zw;
			float f1 = getSize().width / xw;
			float f2 = getSize().height / xw;
			xfac = 0.7f * (f1 < f2 ? f1 : f2) * scalefudge;
		} catch(Exception e) {
			md = null;
			message = e.toString();
		}
		repaint();
	}
	public Mesh2D getMesh() { return mesh; }
	public void start() {
	if (md == null && message == null)
		new Thread(this).start();
	}
	public void stop() {
	}

	public void destroy() {
		removeMouseListener(this);
		removeMouseMotionListener(this);
	}

	public void run() {
		repaint();
	}

	public  void mouseClicked(MouseEvent e) {
	}
	public  void mousePressed(MouseEvent e) {
		if (e.getButton() == e.BUTTON1) left = true;
		else if (e.getButton() == e.BUTTON2) center = true;
		else if (e.getButton() == e.BUTTON3) right = true;
		prevx = e.getX();
		prevy = e.getY();
		e.consume();
	}

	public  void mouseReleased(MouseEvent e) {
		if (e.getButton() == e.BUTTON1) left = false;
		else if (e.getButton() == e.BUTTON2) center = false;
		else if (e.getButton() == e.BUTTON3) right = false;
	}

	public  void mouseEntered(MouseEvent e) {
	}

	public  void mouseExited(MouseEvent e) {
	}

	public  void mouseDragged(MouseEvent e) {
		int x = e.getX();
		int y = e.getY();
//-		tmat.unit();
//-		float xtheta = (prevy - y) * 360.0f / getSize().width;
//-		float ytheta = (x - prevx) * 360.0f / getSize().height;
//-		tmat.xrot(xtheta);
//-		tmat.yrot(ytheta);
//-		amat.mult(tmat);
		if (center || (left && right)) {
			deltaX += -(prevx - x);
			deltaY += -(prevy - y);
		}
		else if (left) {
			tmat.unit();
			float
				xtheta = (prevy - y) * (360.0f / getSize().width),
				ytheta = (x - prevx) * (360.0f / getSize().height);
			tmat.xrot(xtheta);
			tmat.yrot(ytheta);
			amat.mult(tmat);
		}
		else if (right) {
			scale += (prevy - y)/180.0f;
			if (scale < 0.1f)
				scale = 0.1f;
		}

		if (painted) {
			painted = false;
			repaint();
		}
		prevx = x;
		prevy = y;
		e.consume();
	}

	public  void mouseMoved(MouseEvent e) {
	}

	public void paint(Graphics g) {
	if (md != null) {
		md.mat.unit();
		md.mat.translate(-(md.xmin + md.xmax) / 2,
				 -(md.ymin + md.ymax) / 2,
				 -(md.zmin + md.zmax) / 2);
		md.mat.mult(amat);
		md.mat.scale(xfac, -xfac, 16 * xfac / getSize().width);
		md.mat.scale(scale,scale,scale);//new
		md.mat.translate(getSize().width / 2, getSize().height / 2, 8);
		md.mat.translate(deltaX, deltaY, 0.0f);//new
		md.transformed = false;
		md.paint(g);
		setPainted();
	} else if (message != null) {
		g.drawString("Error in model:", 3, 20);
		g.drawString(message, 10, 40);
	}
	}

	private synchronized void setPainted() {
	painted = true;
	notifyAll();
	}
//	private synchronized void waitPainted() {
//	while (!painted)
//		wait();
//	painted = false;
//	}

	public String getAppletInfo() {
		return "Title: Demo \nAuthor: Andrei.V.Smirnov@gmail.com \nAn applet to put a 3D model into a page.";
	}

	public String[][] getParameterInfo() {
		String[][] info = {
			{"model", "path string", "The path to the model to be displayed."},
			{"scale", "float", "The scale of the model.  Default is 1."}
		};
		return info;
	}
}
/** An applet to put a 3D model into a page */
public class Demo extends Applet
implements	ActionListener, ItemListener, Runnable {
	boolean done,running,suspended;
	Thread thread;
	ControlPanel controlPanel;
	public GraphicsPanel graphicsPanel;
	public Conf conf;
	public Mesh2D mesh;
	public void init() {
		thread=null;
		conf=new Conf();
		try {
			conf.url=getDocumentBase();
			conf.mdname = getParameter("model");
			conf.zoom = Float.valueOf(getParameter("zoom")).floatValue();
			String str;
			if ((str=getParameter("xrot")) != null) conf.xrot=(Float.valueOf(str)).floatValue();
			if ((str=getParameter("yrot")) != null) conf.yrot=(Float.valueOf(str)).floatValue();
			if ((str=getParameter("zrot")) != null) conf.zrot=(Float.valueOf(str)).floatValue();
		} catch (Exception e) {
			System.out.println("Failed to read parameters: "+e.toString());
		};
		setLayout(new BorderLayout());
		//setSize(600, 600);
		add("Center", graphicsPanel = new GraphicsPanel(conf));
		add("South", controlPanel = new ControlPanel(this,graphicsPanel));
		mesh=graphicsPanel.getMesh();
		conf.model=controlPanel.selectModel.getSelectedIndex();
		mesh.init(conf.model);
		done=running=suspended=false;
	}
	public void destroy() {
		remove(controlPanel);
		remove(graphicsPanel);
	}
	public void start() {
		done=false;
	}
	public void run() {
		if(suspended) return;
		//((CardLayout)graphicsPanel.cards.getLayout()).next(graphicsPanel.cards);
		//conf.itransparent=conf.itransparent==0?-1:0;
			int nbnodes=mesh.bnodes.number();
			while(!done&&!suspended) {
				mesh.insert(conf.model);
				int nrelax=mesh.mnodes.number()-nbnodes;
				if(conf.model==Mesh2D.model_monte_carlo) {
					mesh.montecarlo.relax(nrelax);
					done=mesh.montecarlo.done;
				}
				else {
					mesh.forcefield.relax(nrelax);
					done=mesh.forcefield.done;
				}
//-				mesh.mnodes.goTag(); mesh.mnodes.goNext();//first internal node
				mesh.mnodes.goTag(); mesh.mnodes.goNext();
				int inode=nbnodes;
				while(!mesh.mnodes.isLast()) {
					Node node=mesh.mnodes.getItem();
					int type=node.type();
					float
						x=node.coordinate(0),
						y=node.coordinate(1),
//+					z=node.coordinate(2);
						z=0.0F;
					graphicsPanel.md.setVert(inode++,x, y, z);
					mesh.mnodes.goNext();
				}
				Node node=mesh.mnodes.getItem();
				int type=node.type();
				float 
					x=node.coordinate(0),
					y=node.coordinate(1),
//+				z=node.coordinate(2);
					z=0.0F;
				graphicsPanel.md.addVert(type, x, y, z);
				//-graphicsPanel.md.transform();
				graphicsPanel.repaint();
			}
			if(done)controlPanel.buttonRun.setText("Done");
	}
	public void actionPerformed(ActionEvent e) {
		String arg = e.getActionCommand();

		if ("Contour".equals(arg)) {
			graphicsPanel.setBoundary();
			graphicsPanel.init();
			mesh=graphicsPanel.getMesh();
			graphicsPanel.md.reset();
			graphicsPanel.md.renderMesh(mesh);
		}
		else {
			int imodel=controlPanel.selectModel.getSelectedIndex();
			if(!running)conf.model=imodel;
			if ("Mesh".equals(arg)) {
		switch(imodel) {
		case Mesh2D.model_force_field: //Force Field:
		case Mesh2D.model_monte_carlo: //Force Field:
			if(conf.model!=imodel) {
				mesh.init(imodel);
				conf.model=imodel;
			}
			if(done) 
				controlPanel.buttonRun.setText("Done");
			else
				if(!running) {
//-					new Thread(this).start();
					controlPanel.buttonRun.setText("Pause");
				}	else { 
					if(suspended) {
					//-	thread.resume(); 
						controlPanel.buttonRun.setText("Pause");
						suspended=false;
					}
					else {
					//-	thread.suspend();
						controlPanel.buttonRun.setText("Resume");
						suspended=true;
					}
				}
				thread=new Thread(this);
				thread.start();
				running=true;
				suspended=false;
			break;
		case Mesh2D.model_propagating_front: //Propagating Front
			conf.model=imodel;
//-		if(conf.model!=imodel) {
//?			init();
//?			mesh=graphicsPanel.getMesh();
//-			graphicsPanel.md.reset();
//-			graphicsPanel.md.renderMesh(mesh);
//-			graphicsPanel.repaint();
//-		}
			mesh.Triangulate(imodel);
			graphicsPanel.md.reset();
			graphicsPanel.md.renderMesh(mesh);
//-			graphicsPanel.md.transform();
			graphicsPanel.repaint();
			break;
		}
		}
		else if ("Relax".equals(arg)) {
			if(mesh==null) return;
			int nrelax=mesh.mnodes.number()-mesh.bnodes.number();
			if(conf.model!=imodel) {
				mesh.init(imodel);
				conf.model=imodel;
			}
			switch(imodel) {
			case Mesh2D.model_propagating_front:
				System.out.println("Propagating Front Model Not Implemented");
				break;
			case Mesh2D.model_monte_carlo:
				mesh.montecarlo.relax(nrelax);
				done=mesh.montecarlo.done;
				setMeshVerts();
				break;
			case Mesh2D.model_force_field:
				mesh.forcefield.relax(nrelax);
				done=mesh.forcefield.done;
				setMeshVerts();
				break;
			} //end switch
			if(done) {
				controlPanel.buttonRun.setText("Done");
				suspended=false;
				running=false;
			}
			else {
				controlPanel.buttonRun.setText("Mesh");
				running=false;
				suspended=false;
			}
		}	else if("Pause".equals(arg)) {
			controlPanel.buttonRun.setText("Resume");
			suspended=true;
		}	else if("Resume".equals(arg)) {
			controlPanel.buttonRun.setText("Pause");
					thread=new Thread(this);
					thread.start();
			running=true;
			suspended=false;
		}
		}
		graphicsPanel.repaint();
	}
	void setMeshVerts() {
		int inode=mesh.bnodes.number();
		if(inode==0) return;
		mesh.mnodes.goTag(); mesh.mnodes.goNext();//first internal node
		do {
			Node node=mesh.mnodes.getItem();
			int type=node.type();
				float 
				x=node.coordinate(0),
				y=node.coordinate(1),
//+			z=node.coordinate(2);
				z=0.0F;
			graphicsPanel.md.setVert(inode++,x, y, z);
			mesh.mnodes.goNext();
		}	while(!mesh.mnodes.isFirst()); 
//-		graphicsPanel.md.transform();
	}
	public void itemStateChanged(ItemEvent e) {
System.out.println("itemStateChanged(ItemEvent e): NOT IMPLEMENTED");///DDD
		//((CardLayout)graphicsPanel.cards.getLayout()).show(graphicsPanel.cards,(String)e.getItem());
	}
	public static void main(String args[]) {
		AppletFrame.startApplet("Demo", "Demo", args);
	}
	public String getAppletInfo() {
		return "An interactive demonstration of some graphics.";
	}
} 

