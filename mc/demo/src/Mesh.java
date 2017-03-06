/*
	Monte-Carlo Mesh Generator
	Author: 
		Andrei Smirnov (andrei.v.smirnov@gmail.com)
		Hanzhout Zhang (hanzhou.zhang@gmail.com
*/

/*
	To get the Render mode switch uncomment all lines with comboBox
*/

import java.applet.Applet;
import java.awt.Image;
import java.awt.Event;
import java.awt.Graphics;
import java.awt.Dimension;
import java.io.*;
import java.net.URL;
import java.util.Hashtable;
import java.awt.image.IndexColorModel;
import java.awt.image.ColorModel;
import java.awt.image.MemoryImageSource;
import java.awt.event.*;
import javax.swing.*;
import java.awt.*;
import java.util.Vector;

/** Wireframe */

class WireFrame 
{
	static int DIM=3;
	int dim,Dim[];
	static Color col[][];
	//Particles Injection data:
	static int nnodes,mnodes,//number of injection points
		inj0,inj1;//injection point and the first particle 
///SIM	static int trapped=0; //number of particels trapped
	static double[] 
		X0,dX,DX,//Domain origin and dimensions
		Xinj=null, //injection coordinates
		Vair=null; //Air velocity
	static float vairscale=3.0F;
///		Vel; //particle velocities
	static double	vdir;// velocity direction
	static int[]	Cinj; //injection points color
	static byte[] grid; //domain geometry
	//Points:
	int nvert, mvert, maxvert;
	float vert[];
	int tvert[];
///	int ZsortMap[];
	int ncon, maxcon;
	int con[];
	byte com[];
	int	nshades;
	boolean transformed;
	Matrix3D mat;
	float xmin, xmax, ymin, ymax, zmin, zmax;

///	private final static int bgGrey = 192;
///	private static int maxr;
///	private int Rl;
///	private int Gl;
///	private int Bl;
///	static 
///	{	data = new byte[R * 2 * R * 2];
///		int mr = 0;
///		for (int Y = 2 * R; --Y >= 0;) 
///		{	int x0 = (int) (Math.sqrt(R * R - (Y - R) * (Y - R)) + 0.5);
///			int p = Y * (R * 2) + R - x0;
///			for (int X = -x0; X < x0; X++) 
///			{	int x = X + hx;
///				int y = Y - R + hy;
///				int r = (int) (Math.sqrt(x * x + y * y) + 0.5);
///				if (r > mr)
///				mr = r;
///				data[p++] = r <= 0 ? 1 : (byte) r;
///			}
///		}
///		maxr = mr;
///	}
///	private final int blend(int fg, int bg, float fgfactor) 
///	{
///		return (int) (bg + (fg - bg) * fgfactor);
///	}

	WireFrame () 
	{
		mat = new Matrix3D ();
		mat.xrot(20);
		mat.yrot(30);
		ncon=0;
		nshades=16;
	}
	WireFrame (Space space, Domain domain, Vector components) 
	{	this();
		dim=space.getDim();
		Dim=domain.getDim();
		int
			ncomp=components.size();
		col = new Color[nshades][ncomp];
		dX = space.getdX();
		DX = domain.getDX();
		X0 = domain.getX0();
		grid=domain.getgrid();
		float[]
			dY = new float[dim],
			Y0 = new float[dim];
		int nv=0;//vertex counter
		vert=null;
		tvert=null;
///		ZsortMap=null;
		nvert=nv;
		//Set components
		for (int icomp=1;icomp<ncomp;icomp++)
		{	Components	comp=(Components)components.elementAt(icomp);
			for (int ishade = 0; ishade < nshades; ishade++) 
			{	double shade=1.0*(1-Math.pow(ishade/((float)nshades-1.0), 2.3));
				int
					grey=170,
					r=(int)((1-shade)*comp.color.getRed()+shade*grey),
					g=(int)((1-shade)*comp.color.getGreen()+shade*grey),
					b=(int)((1-shade)*comp.color.getBlue()+shade*grey);
				col[ishade][icomp] = new Color(r,g,b);
			}
			float[]
				x1 = new float[dim],
				y1 = new float[dim];
			if(comp.boundaryContours==null)continue;
			//ADD BOUNDARY CONTOURS
			for(int idir=0;idir<dim;idir++)
			{	if(comp.boundaryContours[idir]==null)continue;
				BoundaryContour	contour=comp.boundaryContours[idir];
				while (contour!=null)
				{
					x1[0]=(float)contour.getVertexIndX(0);
					x1[1]=(float)contour.getVertexIndY(0);
					x1[2]=(float)contour.getPlaneInd();
					for (int i=0;i<dim;i++)
					{	int	j=(i+dim-idir+2)%dim;
						y1[i]=(float)X0[j]+(float)dX[j]*x1[j];
					}
					addVert(y1[0], y1[1], y1[2]);nv++;
					for(int ic=1;ic<contour.Length();ic++)
					{
						x1[0]=(float)contour.getVertexIndX(ic);
						x1[1]=(float)contour.getVertexIndY(ic);
						for (int i=0;i<dim;i++)
						{	int	j=(i+dim-idir+2)%dim;
							y1[i]=(float)X0[j]+(float)dX[j]*x1[j];
						}
						addVert(y1[0], y1[1], y1[2]);
						add(nv-1,nv,icomp);nv++;
					}
					contour=contour.next;
				}
			}
		}
		//Set particles color (strored as 0-component color)
		for (int ishade = 0; ishade < nshades; ishade++) 
		{	double shade=1.0*(1-Math.pow(ishade/((float)nshades-1.0), 2.3));
			int
				grey=170,
				r=(int)((1-shade)*Cinj[0]+shade*grey),
				g=(int)((1-shade)*Cinj[1]+shade*grey),
				b=(int)((1-shade)*Cinj[2]+shade*grey);
			col[ishade][0] = new Color(r,g,b);
		}
//SIM		{	//ADD PARTICLES
//SIM			float[] xinj = new float[dim];
//SIM			for(int i=0;i<dim;i++)
//SIM				xinj[i]=(float)Math.floor(DX[i]*Xinj[i]);
//SIM			///for(int iinj=0;iinj<nnodes;iinj++)
//SIM			for(int iinj=0;iinj<nnodes+1;iinj++)
//SIM			{	
//SIM				addVert((float)xinj[0], (float)xinj[1], (float)xinj[2]);
//SIM				nv++;
//SIM			}
//SIM			inj0=nvert-nnodes-1; 
//SIM			inj1=inj0+1; 
//SIM			int  jnj0=DIM*inj0,jnj1=DIM*inj1;
//SIM			for (int i=0;i<DIM;i++)
//SIM			{	vert[jnj1+i]=vert[jnj0+i]+vairscale*(float)Vair[i];
//SIM			}
//SIM			add(inj0,inj1,0);
//SIM		}
	}
	/** Create a 3D model by parsing an input stream */
///	WireFrame (InputStream is) throws IOException, FileFormatException 
///	{	this();
///		StreamTokenizer st = new StreamTokenizer(new BufferedReader(new InputStreamReader(is)));
///		st.eolIsSignificant(true);
///		st.commentChar('#');
///	scan:
///	while (true) {
///		switch (st.nextToken()) {
///			default:
///		break scan;
///			case StreamTokenizer.TT_EOL:
///		break;
///			case StreamTokenizer.TT_WORD:
///		if ("v".equals(st.sval)) {
///			double x = 0, y = 0, z = 0;
///			if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
///			x = st.nval;
///			if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
///				y = st.nval;
///				if (st.nextToken() == StreamTokenizer.TT_NUMBER)
///				z = st.nval;
///			}
///			}
///			addVert((float) x, (float) y, (float) z);
///			while (st.ttype != StreamTokenizer.TT_EOL &&
///				st.ttype != StreamTokenizer.TT_EOF)
///			st.nextToken();
///		} else if ("f".equals(st.sval) || "fo".equals(st.sval) || "l".equals(st.sval)) {
///			int start = -1;
///			int prev = -1;
///			int n = -1;
///			while (true)
///			if (st.nextToken() == StreamTokenizer.TT_NUMBER) {
///				n = (int) st.nval;
///				if (prev >= 0)
///				add(prev - 1, n - 1);
///				if (start < 0)
///				start = n;
///				prev = n;
///			} else if (st.ttype == '/')
///				st.nextToken();
///			else
///				break;
///			if (start >= 0)
///			add(start - 1, prev - 1);
///			if (st.ttype != StreamTokenizer.TT_EOL)
///			break scan;
///		} else {
///			while (st.nextToken() != StreamTokenizer.TT_EOL
///				&& st.ttype != StreamTokenizer.TT_EOF);
///		}
///		}
///	}
///	is.close();
///	if (st.ttype != StreamTokenizer.TT_EOF)
///		throw new FileFormatException(st.toString());
///	}
	/** Delete a vertex from this model */
	int delVert(int ivert) 
	{
		if (ivert >= nvert) return nvert;
		nvert--;
		for (int i=ivert;i<nvert;i++)
		{	int j = 3*i, j1=j+3;
			for(int k=0;k<3;k++)
				vert[j+k] = vert[j1+k];
		}
		return nvert;
	}
	/** Add a vertex to this model */
	int addVert(float x, float y, float z) 
	{
		int i = nvert;
		if (i >= maxvert)
			if (vert == null) 
			{
				maxvert = 100;
				vert = new float[maxvert * 3];
			} else 
			{	maxvert *= 2;
				float nv[] = new float[maxvert * 3];
				System.arraycopy(vert, 0, nv, 0, vert.length);
				vert = nv;
			}
		i *= 3;
		vert[i] = x;
		vert[i + 1] = y;
		vert[i + 2] = z;
		if(mvert<=nvert)mvert=nvert+1;
		return nvert++;
	}
	/** Add a line from vertex p1 to vertex p2 */
	void add(int p1, int p2, int icomp) 
	{
		int i = ncon;
		if (p1 >= nvert || p2 >= nvert)
			return;
		if (i >= maxcon)
		if (con == null) 
		{
			maxcon = 100;
			con = new int[maxcon];
			com = new byte[maxcon];
		}	else 
		{	maxcon *= 2;
			int nv[] = new int[maxcon];
			System.arraycopy(con, 0, nv, 0, con.length);
			con = nv;
			byte nc[] = new byte[maxcon];
			System.arraycopy(com, 0, nc, 0, com.length);
			com = nc;
		}
///		if (p1 > p2) 
///		{	int t = p1;
///			p1 = p2;
///			p2 = t;
///		}
		con[i] = (p1 << 16) | p2;
		com[i] = (byte)icomp;
		ncon = i + 1;
	}
	/** Transform all the points in this model */
	void transform() 
	{
		if (transformed || nvert <= 0)
			return;
		if (tvert == null || tvert.length < nvert * 3)
			tvert = new int[nvert*3];
		mat.transform(vert, tvert, nvert);
		transformed = true;
	}
//QS	/* Quick Sort implementation
//QS	*/
//QS	private void quickSort(int a[], int left, int right)
//QS	{
//QS		int leftIndex = left;
//QS		int rightIndex = right;
//QS		int partionElement;
//QS		if ( right > left)
//QS		{	/* Arbitrarily establishing partition element as the midpoint of
//QS				* the array.
//QS				*/
//QS			partionElement = a[ ( left + right ) / 2 ];
//QS
//QS			// loop through the array until indices cross
//QS			while( leftIndex <= rightIndex )
//QS			{
//QS				/* find the first element that is greater than or equal to
//QS				 * the partionElement starting from the leftIndex.
//QS				 */
//QS				while( ( leftIndex < right ) && ( a[leftIndex] < partionElement ) )
//QS					++leftIndex;
//QS				/* find an element that is smaller than or equal to
//QS				 * the partionElement starting from the rightIndex.
//QS				 */
//QS				while( ( rightIndex > left ) &&
//QS					( a[rightIndex] > partionElement ) )
//QS					--rightIndex;
//QS				// if the indexes have not crossed, swap
//QS				if( leftIndex <= rightIndex )
//QS				{
//QS					swap(a, leftIndex, rightIndex);
//QS					++leftIndex;
//QS					--rightIndex;
//QS				}
//QS			}
//QS			/* If the right index has not reached the left side of array
//QS				* must now sort the left partition.
//QS				*/
//QS			if( left < rightIndex )
//QS				quickSort( a, left, rightIndex );
//QS			/* If the left index has not reached the right side of array
//QS				* must now sort the right partition.
//QS				*/
//QS			if( leftIndex < right )
//QS				quickSort( a, leftIndex, right );
//QS		}
//QS	}
//QS	private void swap(int a[], int i, int j)
//QS	{
//QS		int T;
//QS		T = a[i];
//QS		a[i] = a[j];
//QS		a[j] = T;
//QS	}
//QS	/** eliminate duplicate lines */
//QS	void compress() 
//QS	{
//QS		int limit = ncon;
//QS		int c[] = con;
//QS		quickSort(con, 0, ncon - 1);
//QS		int d = 0;
//QS		int pp1 = -1;
//QS		for (int i = 0; i < limit; i++) 
//QS		{	int p1 = c[i];
//QS			if (pp1 != p1) 
//QS			{	c[d] = p1;
//QS				d++;
//QS			}
//QS			pp1 = p1;
//QS		}
//QS		ncon = d;
//QS	}
	/** Paint this model to a graphics context.	It uses the matrix associated
	with this model to map from model space to screen space.
	The next version of the browser should have double buffering,
	which will make this *much* nicer */
	void paint(Graphics g) 
	{
		if (vert == null || nvert <= 0)
			return;
		transform();
		int v[] = tvert;
///		int zs[] = ZsortMap;
///		if (zs == null) 
///		{
///			ZsortMap = zs = new int[nvert];
///			for (int i = nvert; --i >= 0;)
///			zs[i] = i * 3;
///		}
///		/*
///		* I use a bubble sort since from one iteration to the next, the sort
///		* order is pretty stable, so I just use what I had last time as a
///		* "guess" of the sorted order.  With luck, this reduces O(N log N)
///		* to O(N)
///		*/
///		for (int i = nvert - 1; --i >= 0;) 
///		{	boolean flipped = false;
///			for (int j = 0; j <= i; j++) 
///			{	int a = zs[j];
///				int b = zs[j + 1];
///				if (v[a + 2] > v[b + 2]) 
///				{
///					zs[j + 1] = a;
///					zs[j] = b;
///					flipped = true;
///				}
///			}
///			if (!flipped)
///			break;
///		}
		int lg = 0;
		int lim = ncon;
		int c[] = con;
		if (lim <= 0 || nvert <= 0) return;
		for (int i = 0; i < lim; i++) 
		{	int T = c[i];
			int 
				p1=((T >> 16) & 0xFFFF) * 3,
				p2=(T & 0xFFFF) * 3,
				grey = v[p1 + 2] + v[p2 + 2];
			if (grey < 0)grey = 0;
			if (grey >= nshades)	grey = nshades-1;
			if (i*nshades+grey != lg) 
			{	lg = i*nshades+grey;
				g.setColor(col[grey][(int)com[i]]);
			}
			g.drawLine
			(	v[p1], v[p1 + 1],
				v[p2], v[p2 + 1]
			);
		}
		//Display particles:
		{	//Injection point
			////int i=nvert-nnodes,
			int i=inj0,
				p=3*i,
				grey = v[p + 2];
			if (grey < 0)grey = 0;
			if (grey >= nshades)	grey = nshades-1;
			if (i*nshades+grey != lg) 
			{	lg = i*nshades+grey;
				g.setColor(col[grey][0]);
			}
			g.fillOval
			(	v[p]-4, v[p + 1]-4,
				8, 8 // SIZE OF THE POINT SOURCE
			);
		}
		//Other points
		for (int i = nvert-nnodes+2; i < nvert; i++) 
	///	for (int i = inj1; i < nvert; i++) 
		{	int 
				p=3*i,
				grey = v[p + 2];
			if (grey < 0)grey = 0;
			if (grey >= nshades)	grey = nshades-1;
			if (i*nshades+grey != lg) 
			{	lg = i*nshades+grey;
				g.setColor(col[grey][0]);
			}
			g.fillOval
			(	v[p], v[p + 1],
				1, 1 // SIZE OF THE POINT SOURCE
			);
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
	for (int i = nvert * 3; (i -= 3) > 0;) 
	{	float x = v[i];
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
	//PARTICLE ROUTINES:
///SIM	static public void	setVair(double degrees)
///SIM	{	double vel=1.0,
///SIM			radians=Math.PI*degrees/180.0;//Velocity magnitude
///SIM		vdir=degrees;
///SIM		Vair[0]=vel*Math.cos(radians);
///SIM		Vair[1]=vel*Math.sin(radians);
///SIM		Vair[2]=0.0;
///SIM	}
///SIM	static public double[]	getVair()
///SIM	{
///SIM		return Vair;
///SIM	}
///SIM	static public void	setXinj(double xinj, double yinj, double zinj)
///SIM	{
///SIM		Xinj[0]=xinj;
///SIM		Xinj[1]=yinj;
///SIM		Xinj[2]=zinj;
///SIM		for(int i=0;i<DIM;i++)
///SIM		{
///SIM			if(Xinj[i]<0.1)
///SIM			{	if(Xinj[i]<0.0)System.out.println("WARNING: Injectction coordinate "+Xinj[i]+" should be positive. Set to 0.0.\n");
///SIM				Xinj[i]=0.1F;
///SIM			}
///SIM			if(Xinj[i]>0.9)
///SIM			{	if(Xinj[i]>1.0)System.out.println("WARNING: Injectction coordinate "+Xinj[i]+" should be no greater than 1. Set to 1.\n");
///SIM				Xinj[i]=0.9F;
///SIM			}
///SIM		}
///SIM	}
///SIM	static public void initinj
///SIM	(	int n, 
///SIM		int red, 
///SIM		int green, 
///SIM		int blue, 
///SIM		float xinj, 
///SIM		float yinj, 
///SIM		float zinj,
///SIM		float vdir
///SIM	)
///SIM	{
///SIM		if(n<0) n=0;
///SIM		Vair = new double[DIM];
///SIM		Vair[0]=1.0;
///SIM		Vair[1]=0.0;
///SIM		Vair[2]=0.0;
///SIM		mnodes=n;
///SIM		nnodes=mnodes;
///SIM		Cinj = new int[DIM];
///SIM		Cinj[0]=red;
///SIM		Cinj[1]=green;
///SIM		Cinj[2]=blue;
///SIM		for(int i=0;i<3;i++)
///SIM		{
///SIM			if(Cinj[i]<0)
///SIM			{	System.out.println
///SIM				("WARNING: Injectction color "+Cinj[i]+" should be positive. Set to 0.\n");
///SIM				Cinj[i]=0;
///SIM			}
///SIM			if(Cinj[i]>255)
///SIM			{	System.out.println
///SIM				("WARNING: Injectction color "+Cinj[i]+" should be no greater than 255. Set to 255.\n");
///SIM				Cinj[i]=255;
///SIM			}
///SIM		}
///SIM		Xinj = new double[DIM];
///SIM		setXinj((double)xinj,(double)yinj,(double)zinj);
///SIM		Vair = new double[DIM];
///SIM		setVair(vdir);
///SIM	}
	public int getnnodes()
	{
		return nnodes;
	}
///SIM	public void settrapped(int n)
///SIM	{
///SIM		trapped=n;
///SIM	}
///SIM	public void addtrapped(int n)
///SIM	{
///SIM		trapped=n;
///SIM	}
///SIM	public int gettrapped()
///SIM	{
///SIM		return trapped;
///SIM	}
///SIM	public int[] getcolinj()
///SIM	{
///SIM		return Cinj;
///SIM	}
///SIM	static public double[] getXinj()
///SIM	{ return Xinj;
///SIM	}
///SIM	static public double getinj(int i)
///SIM	{	return Xinj[i];
///SIM	}
///SIM	static public double getvdir()
///SIM	{	return vdir;
///SIM	}
///SIM	public void setinj()
///SIM	{
///SIM		trapped=0;
///SIM		nvert=mvert;
///SIM		nnodes=mnodes;
///SIM		for(int i=inj0;i<nvert;i++)
///SIM		///for(int i=inj1;i<nvert;i++)
///SIM		{	int ip=DIM*i;
///SIM			for (int j=0;j<3;j++)
///SIM			{	vert[ip+j]=(float)Math.floor(DX[j]*Xinj[j]);
///SIM			}
///SIM		}
///SIM		int jnj0=DIM*inj0,jnj1=DIM*inj1; 
///SIM		for (int i=0;i<DIM;i++)
///SIM		{	vert[jnj1+i]=vert[jnj0+i]+vairscale*(float)Vair[i];
///SIM		}
///SIM	}
///SIM	public int runnodes()
///SIM	{	int trapped=0,
///SIM			i0=nvert-nnodes+2,
///SIM			ii[] = new int[DIM];
///SIM		float turb=2.0F, dt=0.5F;
///SIM		boolean outside=false;
///SIM		do
///SIM		{	outside=false;
///SIM			for(int i=i0;i<nvert;i++)
///SIM			{	int ip=DIM*i;
///SIM				for (int j=0;j<DIM;j++)
///SIM				{	int k=(int)((vert[ip+j])/dX[j]+0.5);
///SIM					if(k<0) {k=0; outside=true; break;}
///SIM					else if(k>=Dim[j]){k=Dim[j]-1;outside=true;break;}
///SIM					ii[j]=k;
///SIM				}
///SIM				///if(outside){trapped++; continue;}
///SIM				if(outside)
///SIM					if(nvert>delVert(i))
///SIM					{	nnodes--; i0=i; 
///SIM						break;
///SIM					}
///SIM					else 
///SIM					{ 	trapped++;
///SIM						outside=false;
///SIM						continue;
///SIM					}
///SIM				if((int)grid[Domain.I(ii[0],ii[1],ii[2])]!=0)
///SIM				{ trapped++; continue;
///SIM				}
///SIM				for (int j=0;j<DIM;j++)
///SIM				{	float vel=(float)Vair[j]+(float)turb*(-1.0F+2.0F*(float)Math.random());
///SIM					vert[ip+j]+=vel*dt;
///SIM				}
///SIM			}
///SIM		}	while(outside);
///SIM		return trapped;
///SIM	}
}
/** Monte-Carlo Node Placement */
class Nodes
{
	static int NNODES=0;//initial value
	static int	DIM,Dim[];
	static double[] X0,dX,DX;//Domain origin and dimensions
	Space space;
	Domain domain;
	float vert[];
	Atom atoms[];
	int tvert[],
		ZsortMap[],
		nvert, maxvert,
		///nnodes=NNODES, 
		trapped=0;
	double energy=0.0;
	Potential potential;
	static	Hashtable atomTable = new Hashtable();
	static	Atom defaultAtom;
	static
	{
		defaultAtom = new Atom(255, 100, 200);
	}
	Matrix3D mat;
	boolean transformed;
	
	float xmin, xmax, ymin, ymax, zmin, zmax;
	
	Nodes(double sigma, double eta)
	{	mat = new Matrix3D();
		mat.xrot(20);
		mat.yrot(30);
		potential = new Potential(sigma,eta);
	}
	static void Init(Space space)
	{
		DIM = space.getDim();
		atomTable=null;
		// Initialize atom table
		atomTable = new Hashtable();
		for(int icomp=0;icomp<space.vComponents.size();icomp++)
		{	Components	comp=(Components)space.vComponents.elementAt(icomp);
			atomTable.put
			(	comp.name, 
				new Atom
				(	comp.color.getRed(), 
					comp.color.getGreen(), 
					comp.color.getBlue()
				)
			);
		}
	}
	public void init 
	(	Space	space,
		Domain	domain,
		int nnodes
	) throws Exception
	{
		///nnodes=n;
		this.space=space;
		this.domain=domain;
		Dim = domain.getDim();
		dX = space.getdX(); 
		X0 = domain.getX0();
		DX = domain.getDX();
		vert=null;
		tvert=null;
		ZsortMap=null;
		nvert=0;
		///for(int icomp=1;icomp<space.vComponents.size();icomp++)
		int icomp=1;// add the first component at the beginning
		{	Components	comp=(Components)space.vComponents.elementAt(icomp);
			for(int inode=0;inode<nnodes;inode++)
			{	double
				x=X0[0]+0.5*DX[0],
				y=X0[1]+0.5*DX[1],
				z=X0[2]+0.5*DX[2];
				addVert(comp.name, (float) x, (float) y, (float) z);
			}
		}
		set();
	} 
	public void set()
	{
		for(int i=0;i<nvert;i++)
		{	int ip=DIM*i;
			for (int j=0;j<DIM;j++)
			{	vert[ip+j]=(float)(X0[j]+Math.random()*DX[j]);
			}
		}
		//Compute system energy
		energy=0.0;//total system energy (MC)
		for(int i=0;i<nvert;i++)
		{	int ip=DIM*i;
			//Compute the system potential:
			double p=0.0;//old and new particle potentials
			for(int j=0;j<i;j++)
			{	int jp=DIM*j;
				double r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double
						d=vert[jp+k]-vert[ip+k];
					r+=d*d;
				}
				r=Math.sqrt(r);//new distance to the neighbor
				p+=potential.value(r);
			}
			//Skip the current atom and continue with the next one:
			for(int j=i+1;j<nvert;j++)
			{	int jp=DIM*j;
				double r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double
						d=vert[jp+k]-vert[ip+k];
					r+=d*d;
				}
				r=Math.sqrt(r);//distance to the neighbor
				p+=potential.value(r);
			}
			energy+=p;
		}
		energy/=2.0*(double)getnvert();
		trapped=0;
	}
	double getstepsize()
	{	double
			lengthscale=potential.lengthscale(),
			stepsize=0.1*lengthscale; ///should be read in as a parameter
		return stepsize;
	}
	double mcmove()
	{	//MONTE_CARLO NODE PLACEMENT
		double energy=0.0,
			stepsize=getstepsize(),
			smallength=1.e-9*stepsize,
			x[] = new double[DIM], //new position vector
			h[] = new double[DIM]; //dimension indicator: 0 or 1.
		for(int i=0;i<DIM;i++)
			if(DX[i]>smallength)h[i]=1.0;
			else h[i]=0.0;
		for(int i=0;i<nvert;i++)
		{	int ip=DIM*i;
			//Make a trial step:
			boolean move_rejected;
			do//Check if it does not go outside of domain bounds:
			{	move_rejected=false;
				for(int k=0;k<DIM;k++)
				{	///x[k]=vert[ip+k]+lengthscale*(1.0-2.0*Math.random());
					double 
						d=stepsize*(1.0-2.0*Math.random())*h[k],
						r=vert[ip+k]+d;
				///	if(d<X0[k]-smallength||d>X0[k]+DX[k]+smallength)
					if(r<X0[k]||r>X0[k]+DX[k])
					{	move_rejected=true;
						break;
					}
					x[k]=r;
				}
			}	while (move_rejected);
			//Compute the change in system potential:
			double p0=0.0,p=0.0;//old and new particle potentials
			for(int j=0;j<i;j++)
			{	int jp=DIM*j;
				double r0=0.0,r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double
						d0=vert[jp+k]-vert[ip+k], 
						d=vert[jp+k]-x[k];
					r0+=d0*d0;
					r+=d*d;
				}
				r0=Math.sqrt(r0);//old distance to the neighbor
				p0+=potential.value(r0);
				r=Math.sqrt(r);//new distance to the neighbor
				p+=potential.value(r);
			}
			//Skip the current atom and continue with the next one:
			for(int j=i+1;j<nvert;j++)
			{	int jp=DIM*j;
				double r0=0.0,r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double
						d0=vert[jp+k]-vert[ip+k], 
						d=vert[jp+k]-x[k];
					r0+=d0*d0;
					r+=d*d;
				}
				r0=Math.sqrt(r0);//old distance to the neighbor
				p0+=potential.value(r0);
				r=Math.sqrt(r);//distance to the neighbor
				p+=potential.value(r);
			}
			if(p<p0)//accept the step if the potential decreased
			{	for(int k=0;k<DIM;k++)
					vert[ip+k]=(float)x[k];///DDD: Math.floor(DX[j]*Xinj[j]);
				energy+=p;
			}
			else
				energy+=p0;
		}
		return energy/2.0*(double)getnvert();
	}
	double mdmove()
	{
		double
			lengthscale=potential.lengthscale(),
			stepsize=getstepsize(), ///should be read in as a parameter
			energy=0.0,
			x[] = new double[DIM], //new position vector
			F[] = new double[DIM]; //force 
		for(int i=0;i<nvert;i++)
		{	int ip=DIM*i;
			//MOLECULAR DYNAMICS:
			//Compute forces from all particles:
			for(int k=0;k<DIM;k++) F[k]=0.0;
			//Loop over neighbors
			for(int j=0;j<i;j++)
			{	int jp=DIM*j;
				double r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double y=vert[jp+k]-vert[ip+k];
					r+=y*y;
					x[k]=y;
				}
				r=Math.sqrt(r);//distance to the neighbor
				for(int k=0;k<DIM;k++)
					F[k]+=x[k]/r*potential.value(r);//Force 
			}
			//Jump over this node and continue looping
			// over the neighbors (this is to avoid
			// the if-statement inside the loop)
			for(int j=i+1;j<nvert;j++)
			{	int jp=DIM*j;
				double r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double y=vert[jp+k]-vert[ip+k];
					r+=y*y;
					x[k]=y;
				}
				r=Math.sqrt(r);
				for(int k=0;k<DIM;k++)
					F[k]+=x[k]/r*potential.value(r);
			}
			//Move particle under the action of cumulative force:
			double
				relax=1.e-1, 
				f2=0.0,force=0.0,displacement=0.0;
			//Normalize the force
			for (int k=0;k<DIM;k++)
			{	double
					f=F[k];//displacement
				f2+=f*f;
			}
			force=Math.sqrt(f2);
			displacement=relax*lengthscale*force;
			if(displacement>stepsize)
			for(int k=0;k<DIM;k++)
				F[k]*=stepsize/force;
			else
			for(int k=0;k<DIM;k++)
				F[k]/=force;
			boolean move_accepted=true;
			for (int k=0;k<DIM;k++)
			{	double
					r=vert[ip+k]+F[k];
				//Check if the particle is outside of bounds
				if(r<X0[k]||r>X0[k]+DX[k])
				{	move_accepted=false;
					break;
				}
				x[k]=r;
			}
			if(move_accepted)
			for(int k=0;k<DIM;k++)
				vert[ip+k]=(float)x[k];
		}
		return energy/2.0*(double)getnvert();
	}
	double mdmcmove()
	{
		double
			stepsize=getstepsize(), ///should be read in as a parameter
			energy=0.0,
			x[] = new double[DIM], //new position vector
			F[] = new double[DIM]; //force
///			fdir[] = new double[DIM]; //force unit direction vector
		for(int i=0;i<nvert;i++)
		{	int ip=DIM*i;
			//MOLECULAR DYNAMICS/MONTE-CARLO: MDMC
			//Compute forces from all particles
			// as well as the change in system potential:
			for(int k=0;k<DIM;k++) F[k]=0.0;
			double p0=0.0;//old and new particle potential
			//Loop over neighbors
			for(int j=0;j<i;j++)
			{	int jp=DIM*j;
				double r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double y=vert[jp+k]-vert[ip+k];
					r+=y*y;
					x[k]=y;
				}
				r=Math.sqrt(r);//distance to the neighbor
				double p=potential.value(r);
				p0+=p;
				for(int k=0;k<DIM;k++)
					F[k]+=x[k]/r*p;//Force 
			}
			//Jump over this node and continue looping
			// over the neighbors (this is to avoid
			// the if-statement inside the loop)
			for(int j=i+1;j<nvert;j++)
			{	int jp=DIM*j;
				double r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double y=vert[jp+k]-vert[ip+k];
					r+=y*y;
					x[k]=y;
				}
				r=Math.sqrt(r);
				double p=potential.value(r);
				p0+=p;
				for(int k=0;k<DIM;k++)
					F[k]+=x[k]/r*p;
			}
			//Move particle under the action of cumulative force:
			double
				relax=1.e-1, 
				f2=0.0,force=0.0,
				displacement=0.0;
			//Normalize the force
			for (int k=0;k<DIM;k++)
			{	double
					f=F[k];//displacement
				f2+=f*f;
			}
			force=Math.sqrt(f2);
///			//Unit direction vector of the force:
///			for(int k=0;k<DIM;k++)fdir[k]=F[k]/force;
///			displacement=relax*lengthscale*force;
///			if(displacement>stepsize)
///			for(int k=0;k<DIM;k++)
///				F[k]*=stepsize/force;
///			else
///			for(int k=0;k<DIM;k++)
///				F[k]/=force;
			//Move the particle
			for (int k=0;k<DIM;k++)
			{	double
					r=vert[ip+k]-stepsize*Math.random()*F[k]/force;
				//Check if the particle is outside of bounds
				if(r<X0[k])r=X0[k];
				if(r>X0[k]+DX[k])r=X0[k]+DX[k];
				x[k]=r;
			}
			//Compute the change in system potential:
			double p=0.0; //new particle potential
			for(int j=0;j<i;j++)
			{	int jp=DIM*j;
				double r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double
						d=vert[jp+k]-x[k];
					r+=d*d;
				}
				r=Math.sqrt(r);//new distance to the neighbor
				p+=potential.value(r);
			}
			//Skip the current atom and continue with the next one:
			for(int j=i+1;j<nvert;j++)
			{	int jp=DIM*j;
				double r=0.0D;
				for(int k=0;k<DIM;k++)
				{	double
						d=vert[jp+k]-x[k];
					r+=d*d;
				}
				r=Math.sqrt(r);//distance to the neighbor
				p+=potential.value(r);
			}
			if(p<p0)//accept the step if the potential decreased
			{	for(int k=0;k<DIM;k++)
					vert[ip+k]=(float)x[k];
				energy+=p;
			}
			else
				energy+=p0;
			//if(Math.sqrt(d2)<lengthscale/10.0)trapped++;//does not work for large initial separations
		}
		return energy/2.0*(double)getnvert();
	}
	public void run()
	{
		switch(domain.model)
		{
		case Domain.MONTE_CARLO_MODEL:
			energy=mcmove();
			break;
		case Domain.GCMC_MODEL:
			energy=0.0;//total system energy (MC)
			if((energy=mcmove())<=0.0)
			{	//Add node
				int icomp=space.vComponents.size()-1;//select the last component
				Components	comp=(Components)space.vComponents.elementAt(icomp);
				double
					x0=X0[0]+Math.random()*DX[0],
					x1=X0[1]+Math.random()*DX[1],
					x2=X0[2]+Math.random()*DX[2];
				addVert(comp.name, (float) x0, (float) x1, (float) x2);
				ZsortMap=null;
			}
			break;//END GCMC_MODEL
		case Domain.MOLECULAR_DYNAMICS_MODEL: 
			energy=mdmove();
			break;
		case Domain.MDMC_MODEL: 
			energy=mdmcmove();
			break;
		case Domain.MDGCMC_MODEL:
			energy=0.0;//total system energy (MC)
			if((energy=mdmcmove())<=0.0)
			{	//Add node
				int icomp=space.vComponents.size()-1;//select the last component
				Components	comp=(Components)space.vComponents.elementAt(icomp);
				double
					x0=X0[0]+Math.random()*DX[0],
					x1=X0[1]+Math.random()*DX[1],
					x2=X0[2]+Math.random()*DX[2];
				addVert(comp.name, (float) x0, (float) x1, (float) x2);
				ZsortMap=null;
			}
			break;//END GCMC_MODEL
		}
///SIM		int jnj0=DIM*inj0,jnj1=DIM*inj1; 
///SIM		for (int i=0;i<DIM;i++)
///SIM		{	vert[jnj1+i]=vert[jnj0+i]+vairscale*(float)Vair[i];
///SIM		}
	}
	int getnvert()
	{	return nvert;
	}
	int gettrapped()
	{	return trapped;
	}
	double getEnergy()
	{	return energy;
	}
	/** Add a vertex to this model */
	int addVert(String name, float x, float y, float z) 
	{	int i = nvert;
		if (i >= maxvert)
		if (vert == null) 
		{	maxvert = 100;
			vert = new float[maxvert * 3];
			atoms = new Atom[maxvert];
		} else 
		{	maxvert *= 2;
			float nv[] = new float[maxvert * 3];
			System.arraycopy(vert, 0, nv, 0, vert.length);
			vert = nv;
			Atom na[] = new Atom[maxvert];
			System.arraycopy(atoms, 0, na, 0, atoms.length);
			atoms = na;
		}
		Atom a = (Atom) atomTable.get(name);
		if (a == null) a = defaultAtom;
		atoms[i] = a;
		i *= 3;
		vert[i] = x;
		vert[i + 1] = y;
		vert[i + 2] = z;
		return nvert++;
	}
	/** Delete a vertex from this model */
	int delVert(int ivert) 
	{
		if (ivert >= nvert) return nvert;
		nvert--;
		for (int i=ivert;i<nvert;i++)
		{	int j = 3*i, j1=j+3;
			for(int k=0;k<3;k++)
				vert[j+k] = vert[j1+k];
			atoms[i]=atoms[i+1];
		}
		return nvert;
	}
	/** Transform all the points in this model */
	void transform() 
	{	if (transformed || nvert <= 0) return;
		if (tvert == null || tvert.length < nvert * 3)
		tvert = new int[nvert * 3];
		mat.transform(vert, tvert, nvert);
		transformed = true;
	}
	/** Paint this model to a graphics context.  It uses the matrix associated
	with this model to map from model space to screen space.
	The next version of the browser should have double buffering,
	which will make this *much* nicer */
	void paint(Graphics g) 
	{
		if (vert == null || nvert <= 0) return;
		transform();
		int v[] = tvert;
		int zs[] = ZsortMap;
		if (zs == null) 
		{
			ZsortMap = zs = new int[nvert];
			for (int i = nvert; --i >= 0;)
			zs[i] = i * 3;
		}
		/*
		* I use a bubble sort since from one iteration to the next, the sort
		* order is pretty stable, so I just use what I had last time as a
		* "guess" of the sorted order.  With luck, this reduces O(N log N)
		* to O(N)
		*/
		for (int i = nvert - 1; --i >= 0;) 
		{	boolean flipped = false;
			for (int j = 0; j <= i; j++) 
			{	int a = zs[j];
				int b = zs[j + 1];
				if (v[a + 2] > v[b + 2]) 
				{
					zs[j + 1] = a;
					zs[j] = b;
					flipped = true;
				}
			}
			if (!flipped)
			break;
		}
		int lg = 0;
		int lim = nvert;
		Atom ls[] = atoms;
		if (lim <= 0 || nvert <= 0)return;
		for (int i = 0; i < lim; i++) 
		{	int j = zs[i];
			int grey = v[j + 2];
			if (grey < 0)
			grey = 0;
			if (grey > 15)
			grey = 15;
			// g.drawString(names[i], v[j], v[j+1]);
			atoms[j/3].paint(g, v[j], v[j + 1], grey);
			// g.drawImage(iBall, v[j] - (iBall.width >> 1), v[j + 1] -
			// (iBall.height >> 1));
		}
	}
	/** Find the bounding box of this model */
	void findBB() 
	{	if (nvert <= 0)
		return;
///		float v[] = vert;
///		float xmin = v[0], xmax = xmin;
///		float ymin = v[1], ymax = ymin;
///		float zmin = v[2], zmax = zmin;
///		for (int i = nvert * 3; (i -= 3) > 0;) 
///		{	float x = v[i];
///			if (x < xmin)
///			xmin = x;
///			if (x > xmax)
///			xmax = x;
///			float y = v[i + 1];
///			if (y < ymin)
///			ymin = y;
///			if (y > ymax)
///			ymax = y;
///			float z = v[i + 2];
///			if (z < zmin)
///			zmin = z;
///			if (z > zmax)
///			zmax = z;
///		}
///		this.xmax = xmax;
///		this.xmin = xmin;
///		this.ymax = ymax;
///		this.ymin = ymin;
///		this.zmax = zmax;
///		this.zmin = zmin;
		this.xmax = (float)(domain.X0[0]+domain.DX[0]);
		this.xmin = (float)domain.X0[0];
		this.ymax = (float)(domain.X0[1]+domain.DX[1]);;
		this.ymin = (float)domain.X0[1];
		this.zmax = (float)(domain.X0[2]+domain.DX[2]);;
		this.zmin = (float)domain.X0[2];
	}
	class Potential
	{
		public double sigma=1.0, eta=1.0;
		Potential(double sigma, double eta)
		{	this.sigma=sigma;
			this.eta=eta;
		}
		public double value(double r)
		{
			return eta*(Math.pow(sigma/r,12.0) - Math.pow(sigma/r,6));
		}
		public double lengthscale()
		{
			return sigma;
		}
		public double strength()
		{
			return eta;
		}
	}
}

/** The atomized suface */

/** The representation of a 3D suface */
//+class BoundarySurface
//+{
//+	float vert[];
//+	Atom atoms[];
//+	int tvert[];
//+	int ZsortMap[];
//+	int nvert, maxvert;
//+	
//+	static Hashtable atomTable = new Hashtable();
//+	static Atom defaultAtom;
//+	static 
//+	{
//+		defaultAtom = new Atom(255, 100, 200);
//+	}
//+	boolean transformed;
//+	Matrix3D mat;
//+	
//+	float xmin, xmax, ymin, ymax, zmin, zmax;
//+	
//+	BoundarySurface() 
//+	{
//+		mat = new Matrix3D();
//+		mat.xrot(20);
//+		mat.yrot(30);
//+	}
//+	// Create Atom table
//+	static public void initAtoms (Space space) 
//+	{
//+		atomTable=null;
//+		// Initialize atom table
//+		atomTable = new Hashtable();
//+		for(int icomp=0;icomp<space.vComponents.size();icomp++)
//+		{	Components	comp=(Components)space.vComponents.elementAt(icomp);
//+			atomTable.put
//+			(	comp.name, 
//+				new Atom
//+				(	comp.color.getRed(), 
//+					comp.color.getGreen(), 
//+					comp.color.getBlue()
//+				)
//+			);
//+		}
//+	}
//+	/** Create a boundary surface */
//+	public void init 
//+	(	Space	space,
//+		Domain	domain
//+	) throws Exception
//+	{
//+		double
//+			dx[] = space.getdX(),
//+			x0[] = domain.getX0();
//+		vert=null;
//+		tvert=null;
//+		ZsortMap=null;
//+		nvert=0;
//+		for(int icomp=1;icomp<space.vComponents.size();icomp++)
//+		{	Components	comp=(Components)space.vComponents.elementAt(icomp);
//+			if(comp.boundaryCell!=null)
//+			{	int[] ind=comp.boundaryCell.getInd();
//+				for
//+				(	BoundaryCell curr=comp.boundaryCell;
//+					curr!=null; curr=curr.next
//+				)
//+				{	double
//+						x=x0[0]+curr.ind[0]*dx[0],
//+						y=x0[1]+curr.ind[1]*dx[1],
//+						z=x0[2]+curr.ind[2]*dx[2];
//+					addVert(comp.name, (float) x, (float) y, (float) z);
//+				}
//+			}
//+			else
//+				System.out.println("No boundary for component "+icomp);
//+		}
//+	}  // end BoundarySurface ()
//+	/** Add a vertex to this model */
//+	int addVert(String name, float x, float y, float z) 
//+	{	int i = nvert;
//+		if (i >= maxvert)
//+		if (vert == null) 
//+		{	maxvert = 100;
//+			vert = new float[maxvert * 3];
//+			atoms = new Atom[maxvert];
//+		} else 
//+		{	maxvert *= 2;
//+			float nv[] = new float[maxvert * 3];
//+			System.arraycopy(vert, 0, nv, 0, vert.length);
//+			vert = nv;
//+			Atom na[] = new Atom[maxvert];
//+			System.arraycopy(atoms, 0, na, 0, atoms.length);
//+			atoms = na;
//+		}
//+		Atom a = (Atom) atomTable.get(name);
//+		if (a == null) a = defaultAtom;
//+		atoms[i] = a;
//+		i *= 3;
//+		vert[i] = x;
//+		vert[i + 1] = y;
//+		vert[i + 2] = z;
//+		return nvert++;
//+	}
//+	/** Transform all the points in this model */
//+	void transform() 
//+	{	if (transformed || nvert <= 0) return;
//+		if (tvert == null || tvert.length < nvert * 3)
//+		tvert = new int[nvert * 3];
//+		mat.transform(vert, tvert, nvert);
//+		transformed = true;
//+	}
//+	/** Paint this model to a graphics context.  It uses the matrix associated
//+	with this model to map from model space to screen space.
//+	The next version of the browser should have double buffering,
//+	which will make this *much* nicer */
//+	void paint(Graphics g) 
//+	{
//+		if (vert == null || nvert <= 0) return;
//+		transform();
//+		int v[] = tvert;
//+		int zs[] = ZsortMap;
//+		if (zs == null) 
//+		{
//+			ZsortMap = zs = new int[nvert];
//+			for (int i = nvert; --i >= 0;)
//+			zs[i] = i * 3;
//+		}
//+		/*
//+		* I use a bubble sort since from one iteration to the next, the sort
//+		* order is pretty stable, so I just use what I had last time as a
//+		* "guess" of the sorted order.  With luck, this reduces O(N log N)
//+		* to O(N)
//+		*/
//+		for (int i = nvert - 1; --i >= 0;) 
//+		{	boolean flipped = false;
//+			for (int j = 0; j <= i; j++) 
//+			{	int a = zs[j];
//+				int b = zs[j + 1];
//+				if (v[a + 2] > v[b + 2]) 
//+				{
//+					zs[j + 1] = a;
//+					zs[j] = b;
//+					flipped = true;
//+				}
//+			}
//+			if (!flipped)
//+			break;
//+		}
//+		int lg = 0;
//+		int lim = nvert;
//+		Atom ls[] = atoms;
//+		if (lim <= 0 || nvert <= 0)return;
//+		for (int i = 0; i < lim; i++) 
//+		{	int j = zs[i];
//+			int grey = v[j + 2];
//+			if (grey < 0)
//+			grey = 0;
//+			if (grey > 15)
//+			grey = 15;
//+			// g.drawString(names[i], v[j], v[j+1]);
//+			atoms[j/3].paint(g, v[j], v[j + 1], grey);
//+			// g.drawImage(iBall, v[j] - (iBall.width >> 1), v[j + 1] -
//+			// (iBall.height >> 1));
//+		}
//+	}
//+	/** Find the bounding box of this model */
//+	void findBB() 
//+	{	if (nvert <= 0)
//+		return;
//+		float v[] = vert;
//+		float xmin = v[0], xmax = xmin;
//+		float ymin = v[1], ymax = ymin;
//+		float zmin = v[2], zmax = zmin;
//+		for (int i = nvert * 3; (i -= 3) > 0;) 
//+		{	float x = v[i];
//+			if (x < xmin)
//+			xmin = x;
//+			if (x > xmax)
//+			xmax = x;
//+			float y = v[i + 1];
//+			if (y < ymin)
//+			ymin = y;
//+			if (y > ymax)
//+			ymax = y;
//+			float z = v[i + 2];
//+			if (z < zmin)
//+			zmin = z;
//+			if (z > zmax)
//+			zmax = z;
//+		}
//+		this.xmax = xmax;
//+		this.xmin = xmin;
//+		this.ymax = ymax;
//+		this.ymin = ymin;
//+		this.zmax = zmax;
//+		this.zmin = zmin;
//+	}
//+}
/** An applet to put a model into a page */
public class Mesh
//extends JPanel
extends Applet
implements Runnable, MouseListener, MouseMotionListener //-, ItemListener
{
//+	static final int MODE_WIREFRAME = 0;
//+	static final int MODE_ATOMS = 1;
//+	int mode = MODE_WIREFRAME;

	static Space	space;
	static Domain	domain=null;
//-	boolean XON=true, YON=true, ZON=true;
//-	BoundarySurface md;
	boolean painted = true;
	boolean doRun = false;       // True if thread currently paused by user
	float xfac;
	int prevx, prevy;
	float xtheta, ytheta;
	float scalefudge = 1;
	Matrix3D amat = new Matrix3D(), tmat = new Matrix3D();
	String mdname = null;
	String message = null;
	Image backBuffer;
	Graphics backGC = null;
	Dimension backSize;
	Nodes nodes;
	WireFrame wireframe;
///	JComboBox comboBox3DRenderModes = new JComboBox();

	float deltaX = 0.0f, deltaY = 0.0f, scale = 1.0f;
	boolean left = false, center = false, right = false;
    
	JPanel drawingSurface = new JPanel();
///    JButton buttonRotX = new JButton("Rotate X"), buttonRotY = new JButton("Rotate Y");
///	JScrollBar scrollBarX,scrollBarY,scrollBarZ,scrollBarA;
	JLabel labelEnergy = new JLabel("Energy=0.0000e+00");
	Choice choiceModel = new Choice();
	JTextField Nnodes;
	JButton buttonRun = new JButton("Run   ");
	JButton buttonSet = new JButton("Set");
//-	Checkbox checkboxX = new Checkbox("X");
//-	Checkbox checkboxY = new Checkbox("Y");
//-	Checkbox checkboxZ = new Checkbox("Z");

	private synchronized void newBackBuffer() 
	{
		backBuffer = createImage(getSize().width, getSize().height);
		if (backGC != null) 
		{
			backGC.dispose();
		}
		backGC = backBuffer.getGraphics();
		backSize = getSize();
	}
	public void init() 
	{
///		comboBox3DRenderModes.addItem("Wireframe View");
///		comboBox3DRenderModes.addItem("Atom View");
///		comboBox3DRenderModes.addActionListener(new ChangeRenderMode());
///		mode = MODE_WIREFRAME;
		int 
			model=0,
			nnodes=0;//number of injected particles
//SIM			rinj=255, //particles default colors
//SIM			ginj=0,
//SIM			binj=0;
		float 
			xrot0=0.0F,
			yrot0=0.0F,
			zrot0=0.0F;
//SIM			xinj=0.0F,
//SIM			yinj=0.0F,
//SIM			zinj=0.0F,
//SIM			vdir=0.0F;
		double sigma=1.0, eta=1.0;//potential constants
		String str=null;
		mdname = getParameter("model");
		if (mdname == null)
			mdname = "project.config";
		try
		{	str = getParameter("xrot");
			if (str != null) 
				xrot0 = (Float.valueOf(str)).floatValue() ;
			str = getParameter("yrot");
			if (str != null) 
				yrot0 = (Float.valueOf(str)).floatValue() ;
			str = getParameter("zrot");
			if (str != null) 
				zrot0 = (Float.valueOf(str)).floatValue() ;
//SIM			str = getParameter("nnodes");
//SIM			if (str != null) 
//SIM				nnodes = (Integer.valueOf(str)).intValue() ;
//SIM			str = getParameter("xinj");
//SIM			if (str != null) 
//SIM				xinj = (Float.valueOf(str)).floatValue() ;
//SIM			if(xinj<0.0)xinj=0.0F; if(xinj>1.0)xinj=1.0F;
//SIM			str = getParameter("yinj");
//SIM			if (str != null) 
//SIM				yinj = (Float.valueOf(str)).floatValue() ;
//SIM			if(yinj<0.0)yinj=0.0F; if(yinj>1.0)yinj=1.0F;
//SIM			str = getParameter("zinj");
//SIM			if (str != null) 
//SIM				zinj = (Float.valueOf(str)).floatValue() ;
//SIM			if(zinj<0.0)zinj=0.0F; if(zinj>1.0)zinj=1.0F;
//SIM			str = getParameter("rinj");
//SIM			if (str != null) 
//SIM				rinj = (Integer.valueOf(str)).intValue() ;
//SIM			if(rinj<0)rinj=0; if(rinj>255)rinj=255;
//SIM			str = getParameter("ginj");
//SIM			if (str != null) 
//SIM				ginj = (Integer.valueOf(str)).intValue() ;
//SIM			if(ginj<0)ginj=0; if(ginj>255)ginj=255;
//SIM			str = getParameter("binj");
//SIM			if (str != null) 
//SIM				binj = (Integer.valueOf(str)).intValue() ;
//SIM			if(binj<0)binj=0; if(binj>255)binj=255;
//SIM			str = getParameter("vdir");
//SIM			if (str != null) 
//SIM				vdir = (Float.valueOf(str)).floatValue() ;
		}
		catch (NumberFormatException e) 
		{	System.out.println("Failed to read from "+mdname+": "+e.toString());
		}
		String	projectDirectory = "./";
		space= new Space("Demo",1.0,1.0,1.0);
		URL url = this.getCodeBase();
		String url_str = url.toString();
///		System.out.println("URL: "+url_str);
		space.setPath(projectDirectory);
		space.vComponents.clear();
		Components empty = new Components();
		empty.name = "This is an empty component.  It's here only to make things compatible w/the master.";
		empty.color = new Color(255,255,255);
		space.vComponents.add(empty);
		try
		{
			InputStream is = new URL(getDocumentBase(),mdname).openStream();
			StreamTokenizer in = new StreamTokenizer
			(new BufferedReader(new InputStreamReader(is)));
			//Read in the components
			in.nextToken(); //NumComponents
			in.nextToken();
			int numComponents = (int)in.nval;
			for ( int i = 0; i < numComponents; i++ )
			{	Components newDCI = new Components();
				in.nextToken();
				newDCI.name = in.sval;
				in.nextToken();
				int r = (int)in.nval;
				in.nextToken();
				int g = (int)in.nval;
				in.nextToken();
				int b = (int)in.nval;
				newDCI.color = new Color(r, g, b);
				space.vComponents.add(newDCI);
			}
			if (space.vComponents.size() > 1)
			//Read in the domain names and load them as "ghosts"
			in.nextToken(); //NumDomains
			if(!"NumDomains".equals(in.sval))throw new Exception("Was expecting \"NumDomains\"");
			in.nextToken();
			int numDomains = (int)in.nval;///DDD: this should be removed
			                             /// we don't need the domain number
			//Empty all the current domains/connections
			space.RemoveDomains();
			Nodes.Init(space);
			Domain dom;///DDD = new Domain("dummy", 0, 0, 0, space);
			for (int index = 0; index < numDomains; index++)
			{	String name;
				/*
				File f = new File
				(projectDirectory + File.separator + in.sval + ".dom");
				*/
				///if (f.exists()) 
				///{
				in.nextToken(); //"Model"
				if(!"Model".equals(in.sval))throw new Exception("Was expecting \"Model\"");
				in.nextToken();
				model=(int)in.nval;
				in.nextToken();
				name = in.sval;
				{	int NX, NY, NZ;
					double X0,Y0,Z0,DX,DY,DZ,dX,dY,dZ;
				// rEAD Dim[], X0[i], DX[DIM] from config file
					in.nextToken(); //"DIM"
					if(!"DIM".equals(in.sval))throw new Exception("Was expecting \"DIM\"");
					in.nextToken();
					NX=(int)in.nval;
					in.nextToken();
					NY=(int)in.nval;
					in.nextToken();
					NZ=(int)in.nval;
					in.nextToken(); //"X0"
					if(!"X0".equals(in.sval))throw new Exception("Was expecting \"X0\"");
					in.nextToken();
					X0=(double)in.nval;
					in.nextToken();
					Y0=(double)in.nval;
					in.nextToken();
					Z0=(double)in.nval;
///					in.nextToken(); //"DX"
///					if(!"DX".equals(in.sval))throw new Exception("Was expecting \"DX\"");
///					in.nextToken();
///					DX=(double)in.nval;
///					in.nextToken();
///					DY=(double)in.nval;
///					in.nextToken();
///					DZ=(double)in.nval;
					in.nextToken(); //"DX"
					if(!"DX".equals(in.sval))throw new Exception("Was expecting \"DX\"");
					in.nextToken();
//-					DX=XON ? (double)in.nval : 0.0;
					DX=(double)in.nval;
					in.nextToken();
//-					DY=YON ? (double)in.nval : 0.0;
					DY=(double)in.nval;
					in.nextToken();
//-					DZ=ZON ? (double)in.nval : 0.0;
					DZ=(double)in.nval;
					dX=DX/(double)NX;
					dY=DY/(double)NY;
					dZ=DZ/(double)NZ;
					space.setdX(dX,dY,dZ);
					dom = new Domain(name,0,0,0,NX,NY,NZ,space);
					dom.setX0(X0,Y0,Z0);
					dom.setModel(model);
					if (space.first_domain == null)
						space.first_domain = dom;
					else
						space.first_domain.Append(dom);
					in.nextToken(); //"NumParticles"
					if(!"NumParticles".equals(in.sval))throw new Exception("Was expecting \"NumParticles\"");
					in.nextToken(); //"NumParticles"
					nnodes=(int)in.nval;
					in.nextToken(); //"Potential"
					if(!"Potential".equals(in.sval))throw new Exception("Was expecting \"Potential\"");
					in.nextToken(); //"sigma"
					if(!"sigma".equals(in.sval))throw new Exception("Was expecting \"sigma\"");
					in.nextToken();
					sigma=(double)in.nval;
					in.nextToken(); //"eta"
					if(!"eta".equals(in.sval))throw new Exception("Was expecting \"eta\"");
					in.nextToken();
					eta=(double)in.nval; //constants for a potential
				}
///			}
				/*
				else
				JOptionPane.showMessageDialog
				(	null,
					new JLabel("Domain " + in.sval + " doesn't exist.  Did you delete it?"), 
					"Warning!",
					JOptionPane.WARNING_MESSAGE
				);
				*/
			}
			setLayout(new BorderLayout());
			drawingSurface.setVisible(false);
			add(drawingSurface);
			JPanel north = new JPanel();
//-			north.add(checkboxX); checkboxX.addItemListener(this);
//-			north.add(checkboxY); checkboxY.addItemListener(this);
//-			north.add(checkboxZ); checkboxZ.addItemListener(this);
			north.add(labelEnergy);
			add(north, BorderLayout.NORTH);

			JPanel south = new JPanel();
///			buttonRotX.addActionListener(new Rotate());
///			buttonRotY.addActionListener(new Rotate());
///			south.add(comboBox3DRenderModes);
///			south.add(buttonRotX); 
///			south.add(buttonRotY);

//SIM			scrollBarX = new JScrollBar
//SIM			(	JScrollBar.VERTICAL, 1, 1, 1, 
//SIM				dom.Dim[0]-1
//SIM			);
//SIM			scrollBarX.addAdjustmentListener(new ChangeXinj());
//SIM			scrollBarX.setToolTipText("Injection position X");
//SIM			south.add(scrollBarX); 
//SIM
//SIM			scrollBarY = new JScrollBar
//SIM			(	JScrollBar.VERTICAL, 1, 1, 1, 
//SIM				dom.Dim[1]-1
//SIM			);
//SIM			scrollBarY.addAdjustmentListener(new ChangeYinj());
//SIM			scrollBarY.setToolTipText("Injection position Y");
//SIM			south.add(scrollBarY); 
//SIM
//SIM			scrollBarZ = new JScrollBar
//SIM			(	JScrollBar.VERTICAL, 1, 1, 1, 
//SIM				dom.Dim[2]-1
//SIM			);
//SIM			scrollBarZ.addAdjustmentListener(new ChangeZinj());
//SIM			scrollBarZ.setToolTipText("Injection position Z");
//SIM			south.add(scrollBarZ); 
//SIM
//SIM			scrollBarA = new JScrollBar
//SIM			(	JScrollBar.VERTICAL, 1, 1, 1, 
//SIM				360
//SIM			);
//SIM			scrollBarA.addAdjustmentListener(new ChangeVair());
//SIM			scrollBarA.setToolTipText("Wind direction angle");
//SIM			south.add(scrollBarA); 

			south.add(choiceModel);
			choiceModel.addItem("MC");
			choiceModel.addItem("GCMC");
			choiceModel.addItem("MD");
			choiceModel.addItem("MDMC");
			choiceModel.addItem("MDGCMC");
			choiceModel.select(Domain.MONTE_CARLO_MODEL);

			south.add(Nnodes = new JTextField(Integer.toString(nnodes), 4));
			Nnodes.addKeyListener
			(	new KeyAdapter() 
				{	public void keyTyped(KeyEvent e) 
					{	char c = e.getKeyChar();      
						if (!((Character.isDigit(c)) ||
							(c == KeyEvent.VK_BACK_SPACE) ||
							(c == KeyEvent.VK_DELETE))) 
						{
							getToolkit().beep();
							e.consume();
						}
						if(c==KeyEvent.VK_ENTER)
						{	int np=nodes.getnvert();	
							String s = Nnodes.getText();
							int nnodes=0;
							if(s.equals(""))
							{
								Nnodes.setText("0");
								nnodes = 0;
							}
							else
								nnodes = Integer.parseInt(s);
							if(np!=nnodes)
							{	try
								{	double 
										sigma=nodes.potential.lengthscale(),
										eta=nodes.potential.strength();
									nodes = new Nodes(sigma,eta);
									nodes.init(space,domain,nnodes);
									nodes.findBB();
									float xw = nodes.xmax - nodes.xmin;
									float yw = nodes.ymax - nodes.ymin;
									float zw = nodes.zmax - nodes.zmin;
									if (yw > xw)
									xw = yw;
									if (zw > xw)
									xw = zw;
									float f1 = getSize().width / xw;
									float f2 = getSize().height / xw;
									xfac = 0.7f * (f1 < f2 ? f1 : f2) * scalefudge;
								} catch(Exception e1) 
								{	e1.printStackTrace();
									nodes = null;
									message = e1.toString();
									System.out.println(message);
								}
								nodes.set();
								doRun=false;
								buttonRun.setText("Run   ");
								labelEnergy.setText("Energy: "+nodes.getEnergy());
								repaint();
							}
						}
					}
				}
			);

			buttonSet.addActionListener(new Set());
			south.add(buttonSet); 
			buttonRun.addActionListener(new Run());
			south.add(buttonRun); 

			add(south, BorderLayout.SOUTH);
		}
		catch(Exception e)
		{
			System.out.println("exception in init():"+ e);
		}
		try 
		{	str = getParameter("scale");
			if(str!=null)scalefudge = Float.valueOf(str).floatValue();
		}	catch(Exception e) 
		{
			System.out.println("exception in init() valueOf:"+ e);
		}
		amat.xrot(xrot0);
		amat.yrot(yrot0);
		amat.zrot(zrot0);
		resize(getSize().width <= 300 ? 480 : getSize().width,
		getSize().height <= 300 ? 480 : getSize().height);
		newBackBuffer();
		addMouseListener(this);
		addMouseMotionListener(this);
		domain= space.first_domain;
//SIM		domain.Load(url_str);
//SIM		domain.setBoundaryCells
//SIM		(
//SIM			space.vComponents
//SIM		);
//SIM		for 
//SIM		(	int plane_orientation = 0;
//SIM			plane_orientation < space.getDim();
//SIM			plane_orientation++
//SIM		)
//SIM		{	int separation=1;///DDD should be user-defined via gui.
//SIM			domain.setBoundaryContours
//SIM			(
//SIM				plane_orientation, ///DDD: should be user-defined
//SIM				separation, ///DDD: should be user-defined
//SIM				space.vComponents
//SIM			);
//SIM		}
//+		BoundarySurface.initAtoms(space);
//SIM		WireFrame.initinj(nnodes,rinj,ginj,binj,xinj,yinj,zinj,vdir);
//SIM		wireframe = new WireFrame (space, domain, space.vComponents);
//SIM		scrollBarA.setValue((int)(WireFrame.getvdir()));
//SIM		scrollBarX.setValue((int)(WireFrame.getinj(0)*(double)domain.Dim[0]));
//SIM		scrollBarY.setValue((int)(WireFrame.getinj(1)*(double)domain.Dim[1]));
//SIM		scrollBarZ.setValue((int)(WireFrame.getinj(2)*(double)domain.Dim[2]));
		validate();

/// MOVED FROM run():
//SIM		try //Initialize the WireFrame:
//SIM		{	Thread.currentThread().setPriority(Thread.MIN_PRIORITY);
//SIM		//if (domain==null)break;
//SIM///			WireFrame m = new WireFrame (space, domain, space.vComponents);
//SIM///			wireframe = m;
//SIM			WireFrame m = wireframe;
//SIM			m.findBB();
//SIM///QS			m.compress();
//SIM			float xw = m.xmax - m.xmin;
//SIM			float yw = m.ymax - m.ymin;
//SIM			float zw = m.zmax - m.zmin;
//SIM			if (yw > xw)
//SIM			xw = yw;
//SIM			if (zw > xw)
//SIM			xw = zw;
//SIM			float f1 = getSize().width / xw;
//SIM			float f2 = getSize().height / xw;
//SIM			xfac = 0.7f * (f1 < f2 ? f1 : f2) * scalefudge;
//SIM
//SIM		} 
//SIM		catch(Exception e) 
//SIM		{
//SIM			wireframe = null;
//SIM			message = e.toString();
//SIM		}

		try //Initialize Nodes
		{	Thread.currentThread().setPriority(Thread.MIN_PRIORITY);
///			BoundarySurface m = new BoundarySurface ();
///			if(domain!=null)m.init(space,domain);
///			else BoundarySurface.initAtoms(space);
			Nodes m = new Nodes(sigma,eta);
			m.init(space,domain,nnodes);
//			Atom.setPanel(this); // APPLICATION
			Atom.setApplet(this); // APPLET
			nodes = m;
			m.findBB();
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
		} catch(Exception e) 
		{	e.printStackTrace();
			nodes = null;
			message = e.toString();
		}
/// MOVED FROM run
		nodes.set();
		labelEnergy.setText("Energy: "+nodes.getEnergy());
		repaint();
	}
	static public void init(Space space, Domain D)
	{
		domain=D;
//+		BoundarySurface.initAtoms(space);
	}
	public void destroy() 
	{
		removeMouseListener(this);
		removeMouseMotionListener(this);
	}
	public void run() 
	{
//SIM		int
//SIM			np=wireframe.getnnodes()-2,//number of free particles (excludes the injection point)
//SIM			trapped=wireframe.gettrapped();//number of trapped particles
//SIM		do
//SIM		{
//SIM			if(doRun)	
//SIM			{
//SIM				trapped=wireframe.runnodes();
//SIM				np=wireframe.getnnodes()-2; 
//SIM				repaint();
//SIM			}
//SIM		}	while(doRun&&trapped<np);
		int
			iter=0,//iteration counter
			np=nodes.getnvert(),//number of free particles (excludes the injection point)
			trapped=nodes.gettrapped();//DDD: number of trapped particles
		if(doRun)
		{	if(trapped<np)	
			{	do
				{	domain.model=choiceModel.getSelectedIndex();
					nodes.run();
///System.out.println(""+iter+"\t"+nodes.getEnergy());iter++;///DDD
///DDD					trapped=nodes.gettrapped();
					if(nodes.getnvert()!=np)
					{	//Number of nodes has changed
						np=nodes.getnvert();
					}
					labelEnergy.setText("Energy: "+nodes.getEnergy());
					Nnodes.setText(""+nodes.getnvert());
					repaint();
				}	while(doRun&&trapped<np);
				doRun=false;
			}
///			else
///			{	nodes.set();
///				doRun=false;
///				buttonRun.setText("Run   ");
///				repaint();
///			}
		}
		if(trapped>=np)
		{	buttonRun.setText("End   ");
			repaint();
		}
	}
	public void start() 
	{
//-	   if (md == null && message == null)
	   if (message == null)
			new Thread(this).start();
	}
	public void stop() 
	{
	}
	/* event handling */
	public void mouseClicked(MouseEvent e) 
	{
	}
	public void mousePressed(MouseEvent e) 
	{
		if (e.getButton() == e.BUTTON1)
		left = true;
		else if (e.getButton() == e.BUTTON2)
		center = true;
		else if (e.getButton() == e.BUTTON3)
		right = true;
		prevx = e.getX();
		prevy = e.getY();
		e.consume();
	}
	public void mouseReleased(MouseEvent e) 
	{
		if (e.getButton() == e.BUTTON1)
		left = false;
		else if (e.getButton() == e.BUTTON2)
		center = false;
		else if (e.getButton() == e.BUTTON3)
		right = false;
	}
	public void mouseEntered(MouseEvent e) 
	{
	}
	public void mouseExited(MouseEvent e) 
	{
	}
	public void mouseDragged(MouseEvent e) 
	{	int x = e.getX();
		int y = e.getY();

		if (center || (left && right))
		{
			deltaX += -(prevx - x);
			deltaY += -(prevy - y);
		}
		else if (left)
		{
			tmat.unit();
			float xtheta = (prevy - y) * (360.0f / getSize().width);
			float ytheta = (x - prevx) * (360.0f / getSize().height);
			tmat.xrot(xtheta);
			tmat.yrot(ytheta);
			amat.mult(tmat);
		}
		else if (right)
		{
			scale += (prevy - y)/180.0f;
			if (scale < 0.1f)
			    scale = 0.1f;
		}
		if (painted) 
		{	painted = false;
			repaint();
		}
		prevx = x;
		prevy = y;
		e.consume();
	}
	public void mouseMoved(MouseEvent e) 
	{
	}
	public void update(Graphics g) 
	{	if (backBuffer == null)
		g.clearRect(0, 0, getSize().width, getSize().height);
		paint(g);
	}
	public void paint(Graphics g) 
	{
	    //g.clearRect(0,0,getSize().width, getSize().height); //This is ENTIRELY too slow

//SIM		switch(mode)
//SIM		{
//SIM		case MODE_ATOMS:
		if (nodes != null) 
		{	nodes.mat.unit();
			nodes.mat.translate(-(nodes.xmin + nodes.xmax) / 2,
			-(nodes.ymin + nodes.ymax) / 2,
			-(nodes.zmin + nodes.zmax) / 2);
			nodes.mat.mult(amat);
			nodes.mat.scale(xfac, -xfac, 16 * xfac / getSize().width);
			nodes.mat.scale(scale,scale,scale);
			nodes.mat.translate(getSize().width / 2, getSize().height / 2, 8);
///DDD			nodes.mat.translate(getSize().width / 2, getSize().height / 2, 4);
			nodes.mat.translate(deltaX, deltaY, 0.0f);
			nodes.transformed = false;
			if (backBuffer != null) 
			{	if (!backSize.equals(getSize()))
				newBackBuffer();
				backGC.setColor(getBackground());
				backGC.fillRect(0,0,getSize().width,getSize().height);
				nodes.paint(backGC);
				g.drawImage(backBuffer, 0, 0, this);
			}
			else
				nodes.paint(g);
			setPainted();
		} else if (message != null) 
		{	g.drawString("Error in model:", 3, 20);
			g.drawString(message, 10, 40);
		}
//SIM		break;
//SIM		case MODE_WIREFRAME: // WireFrame
		if (false) ///DDD (wireframe != null) 
		{	//Paint scene
			wireframe.mat.unit();
			wireframe.mat.translate(-(wireframe.xmin + wireframe.xmax) / 2,
					-(wireframe.ymin + wireframe.ymax) / 2,
					-(wireframe.zmin + wireframe.zmax) / 2);
			wireframe.mat.mult(amat);
			wireframe.mat.scale(xfac, -xfac, 16 * xfac / getSize().width);
			wireframe.mat.scale(scale,scale,scale);
			wireframe.mat.translate(getSize().width / 2, getSize().height / 2, 8);
			wireframe.mat.translate(deltaX, deltaY, 0.0f);
			wireframe.transformed = false;
			if (backBuffer != null)
			{
				if (!backSize.equals(getSize()))
					newBackBuffer();
				backGC.setColor(getBackground());
				backGC.fillRect(0,0,getSize().width, getSize().height);
				//Clear the screen (very slow)
				wireframe.paint(backGC);
				g.drawImage(backBuffer, 0, 0, this);
			}
			else
				wireframe.paint(g);
			setPainted();
		}
		else if (message != null) 
		{	g.drawString("Error in model:", 3, 20);
			g.drawString(message, 10, 40);
		}
//SIM		break;
//SIM		}
		paintBoundingBox
		(	g, 
			(float)domain.getX0(0), (float)domain.getX0(1), (float)domain.getX0(2),
			(float)domain.getXmax(0), (float)domain.getXmax(1), (float)domain.getXmax(2)
		);
///		comboBox3DRenderModes.repaint(); 
///		buttonRotX.repaint();
///		buttonRotY.repaint();
			buttonSet.repaint();
			buttonRun.repaint();
			labelEnergy.repaint();
			Nnodes.repaint();///DDD
//SIM			scrollBarA.repaint();
//SIM			scrollBarX.repaint();
//SIM			scrollBarY.repaint();
//SIM			scrollBarZ.repaint();
	}
	private synchronized void setPainted() 
	{	painted = true;
		notifyAll();
	}
	private synchronized void waitPainted()
	{	while (!painted)
		{	try
			{
				wait();
			}
			catch (InterruptedException e) 
			{
			}
		}
		painted = false;
	}
	public String getAppletInfo() 
	{
		return "Title: Mesh \nAuthor: A.Smirnov, B.Sowers, H.Zhang \nDomain 3D Viewer.";
	}
	public String[][] getParameterInfo() 
	{	String[][] info = 
		{	{"model", "path string", "The path to the model to be displayed in .xyz format (see http://chem.leeds.ac.uk/Project/MIME.html).  Default is model.obj."
			},
			{"scale", "float", "Scale factor.  Default is 1 (i.e. no scale)."
			}
		};
		return info;
	}
	public void rotx(float angle)
	{
		tmat.unit();
		tmat.xrot(angle);
		amat.mult(tmat);
		repaint();
	}
	public void roty(float angle)
	{
		tmat.unit();
		tmat.yrot(angle);
		amat.mult(tmat);
		repaint();
	}
	//Show Bounding BOX
	public void paintBoundingBox
	(	Graphics g,
		float xmin, float ymin, float zmin,
		float xmax, float ymax, float zmax
	)
	{
		Quad faces[] = new Quad[6];
		Color color = new Color(0,0,0);
		faces[0] = new Quad
		(	xmin, ymin, zmin,
			xmax, ymin, zmin,
			xmax, ymax, zmin,
			xmin, ymax, zmin,
			color
		);
		
		faces[1] = new Quad
		(	xmin, ymin, zmax,
			xmax, ymin, zmax,
			xmax, ymax, zmax,
			xmin, ymax, zmax,
			color
		);
		faces[2] = new Quad
		(	xmin, ymin, zmin,
			xmax, ymin, zmin,
			xmax, ymin, zmax,
			xmin, ymin, zmax,
			color
		);
		faces[3] = new Quad
		(	xmin, ymax, zmin,
			xmax, ymax, zmin,
			xmax, ymax, zmax,
			xmin, ymax, zmax,
			color
		);
		faces[4] = new Quad
		(	xmin, ymin, zmin,
			xmin, ymax, zmin,
			xmin, ymax, zmax,
			xmin, ymin, zmax,
			color
		);
		faces[5] = new Quad
		(	xmax, ymin, zmin,
			xmax, ymax, zmin,
			xmax, ymax, zmax,
			xmax, ymin, zmax,
			color
		);
		faces[0].mat.unit();
		faces[0].mat.translate(-(xmin+xmax)/2,
				   -(ymin+ymax)/2,
				   -(zmin+zmax)/2);
		faces[0].mat.mult(amat);
		faces[0].mat.scale(xfac, -xfac, 16*xfac/getSize().width);
		faces[0].mat.scale(scale,scale,scale);
		faces[0].mat.translate(getSize().width/2, getSize().height/2, 8);
		faces[0].mat.translate(deltaX, deltaY, 0.0f);
		faces[1].mat 
		= faces[2].mat 
		= faces[3].mat 
		= faces[4].mat 
		= faces[5].mat 
		= faces[0].mat;
		
		//Setup the current plane indicator
///		Color sliceColor = new Color(255, 0, 0, 64);
///		Quad slice = new Quad
///		(	0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 
///			sliceColor
///		);
///		float position;
///		switch(DrawArea.iplane)
///		{
///			case 0:
///			position = xmin + ( (xmax - xmin)/(xmin+xmax) * DrawArea.iplanepos[0] );
///			slice = new Quad(position, ymin, zmin,
///				 position, ymax, zmin,
///				 position, ymax, zmax,
///				 position, ymin, zmax,
///				 sliceColor);
///			break;
///			case 1:
///			position = ymin + ( (ymax - ymin)/(ymin+ymax) * DrawArea.iplanepos[1] );
///			slice = new Quad(xmin, position, zmin,
///					 xmax, position, zmin,
///					 xmax, position, zmax,
///					 xmin, position, zmax,
///					 sliceColor);
///			break;
///			case 2:
///			position = zmin + ( (zmax - zmin)/(zmin+zmax) * DrawArea.iplanepos[2] );
///			slice = new Quad(xmin, ymin, position,
///					 xmax, ymin, position,
///					 xmax, ymax, position,
///					 xmin, ymax, position,
///					 sliceColor);
///			break;
///		}
///		slice.mat = faces[0].mat;
///		if (showBoundary)
///		{
		   for (int i = 1; i < 6; i++)
			{
				faces[i].setFilled(false);
				faces[i].transformed = false;
				faces[i].paint(g);
			}
///		}
///		if (showPlanePositionIndicator)
///		{
///			slice.setFilled(true);
///			slice.transformed = false;
///			slice.paint(g);
///		}
	}
///-	public void itemStateChanged(ItemEvent e)
///-	{
///-		Object src = e.getSource();
///-		boolean on = e.getStateChange() == ItemEvent.SELECTED;
///-		if (src == checkboxX) XON = on;
///-		else if (src == checkboxY) YON = on;
///-		else if (src == checkboxZ) ZON = on;
///-		init();
///-		try
///-		{	double 
///-				sigma=nodes.potential.lengthscale(),
///-				eta=nodes.potential.strength();
///-			nodes = new Nodes(sigma,eta);
///-			nodes.init(space,domain,nodes.nnodes);
///-			nodes.findBB();
///-			float xw = nodes.xmax - nodes.xmin;
///-			float yw = nodes.ymax - nodes.ymin;
///-			float zw = nodes.zmax - nodes.zmin;
///-			if (yw > xw)
///-			xw = yw;
///-			if (zw > xw)
///-			xw = zw;
///-			float f1 = getSize().width / xw;
///-			float f2 = getSize().height / xw;
///-			xfac = 0.7f * (f1 < f2 ? f1 : f2) * scalefudge;
///-		} catch(Exception ex) 
///-		{	ex.printStackTrace();
///-			nodes = null;
///-			message = ex.toString();
///-			System.out.println(message);
///-		}
///-	}
	class Run implements ActionListener
	{	public void actionPerformed(ActionEvent evt)
		{
			if (evt.getSource() == buttonRun) 
			{
//SIM				wireframe.setinj();
				if(nodes.getnvert()<=nodes.gettrapped())
					buttonRun.setText("End   ");
				else
				if(doRun)
				{	doRun=false;
					buttonRun.setText("Run   ");
					repaint();
				}
				else
				{
					doRun=true;	
					buttonRun.setText("Pause");
					start();
				}
			}
		}
	}
	class Set implements ActionListener
	{	public void actionPerformed(ActionEvent evt)
		{	if (evt.getSource() == buttonSet) 
			{	int nnodes=nodes.getnvert();
				try
				{	String snodes=Nnodes.getText().trim();
					///nnodes=Integer.parseInt(Nnodes.getText().trim());
					if(!snodes.equals(""))
						nnodes=Integer.parseInt(snodes);
				}
				catch(Exception e) 
				{	nnodes=nodes.getnvert();
				}
				if(nnodes>0)
				{	if(nodes.getnvert()!=nnodes)
					{	try
						{	double 
								sigma=nodes.potential.lengthscale(),
								eta=nodes.potential.strength();
							nodes = new Nodes(sigma,eta);
							nodes.init(space,domain,nnodes);
							nodes.findBB();
							float xw = nodes.xmax - nodes.xmin;
							float yw = nodes.ymax - nodes.ymin;
							float zw = nodes.zmax - nodes.zmin;
							if (yw > xw)
							xw = yw;
							if (zw > xw)
							xw = zw;
							float f1 = getSize().width / xw;
							float f2 = getSize().height / xw;
							xfac = 0.7f * (f1 < f2 ? f1 : f2) * scalefudge;
						} catch(Exception e) 
						{	e.printStackTrace();
							nodes = null;
							message = e.toString();
							System.out.println(message);
						}
					}	
				}
				nodes.set();
				doRun=false;
				buttonRun.setText("Run   ");
///			start();
				labelEnergy.setText("Energy: "+nodes.getEnergy());
				///Nnodes.setText(""+nodes.getnvert());///DDD
				repaint();
			}
		}
	}
//SIM	class ChangeVair implements AdjustmentListener
//SIM	{	public void adjustmentValueChanged(AdjustmentEvent e)
//SIM		{
//SIM			if (WireFrame.getVair()!=null)
//SIM			{
//SIM				WireFrame.setVair(e.getValue());
//SIM				wireframe.setinj();
//SIM				repaint();
//SIM			}
//SIM		}
//SIM	}
//SIM	class ChangeXinj implements AdjustmentListener
//SIM	{	public void adjustmentValueChanged(AdjustmentEvent e)
//SIM		{	if (WireFrame.getXinj()!=null)
//SIM			{	WireFrame.setXinj
//SIM				(	(double)e.getValue()/(double)domain.Dim[0],
//SIM					WireFrame.getinj(1),
//SIM					WireFrame.getinj(2)
//SIM				);
//SIM				wireframe.setinj();
//SIM				repaint();
//SIM			}
//SIM		}
//SIM	}
//SIM	class ChangeYinj implements AdjustmentListener
//SIM	{	public void adjustmentValueChanged(AdjustmentEvent e)
//SIM		{	if (WireFrame.getXinj()!=null)
//SIM			{	WireFrame.setXinj
//SIM				(	WireFrame.getinj(0),
//SIM					(double)e.getValue()/(double)domain.Dim[1],
//SIM					WireFrame.getinj(2)
//SIM				);
//SIM				wireframe.setinj();
//SIM				repaint();
//SIM			}
//SIM		}
//SIM	}
//SIM	class ChangeZinj implements AdjustmentListener
//SIM	{	public void adjustmentValueChanged(AdjustmentEvent e)
//SIM		{	if (WireFrame.getXinj()!=null)
//SIM			{	WireFrame.setXinj
//SIM				(	WireFrame.getinj(0),
//SIM					WireFrame.getinj(1),
//SIM					(double)e.getValue()/(double)domain.Dim[2]
//SIM				);
//SIM				wireframe.setinj();
//SIM				repaint();
//SIM			}
//SIM		}
//SIM	}

///    public class Rotate implements ActionListener
///    {
///	public void actionPerformed(ActionEvent evt)
///	{
///	    if (evt.getSource() == buttonRotX)
///		rotx(10);
///	    if (evt.getSource() == buttonRotY)
///		roty(10);
///	}
///    }
}   // end class Mesh


///	public class ChangeRenderMode implements ActionListener
///	{
///		public void actionPerformed(ActionEvent evt)
///		{
///			mode = comboBox3DRenderModes.getSelectedIndex();
///			run();
///		}
///	}
class Atom 
{
///	private static JPanel panel;
	private static Applet applet;
	private static byte[] data;
	private final static int R = 40;
	private final static int hx = 15;
	private final static int hy = 15;
	private final static int bgGrey = 192;
	private final static int nBalls = 16;
	private static int maxr;
	
	private int Rl;
	private int Gl;
	private int Bl;
	private Image balls[];
	
	static 
	{	data = new byte[R * 2 * R * 2];
		int mr = 0;
		for (int Y = 2 * R; --Y >= 0;) 
		{	int x0 = (int) (Math.sqrt(R * R - (Y - R) * (Y - R)) + 0.5);
			int p = Y * (R * 2) + R - x0;
			for (int X = -x0; X < x0; X++) 
			{	int x = X + hx;
				int y = Y - R + hy;
				int r = (int) (Math.sqrt(x * x + y * y) + 0.5);
				if (r > mr)
				mr = r;
				data[p++] = r <= 0 ? 1 : (byte) r;
			}
		}
		maxr = mr;
	}
///	static void setPanel(JPanel p) 
///	{
///		panel = p;
///	}
	static void setApplet(Applet app) 
	{
		applet = app;
	}
	Atom(int Rl, int Gl, int Bl) 
	{
		this.Rl = Rl;
		this.Gl = Gl;
		this.Bl = Bl;
	}
	private final int blend(int fg, int bg, float fgfactor) 
	{
		return (int) (bg + (fg - bg) * fgfactor);
	}
	private void Setup() 
	{
		balls = new Image[nBalls];
		byte red[] = new byte[256];
		red[0] = (byte) bgGrey;
		byte green[] = new byte[256];
		green[0] = (byte) bgGrey;
		byte blue[] = new byte[256];
		blue[0] = (byte) bgGrey;
		for (int r = 0; r < nBalls; r++) 
		{	float b = (float) (r+1) / nBalls;
			for (int i = maxr; i >= 1; --i) 
			{	float d = (float) i / maxr;
				red[i] = (byte) blend(blend(Rl, 255, d), bgGrey, b);
				green[i] = (byte) blend(blend(Gl, 255, d), bgGrey, b);
				blue[i] = (byte) blend(blend(Bl, 255, d), bgGrey, b);
			}
			IndexColorModel model = new IndexColorModel
			(	8, maxr + 1,
				red, green, blue, 0
			);
			balls[r] = applet.createImage
			(	new MemoryImageSource
				(R*2, R*2, model, data, 0, R*2)
			);
		}
	}
	void paint(Graphics gc, int x, int y, int r) 
	{
		Image ba[] = balls;
		if (ba == null) 
		{
			Setup();
			ba = balls;
		}
		Image i = ba[r];
		int size = 10 + r;
		gc.drawImage(i, x - (size >> 1), y - (size >> 1), size, size, applet);
	}
}
//Representation of a quadrilateral to be used for surface rendering
class Quad
{
	float vert[] = new float[12];
	int tvert[] = new int[12];
	Color color;
	int nvert = 0;
	boolean transformed;
	Matrix3D mat = new Matrix3D();
	boolean isFilled = true;
	
	//Create a quad.  Note that the vertices must be organized in counter-clockwise fashion
	Quad
	(	float v1x, float v1y, float v1z,
		float v2x, float v2y, float v2z,
		float v3x, float v3y, float v3z,
		float v4x, float v4y, float v4z,
		Color color
	)
	{
		addVert(v1x, v1y, v1z);
		addVert(v2x, v2y, v2z);
		addVert(v3x, v3y, v3z);
		addVert(v4x, v4y, v4z);
		this.color = color;
	}
	void setFilled(boolean isFilled)
	{
		this.isFilled = isFilled;
	}
	void transform() 
	{
		if (transformed || nvert <= 0)return;
		if (tvert == null || tvert.length < nvert * 3)
			tvert = new int[nvert*3];
		mat.transform(vert, tvert, nvert);
		transformed = true;
	}
	public void paint(Graphics g)
	{
		if (vert == null || nvert <= 0)return;
		transform();
			int lg = -1;
		int v[] = tvert;
		if (nvert <= 0)
			return;
	//for (int i = 0; i < lim; i++) 
	//{	
	/*
	  int T = c[i];
	  int 
	  p1=((T >> 16) & 0xFFFF) * 3,
	  p2=(T & 0xFFFF) * 3,
	  grey = v[p1 + 2] + v[p2 + 2];
	  if (grey < 0)grey = 0;
	  if (grey >= nshades)	grey = nshades-1;
	  if (i*nshades+grey != lg) 
	  {	lg = i*nshades+grey;
	  g.setColor(col[grey][(int)com[i]]);
	  }
	*/
		g.setColor(color);
		int xPoints[] = { v[0], v[3], v[6], v[9] };
		int yPoints[] = { v[1], v[4], v[7], v[10] };
		int iShades[] = { v[2], v[5], v[8], v[11] };
		
		if (isFilled) 
			g.fillPolygon(xPoints, yPoints, 4);
		else
		{	//g.drawPolygon(xPoints, yPoints, 4);
			Shades.drawShadedLine
			(	g, v[0], v[1], v[3], v[4],  
				Shades.getShadeColor(iShades[0], color), 
				Shades.getShadeColor(iShades[1], color)
			);
			Shades.drawShadedLine
			(	g, v[3], v[4], v[6], v[7],  
				Shades.getShadeColor(iShades[1], color), 
				Shades.getShadeColor(iShades[2], color)
			);
			Shades.drawShadedLine
			(	g, v[6], v[7], v[9], v[10], 
				Shades.getShadeColor(iShades[2], color), 
				Shades.getShadeColor(iShades[3], color)
			);
			Shades.drawShadedLine
			(	g, v[9], v[10], v[0], v[1], 
				Shades.getShadeColor(iShades[3], color), 
				Shades.getShadeColor(iShades[0], color)
			);
		}
	}
	int addVert(float x, float y, float z) 
	{
		int i = nvert;
		i *= 3;
		vert[i] = x;
		vert[i + 1] = y;
		vert[i + 2] = z;
		return nvert++;
	}
}

