
/**
 * Hanzhou Zhang
 *
 * This class implements the Monte Carlo simulation
 * presented in the paper "MONTE-CARLO AND POLYHEDRON-BASED
 * SIMULATIONS I: EXTREMAL STATES OF THE LOGRITHMIC N-BODY
 * PROBLEM ON A SPHERE". 
 * Also see the paper "THE SPHERICAL MODEL OF LOGARITHMIC 
 * POTENTIALS AS EXAMINED BY MONTE CARLO METHODS" for the 
 * detail about the picking of rotation axis and angle.
 * http://www.rpi.edu/~limc
 */

import java.util.*;
import java.lang.Math;

public class MonteCarloSphere
{
	NBodyApplet applet;
	//number of particles on the sphere
	private int n;
	//the first m particles are fixed, don't move them
	private int m;
	//radius of the sphere
	private double radius ;
	//center coordinates of the sphere
	private double[] center;
	//Cartesian coordinates of particles on the sphere
	private double[][] z0;
	//z = z0 - center
	private double[][] z; 

	private double[] axis;
	private double[] buffer;
	private double beta = 10.0; //inverse of temperature
	private double epslon = .02; //how to choose it?

	Random random4r;
	Random random4Axis1;
	Random random4Axis2;
	Random random4Angle;
	Random random4Particle;

	int n_accepted;

	public MonteCarloSphere(int n, int m, double radius,
							double[] center, double[][] z0,
							NBodyApplet applet)
	{
		this.applet = applet;
		this.n = n;
		this.m = m;
		this.radius = radius;
		this.center = center;
		this.z0 = z0;
		axis = new double[3];
		buffer = new double[3];
		z = new double[n][3];
		for(int i=0; i<n; i++)
			for(int j=0; j<3; j++)
				z[i][j]=z0[i][j]-center[j];

		random4r = new Random();
		random4Axis1 = new Random();
		random4Axis2 = new Random();
		random4Angle = new Random();
		random4Particle = new Random();

		/*
		System.out.println("Initial Potential Energy:"+ 
						   potentialEnergy());
		*/
	}


	public MonteCarloSphere(int n, double radius,
							double[] center, double[][] z0,
							NBodyApplet applet)
	{
		//if the number of fixed particle is not set, set it to 0
		this(n, 0, radius, center, z0, applet);
	}

	//scalar product of vectors
	private double dot(double[] a, double[] b)
	{
		return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
	}

	private double dist(double[] a, double[] b)
	{
		double u = a[0]-b[0];
		double v = a[1]-b[1];
		double w = a[2]-b[2];
		return Math.sqrt(u*u+v*v+w*w);
	}

	public double interaction(int j, int k)
	{
		//return -Math.log(Math.abs(1.0-dot(z[j], z[k])));//*10.0;
		return -Math.log(1.0-dot(z[j], z[k])/(radius*radius));
		//double d = dist(z[j], z[k]);
		//return 1.0/(d*d);
	}

	public double potentialEnergy()
	{
		double H = 0.0;

		for (int j=0; j<n; j++)
		{
			for(int k=j+1; k<n; k++)
			{
				H = H + interaction(j, k);
			}
		}
		return H;
	}

	//potential energy from interactions between one particle
	// and the other particles
	public double potentialEnergy(int p)
	{
		double H = 0.0;

		for (int j=0; j<n; j++)
		{
			if(j!=p)
			{
				H = H + interaction(j, p);
			}
		}
		return H;
	}

	// rotate a point around a random selected axis with a random angle
	private void rotate(int p)
	{
		//randomly select a axis:
		//first select a random number from [-1, 1]
		double zi = random4Axis1.nextDouble()*2.0-1.0;
		//then select a angle theta from [0, 2PI]
		double theta = random4Axis2.nextDouble()*2.0*Math.PI;
		//now define the axis as:
		axis[0] = Math.sqrt(1.0-zi*zi)*Math.cos(theta);
		axis[1] = Math.sqrt(1.0-zi*zi)*Math.sin(theta);
		axis[2] = zi;

		//get a random angle ranging from (0, epslon);
		double angle = epslon * random4Angle.nextDouble();
		rotate3dPoint(z[p], axis, angle);

		for(int k=0; k<3; k++)
			z0[p][k]= z[p][k]+center[k];
	}

	public void sweep()
	{
		n_accepted=0;
		for(int i =0; i<n-m; i++)
		{
			//randomly select a particle from [m, n-1]
			//leave the first point z[0] unchanged
			int p = random4Particle.nextInt(n-m)+m;
			//get the potential energy before rotation
			double pe0 = potentialEnergy(p); 
			//save the coordinates of the selected point
			for(int k=0; k<3; k++)
				buffer[k] = z[p][k];

			//rotate the selected point
			rotate(p);

			//apply the changes to zz

			//calculate the change of potential energy:
			double dh = potentialEnergy(p) - pe0;
			//if(dh<0.0) System.out.println("dh="+ dh);

			if(!isAccepted(dh))
			{
				//reverse the change
				for(int k=0; k<3; k++)
					z[p][k]=buffer[k];
			}
		}
		applet.moveParticles();
		double ratio = (double)n_accepted/(n-m);
		//System.out.println("         !!!!!!!!!!!!!!!!!!!!!! "+ratio);

		if(ratio<0.5) beta = 0.95*beta;
		if(ratio>0.6) beta = 1.05*beta;
		//System.out.println("          !!!!!!!!!!!!!!!!!!!!!!!!!!!! "+ratio+
		//			   " "+ beta);
		
	}

	public void iterate(int times)
	{
		for(int i=0; i<times; i++)
			sweep();

		for(int i=0; i<n; i++)
			for(int j=0; j<3; j++)
				z0[i][j]=z[i][j]+center[j];

	}


	public boolean isAccepted(double dh)
	{
		double r = random4r.nextDouble();
		double tmp = Math.exp(-beta*dh);

		//System.out.print(" exp(-beta*dh)="+tmp+" r="+r);
		if(r<tmp)
		{
			n_accepted++;
			//System.out.print("dh="+ dh);
			//System.out.println(" accepted");
			return true;
		}
		else
		{
			//System.out.print("                                          dh="+ dh);
			//System.out.println(" rejected");
			return false;
		}
	}

	public String toString()
	{
	    String s= getClass().getName();

		for(int i=0; i<n; i++)
		{
			/*
			double R = Math.sqrt(z[i][0]*z[i][0] +
								 z[i][1]*z[i][1]+
								 z[i][2]*z[i][2]);
			*/
			s = s + "["+z0[i][0] +
				"," + z0[i][1] +
				"," + z0[i][2] + 
				"] ";
		}
		double pe = potentialEnergy();
		s = s + "Final Potential Energy=" + pe;

		return s;
	}

	//rotate a 3D point p around an axis(origin - axis[]) with a angle

	public  static void rotate3dPoint(double[] p, 
										  double[] axis, 
										  double angle)
	/*

	public  static double[] rotate3dPoint(double[] p, 
										  double[] axis, 
										  double angle)
	*/
	{
		//make the axis to be unit vector
		double len = Math.sqrt(axis[0]*axis[0]+
							   axis[1]*axis[1]+
							   axis[2]*axis[2]);
		axis[0] = axis[0]/len;
		axis[1] = axis[1]/len;
		axis[2] = axis[2]/len;

		//find the quaternion (x, y, z, w):
		double x, y, z, w;
		x = axis[0]*Math.sin(angle/2.0);
		y = axis[1]*Math.sin(angle/2.0);
		z = axis[2]*Math.sin(angle/2.0);
		w = Math.cos(angle/2.0);

		//construct the rotation matrix
		//    M11    M12     M13
		//    M21    M22     M23
		//    M31    M32     M33

		double m11 = 1.0-2.0*(y*y+z*z);
		double m22 = 1.0-2.0*(x*x+z*z);
		double m33 = 1.0-2.0*(x*x+y*y);

		double m12 = 2.0*(x*y-w*z);
		double m21 = 2.0*(x*y+w*z);

		double m13 = 2.0*(x*z+w*y);
		double m31 = 2.0*(x*z-w*y);

		double m23 = 2.0*(y*z-w*x);
		double m32 = 2.0*(y*z+w*x);

		//rotate it:
		double[] p1 = new double[3];
		p1[0] = m11*p[0]+m12*p[1]+m13*p[2];
		p1[1] = m21*p[0]+m22*p[1]+m23*p[2];
		p1[2] = m31*p[0]+m32*p[1]+m33*p[2];

		//return p1;

		for(int i=0; i<3; i++)
			p[i] = p1[i];
	}

	public double[][] getParticles()
	{
		return z;
	}

	public double getPotential()
	{
		return potentialEnergy();
	}
}

