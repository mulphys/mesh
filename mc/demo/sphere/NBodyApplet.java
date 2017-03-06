/**
 * Hanzhou Zhang
 */

import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;
import com.sun.j3d.utils.applet.MainFrame;
import com.sun.j3d.utils.geometry.*;
import com.sun.j3d.utils.universe.*;
import javax.media.j3d.*;
import javax.vecmath.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;
import com.sun.j3d.utils.behaviors.vp.*;

public class NBodyApplet extends JApplet 
implements ActionListener 
{
	BranchGroup scene;
	Transform3D[] v_t3d;
	Sphere[] particles;
	TransformGroup[] v_trans;
	int count; // the number of particles
	int pairs; // pairs of particles
	int iterations;  //the number of iterations for monte carlo simulations
	double[][] x;
    JButton start; 
	JTextField n_ptcl, n_iter;
	JLabel message;
	static float R = 1.0f;
	static double scale=0.5/R;
	MonteCarloSphere mcs;
    Material mat1, altMat;			   
    Appearance app, otherApp;			   
    BoundingSphere worldBounds;
    AlternateAppearance altApp;

    // Globally used colors
	Color3f white = new Color3f(1.0f, 1.0f, 1.0f);
	Color3f red = new Color3f(1.0f, 0.0f, 0.0f);
	Color3f green = new Color3f(0.0f, 1.0f, 0.0f);
	Color3f blue = new Color3f(0.0f, 0.0f, 1.0f);
	Color3f[] colors = {white, red, green, blue};

    private SimpleUniverse u;
    
    public NBodyApplet() 
	{
    }

    public void init() 
	{
		Container contentPane = getContentPane();
	
        Canvas3D c = new Canvas3D(SimpleUniverse.getPreferredConfiguration());
        contentPane.add("Center", c);
 
		// SimpleUniverse is a Convenience Utility class
		u = new SimpleUniverse(c);

		// add mouse behaviors to the viewingPlatform
		ViewingPlatform viewingPlatform = u.getViewingPlatform();

        // This will move the ViewPlatform back a bit so the
        // objects in the scene can be viewed.
        viewingPlatform.setNominalViewingTransform();

		OrbitBehavior orbit = new OrbitBehavior(c,OrbitBehavior.REVERSE_ALL);
		BoundingSphere bounds = new BoundingSphere(new Point3d(0.0, 0.0, 0.0),
												   100.0);
		orbit.setSchedulingBounds(bounds);
		viewingPlatform.setViewPlatformBehavior(orbit);

		// Create GUI
		JPanel p0 = new JPanel();
		JPanel p = new JPanel();
		BoxLayout boxlayout = new BoxLayout(p0, BoxLayout.Y_AXIS);
		JLabel n_particles = new JLabel("# of particles");
		JLabel n_iterations = new JLabel("# of iterations");
		n_ptcl = new JTextField("100", 10);
		n_iter = new JTextField("200", 10);
		start = new JButton("Run");
		message = new JLabel();
		n_ptcl.addActionListener(this);
		n_iter.addActionListener(this);
		n_ptcl.addKeyListener(new KeyAdapter() {
				public void keyTyped(KeyEvent e) {
					char c = e.getKeyChar();      
					if (!((Character.isDigit(c)) ||
						  (c == KeyEvent.VK_BACK_SPACE) ||
						  (c == KeyEvent.VK_DELETE))) {
						getToolkit().beep();
						e.consume();
					}
						}
				});
		n_iter.addKeyListener(new KeyAdapter() {
				public void keyTyped(KeyEvent e) {
					char c = e.getKeyChar();      
					if (!((Character.isDigit(c)) ||
						  (c == KeyEvent.VK_BACK_SPACE) ||
						  (c == KeyEvent.VK_DELETE))) {
						getToolkit().beep();
						e.consume();
					}
						}
				});

		start.addActionListener(this);
		p.add(n_particles);
		p.add(n_ptcl);
		p.add(n_iterations);
		p.add(n_iter);
		p.add(start);
		p0.add(p);
		p0.add(message);

		p0.setLayout(boxlayout);

	
		contentPane.add("South", p0);

        //BranchGroup scene = createSceneGraph();
		BranchGroup sphere = createSceneGraph();
		scene = createParticles();
		u.addBranchGraph(sphere);
        u.addBranchGraph(scene);

    }

    public void destroy() 
	{
		u.cleanup();
    }
    
    BranchGroup createSceneGraph() 
	{
		BranchGroup objRoot = new BranchGroup();

		// Create influencing bounds
		worldBounds = new BoundingSphere(new Point3d( 0.0, 0.0, 0.0 ), //Center
										 1000.0 );     // Extent

		Transform3D t = new Transform3D();
		// move the object upwards
		//t.set(new Vector3f(0.0f, 0.1f, 0.0f));
		// Shrink the object 
		t.setScale(scale);

		TransformGroup trans = new TransformGroup(t);
		trans.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		trans.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
		
		otherApp = new Appearance();
		altMat = new Material();
		altMat.setCapability(Material.ALLOW_COMPONENT_WRITE);
		altMat.setDiffuseColor( new Color3f( 0.0f, 1.0f, 0.0f ) );
		otherApp.setMaterial(altMat);

		altApp = new AlternateAppearance();
		altApp.setAppearance(otherApp);
		altApp.setCapability(AlternateAppearance.ALLOW_SCOPE_WRITE);
		altApp.setCapability(AlternateAppearance.ALLOW_SCOPE_READ);
		altApp.setInfluencingBounds( worldBounds );
		objRoot.addChild(altApp);
	
		Appearance app1 = new Appearance();
		mat1 = new Material();
		mat1.setCapability(Material.ALLOW_COMPONENT_WRITE);
		mat1.setDiffuseColor( new Color3f( 1.0f, 0.0f, 0.0f ) );
		app1.setMaterial(mat1);

		//float R = 1.0f;
		Sphere sphere = new Sphere(R, 
								   Primitive.GENERATE_NORMALS,//generate normals
								   64,         // 64 divisions radially
								   app1);
		trans.addChild(sphere);

		// Add lights
		DirectionalLight light1 = null;
		light1 = new DirectionalLight( );
		light1.setEnable( true );
		light1.setColor( new Color3f(0.2f, 0.2f, 0.2f) );
		light1.setDirection( new Vector3f( 1.0f, 0.0f, -1.0f ) );
		light1.setInfluencingBounds( worldBounds );
		objRoot.addChild( light1 );

		DirectionalLight light2 = new DirectionalLight();
		light2.setEnable(true);
		light2.setColor(new Color3f(0.2f, 0.2f, 0.2f));
		light2.setDirection(new Vector3f(-1.0f, 0.0f, 1.0f));
		light2.setInfluencingBounds(worldBounds);
		objRoot.addChild(light2);
		
		// Add an ambient light to dimly illuminate the rest of
		// the shapes in the scene to help illustrate that the
		// directional lights are being scoped... otherwise it looks
		// like we're just removing shapes from the scene
		AmbientLight ambient = new AmbientLight( );
		ambient.setEnable( true );
		ambient.setColor( new Color3f(1.0f, 1.0f, 1.0f) );
		ambient.setInfluencingBounds( worldBounds );
		objRoot.addChild( ambient );

		objRoot.addChild(trans);

		return objRoot;
    }

    BranchGroup createParticles() 
	{
		BranchGroup scene = new BranchGroup();
		scene.setCapability(BranchGroup.ALLOW_CHILDREN_EXTEND);
		scene.setCapability(BranchGroup.ALLOW_CHILDREN_READ);
		scene.setCapability(BranchGroup.ALLOW_CHILDREN_WRITE);
		scene.setCapability(BranchGroup.ALLOW_DETACH);

		// Create influencing bounds
		worldBounds = new BoundingSphere(new Point3d( 0.0, 0.0, 0.0 ), //Center
										 1000.0 );     // Extent

		Transform3D t = new Transform3D();
		// move the object upwards
		//t.set(new Vector3f(0.0f, 0.1f, 0.0f));
		// Shrink the object 
		t.setScale(scale);

		TransformGroup trans = new TransformGroup(t);
		trans.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
		trans.setCapability(TransformGroup.ALLOW_TRANSFORM_READ);

		otherApp = new Appearance();
		altMat = new Material();
		altMat.setCapability(Material.ALLOW_COMPONENT_WRITE);
		altMat.setDiffuseColor( new Color3f( 0.0f, 1.0f, 0.0f ) );
		otherApp.setMaterial(altMat);
		/*
		altApp = new AlternateAppearance();
		altApp.setAppearance(otherApp);
		altApp.setCapability(AlternateAppearance.ALLOW_SCOPE_WRITE);
		altApp.setCapability(AlternateAppearance.ALLOW_SCOPE_READ);
		altApp.setInfluencingBounds( worldBounds );
		scene.addChild(altApp);
	
		Appearance app1 = new Appearance();
		mat1 = new Material();
		mat1.setCapability(Material.ALLOW_COMPONENT_WRITE);
		mat1.setDiffuseColor( new Color3f( 1.0f, 0.0f, 0.0f ) );
		app1.setMaterial(mat1);

		float R = 1.0f;
		Sphere sphere = new Sphere(R, 
								   Primitive.GENERATE_NORMALS,//generate normals
								   64,         // 64 divisions radially
								   app1);
		trans.addChild(sphere);
		*/
		
		//randomly generate coordinates for small spheres

		String s = n_ptcl.getText();
		if(s.equals(""))
		{
			n_ptcl.setText("100");
			count = 100;
		}
		else
			count = Integer.parseInt(s);

		//System.out.println(s+" "+count);
		pairs = count*(count-1)/2;

		x = new double[count][3];
		particles = new Sphere[count];
		v_t3d = new Transform3D[count];
		v_trans = new TransformGroup[count];
		Vector3d vec = new Vector3d( );

		double theta, phi;
		float r = (float)(0.01/scale);
		for ( int i = 0; i < count; i++ )
		{
			theta = Math.random()*Math.PI;
			phi = Math.random()*Math.PI*2.0;
			x[i][0] = R*Math.sin(theta)*Math.cos(phi);
			x[i][1] = R*Math.sin(theta)*Math.sin(phi);
			x[i][2] = R*Math.cos(theta);

			vec.set( x[i][0]*scale,x[i][1]*scale,x[i][2]*scale);
			v_t3d[i] = new Transform3D();
			v_t3d[i].setScale(scale);
			v_t3d[i].setTranslation( vec );
			v_trans[i] = new TransformGroup(v_t3d[i]);
			v_trans[i].setCapability(TransformGroup.ALLOW_TRANSFORM_READ);
			v_trans[i].setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
			particles[i] = new Sphere(r,     // sphere radius
									  Primitive.GENERATE_NORMALS,// generate normals
									  16,
									  otherApp );  // it's appearance
			v_trans[i].addChild(particles[i]);
			scene.addChild(v_trans[i]);
		}

		double[] center = {0,0,0};
		mcs = new MonteCarloSphere(count, R, center, x, this);
		double potential = mcs.getPotential();
		double ppp = potential/pairs;
		message.setText("System potential: " + potential +
						"          Potential/pair: " + ppp);

		// Add lights
		DirectionalLight light1 = null;
		light1 = new DirectionalLight( );
		light1.setEnable( true );
		light1.setColor( new Color3f(0.2f, 0.2f, 0.2f) );
		light1.setDirection( new Vector3f( 1.0f, 0.0f, -1.0f ) );
		light1.setInfluencingBounds( worldBounds );
		scene.addChild( light1 );

		DirectionalLight light2 = new DirectionalLight();
		light2.setEnable(true);
		light2.setColor(new Color3f(0.2f, 0.2f, 0.2f));
		light2.setDirection(new Vector3f(-1.0f, 0.0f, 1.0f));
		light2.setInfluencingBounds(worldBounds);
		scene.addChild(light2);
		
		// Add an ambient light to dimly illuminate the rest of
		// the shapes in the scene to help illustrate that the
		// directional lights are being scoped... otherwise it looks
		// like we're just removing shapes from the scene
		AmbientLight ambient = new AmbientLight( );
		ambient.setEnable( true );
		ambient.setColor( new Color3f(1.0f, 1.0f, 1.0f) );
		ambient.setInfluencingBounds( worldBounds );
		scene.addChild( ambient );

		scene.addChild(trans);

		return scene;
    }

    public void actionPerformed(ActionEvent e) {

		Object target = e.getSource();
		if(target==start)
		{
			MCthread mct = new MCthread();
			mct.start();
		}

		if(target==n_ptcl)
		{
			scene.detach();
			scene = createParticles();
			u.addBranchGraph(scene);
		}

    }
			   
	private class MCthread extends Thread
	{
		public MCthread()
		{
		}
		
		public void run()
		{
			//System.out.println(mcs);
			String s = n_iter.getText();
			if(s.equals(""))
			{
				n_iter.setText("200");
				iterations = 200;
			}
			else
				iterations = Integer.parseInt(s);

			n_iter.setEnabled(false);
			n_ptcl.setEnabled(false);
			for(int i=0; i<iterations; i++){
				mcs.iterate(1);

				double potential = mcs.getPotential();
				//potential per pair
				double ppp = potential/pairs;
				message.setText("System potential: " + potential +
							"          Potential/pair: " + ppp);
			}
			n_ptcl.setEnabled(true);
			n_iter.setEnabled(true);
			writeParticles();
		}
	}
			   
    public static void main(String[] args) {
	Frame frame = new MainFrame(new NBodyApplet(), 800, 800);
    }

	public void moveParticles()
	{
		x = mcs.getParticles();
		Vector3d vec = new Vector3d( );
		for(int i=0; i<count; i++)
		{
			vec.set( x[i][0]*scale,x[i][1]*scale,x[i][2]*scale);
			v_t3d[i].setTranslation( vec );
			v_trans[i].setTransform(v_t3d[i]);
		}

	}

	public void writeParticles()
	{
		x = mcs.getParticles();
		for(int i=0; i<count; i++)
		{
			System.out.println(i+"\t"+x[i][0]+"\t"+x[i][1]+"\t"+x[i][2]+"\t");
		}
	}
}			   
