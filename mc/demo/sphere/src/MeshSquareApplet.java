/**
 * Hanzhou Zhang
 */

import java.applet.Applet;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;
import javax.swing.border.*;

public class MeshSquareApplet extends JApplet 
implements ActionListener 
{
    JButton start; 
    public MeshSquareApplet() 
	{
    }

    public void init() 
	{
		Container contentPane = getContentPane();
	
		// Create GUI
		JPanel p0 = new JPanel();
		JPanel p = new JPanel();
		BoxLayout boxlayout = new BoxLayout(p0, BoxLayout.Y_AXIS);
		start = new JButton("Run");
		start.addActionListener(this);
		p.add(start);
		p0.add(p);
		p0.setLayout(boxlayout);

		contentPane.add("South", p0);
    }

	public void paint( Graphics g ) {

//      g.setColor( Color.black);
//      //g.drawRect( 10, 20, 100, 15 );
//      g.fillRect( 0, 0, 800, 400 );
//      g.setColor( Color.white);
//	  g.fillOval(10,10,3,3);
	}

    public void destroy() 
	{

    }
    

    public void actionPerformed(ActionEvent e) {

		Object target = e.getSource();
//		if(target==start)
//		{
//			MCthread mct = new MCthread();
//			mct.start();
//		}
//
//		if(target==n_ptcl)
//		{
//			scene.detach();
//			scene = createParticles();
//			u.addBranchGraph(scene);
//		}

    }
			   
	private class MCthread extends Thread
	{
		public MCthread()
		{
		}
		
		public void run()
		{
			//System.out.println(mcs);
		}
	}
			   
    public static void main(String[] args) 
	{
		Frame frame = new MainFrame(new MeshSquareApplet(), 800, 400);
    }

}			   
