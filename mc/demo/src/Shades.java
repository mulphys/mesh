
import javax.swing.*;
import java.awt.*;

public class Shades
{
	public static int nshades = 16;
	public static int getShades()
	{	return nshades;
	}
	static Color getShadeColor(int ishade, Color color)
	{
		if (ishade < 0)ishade = 0;
		if (ishade >= nshades)ishade = nshades-1;
		double shade=1.0*(1-Math.pow(ishade/((float)nshades-1.0),2.3));
		int
			grey=170,
			r=(int)((1-shade)*color.getRed()+shade*grey),
			g=(int)((1-shade)*color.getGreen()+shade*grey),
			b=(int)((1-shade)*color.getBlue()+shade*grey);
		return new Color(r,g,b);
	}
	static public void drawShadedLine
	(	Graphics g, 
		int x1, int y1, int x2, int y2, 
		Color colorBegin, Color colorEnd
	)
	{	int deltaX = ((x2 - x1) >= 0) ? (x2 - x1) : -(x2 - x1);
		int deltaY = ((y2 - y1) >= 0) ? (y2 - y1) : -(y2 - y1);
		int x = x1, y = y1;
		int xinc1, xinc2, yinc1, yinc2;
		int num, den, numAdd, numPixels;
		//Determine if the line is going left/right
		if (x2 >= x1)
			xinc1 = xinc2 = 1;
		else
			xinc1 = xinc2 = -1;
		//Determine if the line is going up/down
		if (y2 >= y1)
			yinc1 = yinc2 = 1;
		else
			yinc1 = yinc2 = -1;
		//If the slope is less than one, we increment the x values by one 
		// and increment the y value by the slope each iteration
		if (deltaX >= deltaY)
		{
			xinc1 = yinc2 = 0;
			den = deltaX;
			num = deltaX / 2;
			numAdd = deltaY;
			numPixels = deltaX;
		}
		//If the slope is greater than one, we increment the y values by one 
		// and increment the x values by the slope each iteration
		else
		{
			xinc2 = yinc1 = 0;
			den = deltaY;
			num = deltaY / 2;
			numAdd = deltaX;
			numPixels = deltaY;
		}
		//Use floats to make sure the color incremenents properly. 
		// However, we will have to typecast them to ints when setting the actual
		//color
		float rCurr = (float)colorBegin.getRed();
		float gCurr = (float)colorBegin.getGreen();
		float bCurr = (float)colorBegin.getBlue();
		float rInc = (float)(colorEnd.getRed() - colorBegin.getRed())/(float)numPixels;
		float gInc = (float)(colorEnd.getGreen() - colorBegin.getGreen())/(float)numPixels;
		float bInc = (float)(colorEnd.getBlue() - colorBegin.getBlue())/(float)numPixels;
		for (int currPixel = 0; currPixel < numPixels; currPixel++)
		{	g.setColor(new Color((int)rCurr, (int)gCurr, (int)bCurr));		
			g.drawLine(x, y, x, y);
			num += numAdd;
			if (num >= den)
			{	num -= den;
				x += xinc1;
				y += yinc1;
			}
			x += xinc2;
			y += yinc2;	
			rCurr += rInc; gCurr += gInc; bCurr += bInc;
		}
	}
}
