package run;

import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.min;
import display.MultiplePlot;
import display.SinglePlot;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PointsView;
import grabbers.DataGrabber;

import java.awt.image.IndexColorModel;
import java.io.File;

import javax.swing.SwingUtilities;

import color2D.ColorMap2D;

import somseis.SOM;
import somseis.SOM2;

import morletwavelettransform.CenterFreqs;
import morletwavelettransform.MorletTransform;
import attributes.SimpleAttributes;

public class Run2Dsom2Ddata {
	public static void main(String[] args) {
		File file = new File("C://Users/Chris/Documents/!gb.seismic.binary");
		
		String sFile = file.getPath();

		System.out.println(sFile);
		
		 float[][] rawData = DataGrabber.grab2DDataFromFile(sFile,2142,1500);
		  //File Information startCDP = 0; endCDP = 2141; startSample = 0;
		  endSample = 1499; dt = .004;
		  
		  //Frequencies to analyze and number of filters 
		  fmin = 8;//Hz 
		  fmax = 30.0;//Hz 
		  numfilters = 2;
		  
		  //SOM 
		  int numIter = 100000; 
		  int numAttributes = numfilters; 
		  int n2 =8; 
		  int n1 = 8; 
		  SOM2 som = new SOM2(numIter, n2, n1, numAttributes);
		  
		  
		  totalsamples = endSample - startSample + 1; totaltraces = endCDP -
		  startCDP + 1;
		  
		  final float[][] shortData = DataGrabber.shorten2DData(rawData,
		  startCDP, endCDP, startSample, endSample);
		  
		  
		  Sampling time = new Sampling(totalsamples, .004, 0); 
		  MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
		  Sampling freq = MorletTransform.getLogFrequencySampling(dt);
		  
		  
		  
		
		  float[][][] subbands = gf.apply(shortData);
		  System.out.println("Morlet Complete");
		  
		  
		  float[][][] amp = SimpleAttributes.findAmplitude(subbands);
		  
		  float[][][] somAmp = new float[totaltraces][time.getCount()][numfilters]; 
		  for (int i=0; i<totaltraces; ++i){ 
			  for (int j=0; j<time.getCount(); ++j){ 
				  for (int k=0; k<numfilters; ++k){ 
					  somAmp[i][j][k] = amp[i][k][j]; 
				  }
			  }
		  }
				
		  som.train2D(amp);
			  
		  float[][] classData = som.classify2DData(amp);
		  cdp = new Sampling(totaltraces, 1.0,startCDP); 
		  time = new Sampling(totalsamples, dt, startSample);
		  
		  MultiplePlot mp1 = new MultiplePlot(time, cdp, shortData,cdp,classData); 
		  mp1.setNearestInterpolation(1);
		  mp1.setTitle("Categorized Data"); 
		  mp1.setYLabel("Time (s)");
		  mp1.setXLabel(0, "CDP Number"); 
		  mp1.setXLabel(1, "CDP Number");
		  
		  ColorMap cm = ColorMap2D.getColorMapSingleTransparent(n1, n2, 10, 11);
		 
		  mp1.setColorModel(1, cm.getColorModel()); 
		  mp1.setColorBar("Class Number");
	}
		  
		 
		  

	private static int numfilters, totalsamples, totaltraces, startCDP, endCDP,
			startSample, endSample;
	private static double dt, fmin, fmax;
	private static Sampling cdp, time, freq;
	private static float[][][] conAmp, frtrtiamp;
}
