package run;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteOrder;

import morletwavelettransform.MorletTransform;
import somseis.SOM2;
import attributes.SimpleAttributes;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.sgl.ImagePanelGroup;
import edu.mines.jtk.sgl.SimpleFrame;
import grabbers.DataGrabber;

public class Write3Frequencies {
	public static void main(String[] args) {
		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");
		String sFile = file.getPath();

		float[][][] rawData = DataGrabber.grab3DDataFromFile(sFile, 161, 357,
				401);

		// File Information
		int sN3 = 0;
		int eN3 = 160;

		int sN2 = 0;
		int eN2 = 356;

		int sN1 = 0;
		int eN1 = 400;
		double dt = .004;

		int n1 = eN1 - sN1 + 1;
		int n2 = eN2 - sN2 + 1;
		int n3 = eN3 - sN3 + 1;
		int totaltraces = n3 * n2;

		Sampling sampN3 = new Sampling(n3, 1.0f, sN3);
		Sampling sampN2 = new Sampling(n2, 1.0f, sN2);

		// Frequencies to analyze and number of filters
		double fmin = 8;// Hz
		double fmax = 50.0;// Hz
		int numfilters = 3;

		Sampling time = new Sampling(n1, dt, 0);
		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
		Sampling freq = MorletTransform.getLogFrequencySampling(dt);

		float[][][][] subbands = gf.apply(rawData);
		System.out.println("Morlet Complete");

		
		float[][][][] amp = SimpleAttributes.findAmplitude(subbands);
		System.out.println("Amplitude Complete");
		
		float[][][][] ampF = new float[freq.getCount()][n3][n2][time.getCount()]; 
		File mfile = new
		File("C:/Users/Chris/Data/FrequencyData/MorletSamp.dat");
		gf.writeSampling(mfile, gf.getLogFrequencySampling(dt));
		
		ArrayFile af; File fileF0 = new
		File("C:/Users/Chris/Data/FrequencyData/frequency0.dat"); try { af =
		new ArrayFile(fileF0, "rw",ByteOrder.LITTLE_ENDIAN,
		ByteOrder.BIG_ENDIAN);
		
		af.writeFloats(amp[0]); System.out.println("Frequency0 Written");
		
		
		
		
		} catch (IOException e) { System.out.println("FileNotFound"); }
		
		
		File fileF1 = new
		File("C:/Users/Chris/Data/FrequencyData/frequency1.dat"); try { af =
		new ArrayFile(fileF1, "rw",ByteOrder.LITTLE_ENDIAN,
		ByteOrder.BIG_ENDIAN);
		
		af.writeFloats(amp[1]); System.out.println("Frequency1 Written");
		
		
		 
		 
		} catch (IOException e) { System.out.println("FileNotFound"); }
		
		File fileF2 = new
		File("C:/Users/Chris/Data/FrequencyData/frequency2.dat"); try { af =
		new ArrayFile(fileF2, "rw",ByteOrder.LITTLE_ENDIAN,
		ByteOrder.BIG_ENDIAN);
		
		af.writeFloats(amp[2]); System.out.println("Frequency2 Written");
		
		
		} catch (IOException e) { System.out.println("FileNotFound"); }
		
	}
}
