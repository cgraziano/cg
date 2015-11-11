package run;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteOrder;

import javax.swing.SwingUtilities;

import morletwavelettransform.MorletTransform;
import somseis.SOM2;
import attributes.SimpleAttributes;
import color2D.ColorMap2D;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.sgl.ImagePanelGroup;
import edu.mines.jtk.sgl.SimpleFrame;
import grabbers.DataGrabber;

public class WriteClassData {
	public static void main(String[] args) {

		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");

		String sFile = file.getPath();

		float[][][] rawData = DataGrabber.grab3DDataFromFile(sFile, 161, 357,
				401);

		File fileR = new File(
				"C:/Users/Chris/Data/ClassData/Norm8to60Hz3F4x4som.dat");

		// File Information
		sN3 = 0;
		eN3 = 160;

		sN2 = 0;
		eN2 = 356;

		sN1 = 0;
		eN1 = 400;
		dt = .004;

		n1 = eN1 - sN1 + 1;
		n2 = eN2 - sN2 + 1;
		n3 = eN3 - sN3 + 1;

		Sampling sampN3 = new Sampling(n3, 1.0f, sN3);
		Sampling sampN2 = new Sampling(n2, 1.0f, sN2);

		// Frequencies to analyze and number of filters
		fmin = 8;// Hz
		fmax = 60.0;// Hz
		numfilters = 3;

		// SOM
		int numIter = 100000;
		int numAttributes = numfilters;
		int n2SOM = 4;
		int n1SOM = 4;
		int nClass = n2SOM * n1SOM;
		SOM2 som = new SOM2(numIter, n2SOM, n1SOM, numAttributes);

		// final float[][][] shortData = DataGrabber.shorten3DData(rawData, sN3,
		// eN3, sN2, eN2, sN1, eN1);
		System.out.println("Trimmed Data Built");

		Sampling time = new Sampling(n1, dt, 0);
		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
		Sampling freq = MorletTransform.getLogFrequencySampling(dt);

		float[][][][] subbands = gf.apply(rawData);
		System.out.println("Morlet Complete");

		float[][][][] amp = SimpleAttributes.findAmplitude(subbands);
		System.out.println("Amplitude Complete");

		float[][][][] somAmp = new float[n3][n2][time.getCount()][freq
				.getCount()];

		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int i1 = 0; i1 < time.getCount(); ++i1) {
					for (int f = 0; f < freq.getCount(); ++f) {
						somAmp[i3][i2][i1][f] = amp[f][i3][i2][i1];
					}
				}
			}
		}
		
		som.train3D(somAmp);
		float[][][] classData = som.classify3DData(somAmp);

		ArrayFile af;
		try {
			af = new ArrayFile(fileR, "rw", ByteOrder.LITTLE_ENDIAN,
					ByteOrder.BIG_ENDIAN);
			try {
				af.writeFloats(classData);
				System.out.println("ClassData Written");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} catch (FileNotFoundException e) {
			System.out.println("FileNotFound");
		}
	}

	private static int numfilters;
	private static int sN1, eN1, n1, sN2, eN2, n2, sN3, eN3, n3;
	private static double dt, fmin, fmax;

}
