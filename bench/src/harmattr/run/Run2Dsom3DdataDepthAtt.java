package run;

import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.min;

import java.awt.Color;
import java.awt.image.IndexColorModel;
import java.io.File;

import javax.swing.SwingUtilities;

import color2D.ColorMap2D;

import morletwavelettransform.CenterFreqs;
import morletwavelettransform.MorletTransform;
import somseis.SOM2;
import attributes.SimpleAttributes;
import display.MultiplePlot;
import display.SinglePlot;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.sgl.ImagePanelGroup;
import edu.mines.jtk.sgl.SimpleFrame;
import grabbers.DataGrabber;

public class Run2Dsom3DdataDepthAtt {
	public static void main(String[] args) {

		File file = new File("C:/Users/Chris/401_4_600TPData/tpsz.binary");
		String sFile = file.getPath();

		float[][][] rawData = DataGrabber.grab3DDataFromFile(sFile, 161, 357,
				401);

		// File Information
		sN3 = 0;
		eN3 = 160;

		sN2 = 0;
		eN2 = 356;

		sN1 = 0;
		eN1 = 400;
		dt = .002;

		n1 = eN1 - sN1 + 1;
		n2 = eN2 - sN2 + 1;
		n3 = eN3 - sN3 + 1;
		totaltraces = n3 * n2;

		Sampling sampN3 = new Sampling(n3, 1.0f, sN3);
		Sampling sampN2 = new Sampling(n2, 1.0f, sN2);

		// Frequencies to analyze and number of filters
		fmin = 8;// Hz
		fmax = 30.0;// Hz
		numfilters = 2;

		// SOM
		int numIter = 100000;
		int numAttributes = numfilters;
		int n2SOM = 4;
		int n1SOM = 4;
		SOM2 som = new SOM2(numIter, n2SOM, n1SOM, numAttributes);

		final float[][][] shortData = DataGrabber.shorten3DData(rawData, sN3,
				eN3, sN2, eN2, sN1, eN1);
		System.out.println("Trimmed Data Built");

		ImagePanelGroup img3 = new ImagePanelGroup(rawData);
		SimpleFrame sf3 = new SimpleFrame();
		sf3.addImagePanels(img3);
		
		Sampling time = new Sampling(n1, dt, 0);
		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
		Sampling freq = MorletTransform.getLogFrequencySampling(dt);

		

		float[][][][] subbands = gf.apply(shortData);
		System.out.println("Morlet Complete");

		float[][][][] amp = SimpleAttributes.findAmplitude(subbands);
		System.out.println("Amplitude Complete");

		
		int sn4 = amp.length;
		int sn3 = amp[0].length;
		int sn2 = amp[0][0].length;
		int sn1 = amp[0][0][0].length;
		float[][][][] somAmp = new float[sn4][sn3][sn1][sn2 + 1];
		for (int i4 = 0; i4 < sn4; ++i4) {
			for (int i3 = 0; i3 < sn3; ++i3) {
				for (int i2 = 0; i2 < sn2; ++i2) {
					for (int i1 = 0; i1 < sn1; ++i1) {
						somAmp[i4][i3][i1][i2] = amp[i4][i3][i2][i1];
						somAmp[i4][i3][i1][sn2] = i1;

					}
				}
			}
		}
		 

		som.train3D(amp);

		for (int i = 0; i < n2SOM; ++i) {
			for (int j = 0; j < n1SOM; ++j) {
				System.out.println("nodes = " + som.node(i, j));
			}
		}

		float[][][] classData = som.classify3DData(amp);
		ColorMap cm = ColorMap2D.getColorMap(n1SOM, n2SOM, 10);// new
														// ColorMap(ColorMap.HUE_BLUE_TO_RED);
		ImagePanelGroup img1 = new ImagePanelGroup(classData);
		img1.setColorModel(cm.getColorModel());

		ImagePanelGroup img2 = new ImagePanelGroup(shortData);
		// img2.setColorModel(cm.getColorModel());
		SimpleFrame sf1 = new SimpleFrame();
		sf1.addImagePanels(img1);
		sf1.addImagePanels(img2);

		PlotPanelPixels3 ppp31 = new PlotPanelPixels3(
				PlotPanelPixels3.Orientation.X1RIGHT_X2UP,
				PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM, time, sampN2,
				sampN3, classData);
		ppp31.setLineColor(Color.GREEN);
		ppp31.addColorBar();
		int slipor1 = (int) ((eN3 - sN3) / 2.0f);
		int sliro1 = (int) ((eN2 - sN2) / 2.0f);
		int slivp1 = (int) ((eN1 - sN1) / 2.0f);
		System.out.println("Por = " + slipor1 + "Ro = " + sliro1 + "Vp = "
				+ slivp1);
		ppp31.setSlices(slivp1, sliro1, slipor1);

		ppp31.setColorModel(cm.getColorModel());
		PlotFrame pf1 = new PlotFrame(ppp31);
		pf1.setVisible(true);

		
	}


	private static int numfilters;
	private static int sN1, eN1, n1, sN2, eN2, n2, sN3, eN3, n3, totaltraces;
	private static double dt, fmin, fmax;
}
