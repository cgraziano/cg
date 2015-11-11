package run;

import java.io.File;

import morletwavelettransform.MorletTransform;
import attributes.SimpleAttributes;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import grabbers.DataGrabber;

public class PlotSweep {
	public static void main(String[] args) {

		File file = new File("sweep5-85Hzdt=.004n=300,.binary");
		String sFile = file.getPath();

		float[][] sweep = DataGrabber.grab2DDataFromFile(sFile, 2142, 300);

		System.out.println(sweep.length);
		System.out.println(sweep[0].length);
		// File Information
		int sN2 = 0;
		int eN2 = 2141;

		int sN1 = 0;
		int eN1 = 299;
		double dt = .004;

		int n1 = eN1 - sN1 + 1;
		int n2 = eN2 - sN2 + 1;

		Sampling sampN2 = new Sampling(n2, 1.0f, sN2);
		Sampling time = new Sampling(n1, dt, 0);

		// Frequencies to analyze and number of filters
		double fmin = 10;// Hz
		double fmax = 50;// Hz
		int numfilters = 10;

		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
		Sampling freq = MorletTransform.getLogFrequencySampling(dt);

		float[][][] subbands = gf.apply(sweep);
		System.out.println("Morlet Complete");

		float[][][] amp = SimpleAttributes.findAmplitude(subbands);

		// Plotting purposes
		float[][][] plotFreq = new float[n2][numfilters][n1];
		for (int f = 0; f < numfilters; ++f) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int tr = 0; tr < n1; ++tr) {
					plotFreq[i2][f][tr] = amp[f][i2][tr];
				}
			}
			System.out.println("Morlet Frequeny " + f + " done");
		}
		System.out.println("Amplitude Complete");
		PlotPanel pp = new PlotPanel(1, 2, PlotPanel.Orientation.X1DOWN_X2RIGHT);
		pp.addPoints(0, 0, time, sweep[0]);

		System.out.println(amp[0].length);
		System.out.println(amp.length);

		PixelsView pv = new PixelsView(time, freq, plotFreq[0]);
		pv.setColorModel(ColorMap.HUE_RED_TO_BLUE);
		pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv.setInterpolation(PixelsView.Interpolation.NEAREST);
		pp.addTiledView(0, 1, pv);
		pp.setVLabel("Time (s)");
		pp.setHLabel(0, "Seismic Amplitude");
		pp.setHLabel(1, "Frequency (Hz)");
		pp.addColorBar("Morlet Amplitude");
		pp.getMosaic().setWidthElastic(0, 60);
		pp.getMosaic().setWidthElastic(1, 50);

		pp.setHInterval(1, .2);

		PlotFrame pf = new PlotFrame(pp);
		pf.setFontSizeForSlide(.9, .9);
		pf.setVisible(true);
		pf.paintToPng(1000, 5, "SingleTraceSweep_10-50HzSweep.png");
	}
}
