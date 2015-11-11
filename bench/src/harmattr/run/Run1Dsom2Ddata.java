package run;

import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.min;

import java.awt.image.IndexColorModel;
import java.io.File;

import javax.swing.SwingUtilities;

import morletwavelettransform.CenterFreqs;
import morletwavelettransform.MorletTransform;
import somseis.SOM1;
import somseis.SOM2;
import attributes.SimpleAttributes;
import display.MultiplePlot;
import display.SinglePlot;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import grabbers.DataGrabber;

public class Run1Dsom2Ddata {
	public static void main(String[] args) {
		File file = new File("C://Users/Chris/Documents/!gb.seismic.binary");
	
		String sFile = file.getPath();

		float[][] rawData = DataGrabber.grab2DDataFromFile(sFile, 2142, 1500);
		// File Information
		startCDP = 50;
		endCDP = 52;
		startSample = 0;
		endSample = 500;
		dt = .004;

		// Frequencies to analyze and number of filters
		fmin = 8;// Hz
		fmax = 30.0;// Hz
		numfilters = 2;

		// SOM
		int numIter = 100000;
		int numAttributes = numfilters;

		int n1 = 24;
		SOM1 som = new SOM1(numIter, n1, numAttributes);

		totalsamples = endSample - startSample + 1;
		totaltraces = endCDP - startCDP + 1;

		final float[][] shortData = DataGrabber.shorten2DData(rawData,
				startCDP, endCDP, startSample, endSample);

		Sampling time = new Sampling(totalsamples, .004, 0);
		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);
		Sampling freq = MorletTransform.getLogFrequencySampling(dt);


		float[][][] subbands = gf.apply(shortData);
		System.out.println("Morlet Complete");

		float[][][] amp = SimpleAttributes.findAmplitude(subbands);

		frtrtiamp = new float[freq.getCount()][totaltraces][time.getCount()];

		float[][][] somAmp = new float[totaltraces][time.getCount()][numfilters];
		for (int i = 0; i < totaltraces; ++i) {
			for (int j = 0; j < time.getCount(); ++j) {
				for (int k = 0; k < numfilters; ++k) {
					somAmp[i][j][k] = amp[i][k][j];
				}
			}
		}

		som.train2D(somAmp);

		for (int i = 0; i < n1; ++i) {
			System.out.println("nodes = " + som.node(i));
		}
		float[][] classData = som.classify2DData(somAmp);
		cdp = new Sampling(totaltraces, 1.0, startCDP);
		time = new Sampling(totalsamples, dt, startSample);
		MultiplePlot mp1 = new MultiplePlot(time, cdp, shortData, cdp,
				classData);
		mp1.setNearestInterpolation(1);
		mp1.setTitle("Categorized Data");
		mp1.setYLabel("Time (s)");
		mp1.setXLabel(0, "CDP Number");
		mp1.setXLabel(1, "CDP Number");

		ColorMap cm = som.getColorMap(0, n1);// new
												// ColorMap(ColorMap.HUE_BLUE_TO_RED);
		// mp1.setXLabel(3, "CDP Number for "+freq[2]+" Hz");
		mp1.setFont(14);
		mp1.setColorModel(1, cm.getColorModel());
		// mp1.setColorModel(2, cm.getColorModel());
		// mp1.setColorModel(3, cm.getColorModel());
		mp1.setColorBar("Category Number");

		MultiplePlot mpT = new MultiplePlot(time, cdp, shortData, cdp,
				classData, cm, .5f);
		mpT.setNearestInterpolation(1);
		mpT.setTitle("Categorized Data");
		mpT.setYLabel("Time (s)");
		mpT.setXLabel(0, "CDP Number");

		// mp1.setXLabel(3, "CDP Number for "+freq[2]+" Hz");
		mpT.setFont(14);
		// mpT.setColorModel(1, cm.getColorModel());
		// mp1.setColorModel(2, cm.getColorModel());
		// mp1.setColorModel(3, cm.getColorModel());
		mpT.setColorBar("Category Number");

		// plot
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				// plotStacked(shortData, amplitudeStack, phase);
				plotFrequencies(shortData, frtrtiamp);
			}
		});
	}

	public static void plotFrequencies(float[][] orig, float[][][] amp) {
		float fontsize = 14;
		cdp = new Sampling(totaltraces, 1.0, startCDP);
		time = new Sampling(totalsamples, dt, startSample);
		amp[0][amp[0].length - 1][amp[0][0].length - 1] = max(amp);
		amp[1][amp[0].length - 1][amp[0][0].length - 1] = max(amp);
		// amp[2][amp[0].length-1][amp[0][0].length-1] = max(amp);

		ColorMap cm = new ColorMap(ColorMap.HUE_BLUE_TO_RED);
		cm.setValueRange(min(amp), max(amp));

		float[] freq = div(
				CenterFreqs.calcFloatCenterFreqs(fmin, fmax, dt, numfilters),
				(float) dt);

		// MultiplePlot mp1 = new MultiplePlot(time, cdp, orig,cdp, amp[0], cdp,
		// amp[1], cdp, amp[2]);
		MultiplePlot mp1 = new MultiplePlot(time, cdp, orig, cdp, amp[0], cdp,
				amp[1]);
		mp1.setNearestInterpolation(1);
		mp1.setTitle("Morlet Amplitudes Separated by Frequency");
		mp1.setYLabel("Time (s)");
		mp1.setXLabel(0, "CDP Number for Original Data");
		mp1.setXLabel(1, "CDP Number for " + freq[0] + " Hz");
		mp1.setXLabel(2, "CDP Number for " + freq[1] + " Hz");
		// mp1.setXLabel(3, "CDP Number for "+freq[2]+" Hz");
		mp1.setFont(fontsize);
		mp1.setColorModel(1, cm.getColorModel());
		mp1.setColorModel(2, cm.getColorModel());
		// mp1.setColorModel(3, cm.getColorModel());
		mp1.setColorBar("Amplitude");

	}

	public static void plotStacked(float[][] origdata, float[][] amplitude,
			float[][] phase) {

		cdp = new Sampling(totaltraces, 1.0, startCDP);
		time = new Sampling(totalsamples, dt, startSample);
		freq = MorletTransform.getLogFrequencySampling(dt);
		System.out.println("*****************");
		float fontsize = 14;

		IndexColorModel cmPhase = ColorMap.getHue(0.0, 1.0);
		IndexColorModel cmAmplitude = ColorMap.getHueBlueToRed();

		SinglePlot sp = new SinglePlot(origdata);
		sp.setTitle("Hi");

		MultiplePlot mp1 = new MultiplePlot(time, cdp, origdata, freq,
				amplitude);
		mp1.setColorModel(1, cmAmplitude);
		mp1.setNearestInterpolation(1);
		mp1.setColorBar("Amplitude");
		mp1.setTitle("Gabor-Morlet Amplitude Scalogram");
		mp1.setYLabel("Time (s)");
		mp1.setXLabel(0, "Amplitude");
		mp1.setXLabel(1, "log10[Frequency(Hz)]");
		mp1.createPNG(1000, 5,
				"Amplitude, Viking-Graben, 5-85 Hz, 1500 samples.png");
		mp1.setFont(fontsize);

		MultiplePlot mp2 = new MultiplePlot(time, cdp, origdata, freq, phase);
		mp2.setColorModel(1, cmPhase);
		mp2.setNearestInterpolation(1);
		mp2.setColorBar("Phase (degrees)");
		mp2.setTitle("Gabor-Morlet Phase Scalogram");
		mp2.setYLabel("Time (s)");
		mp2.setXLabel(0, "Amplitude");
		mp2.setXLabel(1, "log10[Frequency(Hz)]");
		mp2.createPNG(1000, 5,
				"Phase, Viking-Graben, 5-85 Hz, 1500 samples.png");
		mp2.setFont(fontsize);

	}

	private static int numfilters, totalsamples, totaltraces, startCDP, endCDP,
			startSample, endSample;
	private static double dt, fmin, fmax;
	private static Sampling cdp, time, freq;
	private static float[][][] conAmp, frtrtiamp;
}
