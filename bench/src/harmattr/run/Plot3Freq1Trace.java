package run;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.swing.SwingUtilities;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;

import grabbers.DataGrabber;

public class Plot3Freq1Trace {
	public static void main(String[] args) {
		int nf = 3;
		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");
		File fileF0 = new File(
				"C:/Users/Chris/Data/FrequencyData/frequency0.dat");
		File fileF1 = new File(
				"C:/Users/Chris/Data/FrequencyData/frequency1.dat");
		File fileF2 = new File(
				"C:/Users/Chris/Data/FrequencyData/frequency2.dat");

		File mfile = new File(
				"C:/Users/Chris/Data/FrequencyData/MorletSamp.dat");
		ArrayInputStream ais = null;
		try {
			ais = new ArrayInputStream(mfile);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		float[] mSampling = new float[3];
		try {
			ais.readFloats(mSampling);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		Sampling sMorlet = new Sampling((int) mSampling[0],
				(double) mSampling[1], mSampling[2]);
		System.out.println("count = " + sMorlet.getCount() + " delta = "
				+ sMorlet.getDelta() + " first = " + sMorlet.getFirst());

		float[][][] rawData = DataGrabber.grab3DDataFromFile(file.getPath(),
				161, 357, 401);
		float[][][] freq0 = DataGrabber.grab3DDataFromFile(fileF0.getPath(),
				161, 357, 401);
		float[][][] freq1 = DataGrabber.grab3DDataFromFile(fileF1.getPath(),
				161, 357, 401);
		float[][][] freq2 = DataGrabber.grab3DDataFromFile(fileF2.getPath(),
				161, 357, 401);
		int n3 = freq0.length;
		int n2 = freq0[0].length;
		int n1 = freq0[0][0].length;

		float[][][][] freqAll = new float[nf][n3][n2][n1];
		freqAll[0] = freq0;
		freqAll[1] = freq1;
		freqAll[2] = freq2;

		float[][] ampTrace = new float[nf][n1];
		float[] trace = new float[n1];
		for (int f = 0; f < nf; ++f) {
			ampTrace[f] = freqAll[f][80][150];
		}
		for (int i = 0; i < n1; ++i) {
			trace[i] = rawData[80][150][i];
		}
		System.out.println("Inline = " + 80 * .025 + " km");
		System.out.println("Crossline = " + 150 * .025 + " km");

		
		Sampling depth = new Sampling(n1, .004f, .6);// in km

		
		PlotPanel pp = new PlotPanel(1, 2, PlotPanel.Orientation.X1DOWN_X2RIGHT);
		pp.addPoints(0, 0, depth, trace);

		PixelsView pv = new PixelsView(depth, sMorlet, ampTrace);
		pv.setColorModel(ColorMap.HUE_RED_TO_BLUE);
		pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv.setInterpolation(PixelsView.Interpolation.NEAREST);
		pp.addTiledView(0, 1, pv);
		pp.setVLabel("Depth (km)");
		pp.setHLabel(0, "Seismic Amplitude");
		pp.setHLabel(1, "Wave Number" + "\n" + " (cycles/km)");
		pp.addColorBar("Morlet Amplitude");
		pp.getMosaic().setWidthElastic(0, 50);
		pp.getMosaic().setWidthElastic(1, 60);

		pp.setHInterval(1, 2);

		PlotFrame pf = new PlotFrame(pp);
		pf.setFontSizeForSlide(.9, .9);
		pf.setVisible(true);
		pf.paintToPng(1000, 5,
				"SingleTraceiL80xL150FrequencyDecomposition8Hz21.9Hz60.0Hz");
		System.out.println("Paint Complete");

	}
}
