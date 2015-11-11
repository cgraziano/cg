package run;

import java.awt.Color;
import java.awt.image.IndexColorModel;
import java.io.File;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanelPixels3;

import grabbers.DataGrabber;

public class Plot3Frequencies {
	public static void main(String[] args) {

		int nf = 3;
		File fileF0 = new File(
				"C:/Users/Chris/Data/FrequencyData/frequency0.dat");
		File fileF1 = new File(
				"C:/Users/Chris/Data/FrequencyData/frequency1.dat");
		File fileF2 = new File(
				"C:/Users/Chris/Data/FrequencyData/frequency2.dat");

		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");

		int sliX = 150;
		int sliIn = 80;
		int sliTime = 200;

		float[][][] seisData = DataGrabber.grab3DDataFromFile(file.getPath(),
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

		Sampling xLine = new Sampling(n3, .025f, 0);
		Sampling iLine = new Sampling(n2, .025f, 0);
		Sampling depth = new Sampling(n1, .004f, .6);// in km

		float[][] seisTimeSlice = new float[n3][n2];
		float[][] seisxLineSlice = new float[n3][n1];
		float[][] seisiLineSlice = new float[n2][n1];

		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				seisTimeSlice[i3][i2] = seisData[i3][i2][sliTime];
			}
		}

		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i1 = 0; i1 < n1; ++i1) {
				seisxLineSlice[i3][i1] = seisData[i3][sliX][i1];
			}
		}

		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i1 = 0; i1 < n1; ++i1) {
				seisiLineSlice[i2][i1] = seisData[sliIn][i2][i1];
			}
		}

		PixelsView pvT = new PixelsView(iLine, xLine, seisTimeSlice);
		pvT.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
		ColorMap cm = pvT.getColorMap();
		IndexColorModel icm = cm.setAlpha(cm.getColorModel(), 127);
		pvT.setColorModel(icm);

		PixelsView pvX = new PixelsView(depth, xLine, seisxLineSlice);
		pvX.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pvX.setColorModel(icm);

		PixelsView pvI = new PixelsView(depth, iLine, seisiLineSlice);
		pvI.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pvI.setColorModel(icm);

		PlotPanelPixels3 ppp31 = new PlotPanelPixels3(
				PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
				PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM, depth, iLine,
				xLine, freqAll);
		ppp31.setLineColor(Color.YELLOW);
		ppp31.addTiledView(0, 0, pvT);
		ppp31.addTiledView(1, 1, pvX);
		ppp31.addTiledView(1, 0, pvI);
		ppp31.getMosaic().setHeightElastic(0, 35);
		ppp31.getMosaic().setWidthElastic(1, 30);
		ppp31.getMosaic().setHeightElastic(1, 70);
		ppp31.getMosaic().setWidthElastic(0, 70);
		ppp31.setVInterval(0, 2);
		ppp31.setHInterval(1, 2);
		ppp31.setHInterval(0, 2);
		ppp31.setInterpolation(PixelsView.Interpolation.NEAREST);

		ppp31.setSlices(sliTime, sliX, sliIn);
		ppp31.setLabel2("Crossline (km)");
		ppp31.setLabel3("Inline (km)");
		ppp31.setLabel1("Depth (km)");
		ppp31.addColorBar("Blank");
		PlotFrame pf1 = new PlotFrame(ppp31);
		pf1.setFontSizeForSlide(.9, .9);
		pf1.setVisible(true);
		pf1.paintToPng(1000, 5, "3Frequencies8RHz21.9GHz60BHzSlide.png");

	}
}
