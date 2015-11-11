package run;

import java.awt.Color;
import java.awt.image.IndexColorModel;
import java.io.File;

import javax.swing.SwingUtilities;

import color2D.ColorMap2D;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import grabbers.DataGrabber;

public class PlotOriginalSeismic {
	public static void main(String[] args) {
		File fileClass = new File(
				"C:/Users/Chris/Data/ClassData/Norm8to60Hz3F4x4som.dat");
		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");

		final int n1SOM = 4;
		final int n2SOM = 4;

		final int sliIn = 80;
		final int sliX = 150;
		final int sliTime = 200;

		final float[][][] seisData = DataGrabber.grab3DDataFromFile(
				file.getPath(), 161, 357, 401);
		float[][][] classData = DataGrabber.grab3DDataFromFile(
				fileClass.getPath(), 161, 357, 401);

		int n3 = classData.length;
		int n2 = classData[0].length;
		int n1 = classData[0][0].length;

		final Sampling xLine = new Sampling(n3, .025f, 0);
		final Sampling iLine = new Sampling(n2, .025f, 0);
		final Sampling depth = new Sampling(n1, .004f, .6);

		final float[][] classTimeSlice = new float[n2][n3];
		final float[][] classxLineSlice = new float[n1][n3];
		final float[][] classiLineSlice = new float[n2][n1];

		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				classTimeSlice[i2][i3] = classData[i3][i2][sliTime];
			}
		}

		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i1 = 0; i1 < n1; ++i1) {
				classxLineSlice[i1][i3] = classData[i3][sliX][i1];
			}
		}

		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i1 = 0; i1 < n1; ++i1) {
				classiLineSlice[i2][i1] = classData[sliIn][i2][i1];

			}
		}

		SwingUtilities.invokeLater(new Runnable() {
			public void run() {

				
				ColorMap cm = ColorMap2D.get16BoxColorMap();
				IndexColorModel icm = cm.getColorModel();
				

				PixelsView pvCT = new PixelsView(xLine, iLine, classTimeSlice);
				pvCT.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
				pvCT.setInterpolation(PixelsView.Interpolation.NEAREST);
				pvCT.setClips(-.5f, 15.5f);
				pvCT.setColorModel(icm);

				PixelsView pvCX = new PixelsView(xLine, depth, classxLineSlice);
				pvCX.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
				pvCX.setInterpolation(PixelsView.Interpolation.NEAREST);
				pvCX.setClips(-.5f, 15.5f);
				pvCX.setColorModel(icm);

				PixelsView pvCI = new PixelsView(depth, iLine, classiLineSlice);
				pvCI.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
				pvCI.setInterpolation(PixelsView.Interpolation.NEAREST);
				pvCI.setClips(-.5f, 15.5f);
				pvCI.setColorModel(icm);

				PlotPanelPixels3 ppp31 = new PlotPanelPixels3(
						PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
						PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM, depth,
						iLine, xLine, seisData);
			
				ppp31.getMosaic().setHeightElastic(0, 35);
				ppp31.getMosaic().setWidthElastic(1, 30);
				ppp31.getMosaic().setHeightElastic(1, 70);
				ppp31.getMosaic().setWidthElastic(0, 70);
				ppp31.setVInterval(0, 2);
				ppp31.setHInterval(1, 2);
				ppp31.setHInterval(0, 2);
				ppp31.setInterpolation(PixelsView.Interpolation.NEAREST);


				ppp31.setLineColor(Color.YELLOW);

				ppp31.setSlices(sliTime, sliX, sliIn);

				ppp31.setLabel2("Crossline (km)");
				ppp31.setLabel3("Inline (km)");
				ppp31.setLabel1("Depth (km)");
				ppp31.addColorBar("Seismic Amplitude");
				PlotFrame pf1 = new PlotFrame(ppp31);
				pf1.setFontSizeForSlide(.9, .9);
				pf1.setVisible(true);
				
				pf1.paintToPng(1000, 5, "OriginalSeismicSlide.png");

				// Create Color Box
				float[][] cbox = new float[n1SOM][n2SOM];
				for (int i2 = 0; i2 < n2SOM; ++i2) {
					for (int i1 = 0; i1 < n1SOM; ++i1) {
						cbox[n2SOM - 1 - i2][i1] = i2 * n1SOM + i1;
						System.out.println(cbox[i2][i1]);
					}
				}
				PixelsView pv = new PixelsView(cbox);
				pv.setColorModel(cm.getColorModel());
				pv.setInterpolation(PixelsView.Interpolation.NEAREST);
				PlotPanel pp = new PlotPanel();
				pp.addTiledView(pv);
				PlotFrame pf2 = new PlotFrame(pp);
				pf2.setVisible(true);
				pf2.paintToPng(1000, 5, "BoxColorBar.png");

			}
		});
	}
}
