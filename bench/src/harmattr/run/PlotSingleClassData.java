package run;

import java.awt.Color;
import java.awt.image.IndexColorModel;
import java.io.File;

import color2D.ColorMap2D;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import grabbers.DataGrabber;

public class PlotSingleClassData {
	public static void main(String[] args) {
		int nf = 3;
		File fileClass = new File(
				"C:/Users/Chris/Data/ClassData/Norm8to60Hz3F4x4som.dat");
		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");

		int n1SOM = 4;
		int n2SOM = 4;

		int sliIn = 80;
		int sliX = 150;
		int sliTime = 200;

		float[][][] seisData = DataGrabber.grab3DDataFromFile(file.getPath(),
				161, 357, 401);
		float[][][] classData = DataGrabber.grab3DDataFromFile(
				fileClass.getPath(), 161, 357, 401);

		int n3 = classData.length;
		int n2 = classData[0].length;
		int n1 = classData[0][0].length;

		Sampling xLine = new Sampling(n3, .025f, 0);
		Sampling iLine = new Sampling(n2, .025f, 0);
		Sampling depth = new Sampling(n1, .004f, .6);// in km

		float[][] classTimeSlice = new float[n2][n3];
		float[][] classxLineSlice = new float[n1][n3];
		float[][] classiLineSlice = new float[n2][n1];

		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				classTimeSlice[i2][i3] = classData[i3][i2][sliTime];
				// classxLineSlice[i3][i]
			}
		}

		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i1 = 0; i1 < n1; ++i1) {
				classxLineSlice[i1][i3] = classData[i3][sliX][i1];
				// classxLineSlice[i3][i]
			}
		}

		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i1 = 0; i1 < n1; ++i1) {
				classiLineSlice[i2][i1] = classData[sliIn][i2][i1];
				// classxLineSlice[i3][i]
			}
		}

		for (int i = 0; i < n1SOM * n2SOM; ++i) {
			
			ColorMap cm = ColorMap2D.get16BoxColorMapSingleTrans(i);
			PixelsView pvCT = new PixelsView(xLine, iLine, classTimeSlice);
			pvCT.setInterpolation(PixelsView.Interpolation.NEAREST);
			pvCT.setClips(-.5f, 15.5f);
			pvCT.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
			pvCT.setColorModel(cm.getColorModel());

			PixelsView pvCX = new PixelsView(xLine, depth, classxLineSlice);
			pvCX.setInterpolation(PixelsView.Interpolation.NEAREST);
			pvCX.setClips(-.5f, 15.5f);
			pvCX.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
			pvCX.setColorModel(cm.getColorModel());

			PixelsView pvCI = new PixelsView(depth, iLine, classiLineSlice);
			pvCI.setInterpolation(PixelsView.Interpolation.NEAREST);
			pvCI.setClips(-.5f, 15.5f);
			pvCI.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
			pvCI.setColorModel(cm.getColorModel());

			PlotPanelPixels3 ppp31 = new PlotPanelPixels3(
					PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
					PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM, depth, iLine,
					xLine, seisData);
			ppp31.addTiledView(0, 0, pvCT);
			ppp31.addTiledView(1, 1, pvCX);
			ppp31.addTiledView(1, 0, pvCI);
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
			ppp31.addColorBar("Class");
			// ppp31.setColorModel(cm.getColorModel());
			PlotFrame pf1 = new PlotFrame(ppp31);
			pf1.setFontSizeForSlide(.9, .9);
			pf1.setVisible(true);
			pf1.paintToPng(1000, 5, i + "of16Class4x4Som3F8-80HzSlide.png");
		}
		System.out.println("Paint Complete");
		
	}
}
