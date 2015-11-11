package run;

import static edu.mines.jtk.util.ArrayMath.div;
import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.mosaic.*;

import static edu.mines.jtk.util.ArrayMath.max;
import static edu.mines.jtk.util.ArrayMath.min;

import java.awt.Color;
import java.awt.image.IndexColorModel;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteOrder;

import javax.swing.SwingUtilities;

import color2D.ColorMap2D;

import jogamp.graph.math.plane.Crossing;

import morletwavelettransform.CenterFreqs;
import morletwavelettransform.MorletTransform;
import somseis.SOM2;
import attributes.SimpleAttributes;
import display.MultiplePlot;
import display.SinglePlot;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.sgl.ImagePanelGroup;
import edu.mines.jtk.sgl.SimpleFrame;
import edu.mines.jtk.util.ArrayMath;
import grabbers.DataGrabber;

public class WellLogCrossPlots {
	public static void main(String[] args) {
		// File file = new File("C://Users/Chris/Documents/!gb.seismic.binary");
		// File file = new
		// File("C://Users/Chris/workspace/Harmonic_Attributes/tpst.binary");
		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");

		File gr_file = new File("C:/Users/Chris/Data/401_4_600TPData/tpgg.dat");
		File por_file = new File("C:/Users/Chris/Data/401_4_600TPData/tpgp.dat");
		File den_file = new File("C:/Users/Chris/Data/401_4_600TPData/tpgd.dat");
		File vel_file = new File("C:/Users/Chris/Data/401_4_600TPData/tpgv.dat");

		// File file = new
		// File("C://Users/Chris/workspace/Harmonic_Attributes/sweep1.binary");
		// File file = new
		// File("C://Users/Chris/workspace/Harmonic_Attributes/impulseat150.binary");
		String sFile = file.getPath();
		String sFilegr = gr_file.getPath();
		String sFilepor = por_file.getPath();
		String sFileden = den_file.getPath();
		String sFilevel = vel_file.getPath();

		float[][][] rawData = DataGrabber.grab3DDataFromFile(sFile, 161, 357,
				401);
		float[][][] gr = DataGrabber.grab3DDataFromFile(sFilegr, 161, 357, 401);
		float[][][] por = DataGrabber.grab3DDataFromFile(sFilepor, 161, 357,
				401);
		float[][][] den = DataGrabber.grab3DDataFromFile(sFileden, 161, 357,
				401);
		float[][][] vel = DataGrabber.grab3DDataFromFile(sFilevel, 161, 357,
				401);

		int n2SOM = 4;
		int n1SOM = 4;
		int nClass = n2SOM * n1SOM;
		
		File fileClass = new File("C:/Users/Chris/Data/ClassData/Norm8to60Hz3F4x4som.dat");
		float[][][] classData = DataGrabber.grab3DDataFromFile(fileClass.getPath(), 161, 357, 401);
		
/*
		float[] grMean = mean(nClass, classData, gr);
		float[] porMean = mean(nClass, classData, por);
		float[] denMean = mean(nClass, classData, den);
		float[] velMean = mean(nClass, classData, vel);
		crossPlot("Gamma Ray", " (API)", gr, classData, grMean);
		crossPlot("Porosity", " ",por, classData, porMean);
		crossPlot("Density"," (g/cc)", den, classData, denMean);
		crossPlot("Velocity"," (km/s)", vel, classData, velMean);
		*/
		//Plotting well logs on seisimc
		int sliIn = 80;
		int sliX = 150;
		int sliTime = 200;

		float[][][] seisData = DataGrabber.grab3DDataFromFile(file.getPath(),
				161, 357, 401);
		float[][][] wellData = DataGrabber.grab3DDataFromFile(
				sFilegr, 161, 357, 401);

		int n3 = wellData.length;
		int n2 = wellData[0].length;
		int n1 = wellData[0][0].length;

		Sampling xLine = new Sampling(n3, .025f, 0);
		Sampling iLine = new Sampling(n2, .025f, 0);
		Sampling depth = new Sampling(n1, .004f, .6);// in km

		float[][] classTimeSlice = new float[n2][n3];
		float[][] classxLineSlice = new float[n1][n3];
		float[][] classiLineSlice = new float[n2][n1];
		
		int size = 3;
		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				classTimeSlice[i2][i3] = wellData[i3][i2][sliTime];
				/*if (classTimeSlice[i2][i3] != 0.0f) {
					for (int j1=-size; j1<=size; ++j1) {
						for (int j2=-size; j2<=size; ++j2)
							if (0 < i2+j1 && i2+j1 < n2 && 0 < i3+j2 && i3+j2 < n3) {
								classTimeSlice[i2+j1][i3+j2] = classTimeSlice[i2][i3];
								

							}
					}
				}*/
			}
		}
		
		
		

		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i1 = 0; i1 < n1; ++i1) {
				classxLineSlice[i1][i3] = wellData[i3][sliX][i1];
			}
		}

		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i1 = 0; i1 < n1; ++i1) {
				classiLineSlice[i2][i1] = wellData[sliIn][i2][i1];
			}
		}
		float[] alpha = new float[256];
		int na = 256;
		alpha[0] = 255;
		for (int ia=1; ia<na; ++ia) 
			alpha[ia] = 135;
		
		                         
		
		IndexColorModel icm1 = ColorMap.HUE_BLUE_TO_RED;
		IndexColorModel icm2 = ColorMap.setAlpha(icm1, alpha);
		ColorMap cm = new ColorMap(100,200,icm2);
		
		PixelsView pvCT = new PixelsView(xLine, iLine, classTimeSlice);
		pvCT.setInterpolation(PixelsView.Interpolation.NEAREST);
		pvCT.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pvCT.setClips(.001f, 100f);
		pvCT.setColorModel(cm.getColorModel());
		

		PixelsView pvCX = new PixelsView(xLine, depth, classxLineSlice);
		pvCX.setInterpolation(PixelsView.Interpolation.NEAREST);
		pvCX.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
		pvCX.setColorModel(cm.getColorModel());

		PixelsView pvCI = new PixelsView(depth, iLine, classiLineSlice);
		pvCI.setInterpolation(PixelsView.Interpolation.NEAREST);
		pvCI.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pvCI.setColorModel(cm.getColorModel());
		
		PlotPanelPixels3 ppp31 = new PlotPanelPixels3(
				PlotPanelPixels3.Orientation.X1DOWN_X2RIGHT,
				PlotPanelPixels3.AxesPlacement.LEFT_BOTTOM, depth, iLine,
				xLine, seisData);
		ppp31.setColorModel(ColorMap.GRAY);

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
		PlotFrame pf1 = new PlotFrame(ppp31);
		pf1.setFontSizeForSlide(.9, .9);
		pf1.setVisible(true);
		//pf1.paintToPng(1000, 5, i + "of16Class4x4Som3F8-80HzSlide.png");
		//pf1.paintToPng(1000, 5, "grWellLocations.png");


	
	}

	private static float[] mean(int nClass, float[][][] classData,
			float[][][] wellData) {
		int n1 = wellData[0][0].length;
		int n2 = wellData[0].length;
		int n3 = wellData.length;
		float nullVal = 0.0f;
		double[] classSum = new double[nClass];
		double[] classCount = new double[nClass];
		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int i1 = 0; i1 < n1; ++i1) {
					if (wellData[i3][i2][i1] != nullVal) {
						int ic = (int) classData[i3][i2][i1];
						classSum[ic] += wellData[i3][i2][i1];
						classCount[ic] += 1.0;
					}
				}
			}
		}
		float[] classMean = new float[nClass];
		for (int ic = 0; ic < nClass; ++ic) {
			classMean[ic] = (float) (classSum[ic] / classCount[ic]);
		}
		return classMean;
	}

	private static void reclassify(int nClass, float[][][] classData,
			float[][][] wellData) {
		int n1 = wellData[0][0].length;
		int n2 = wellData[0].length;
		int n3 = wellData.length;
		float[] classMean = mean(nClass, classData, wellData);
		int[] ci = new int[nClass];
		int[] cj = new int[nClass]; // sorted class index
		for (int ic = 0; ic < nClass; ++ic) {
			ci[ic] = ic;
		}
		ArrayMath.quickIndexSort(classMean, ci);
		float[] sclassMean = new float[nClass];
		for (int ic = 0; ic < nClass; ++ic) {
			cj[ci[ic]] = ic;
			sclassMean[ic] = (float) classMean[ci[ic]];
		}
		for (int i3 = 0; i3 < n3; ++i3) {
			for (int i2 = 0; i2 < n2; ++i2) {
				for (int i1 = 0; i1 < n1; ++i1) {
					classData[i3][i2][i1] = cj[(int) classData[i3][i2][i1]];
				}
			}
		}
	}

	private static void crossPlot(String typeWellData, String units, float[][][] gr,
			float[][][] classData, float[] classMean) {
		int counter = 0;
		for (int i3 = 0; i3 < gr.length; ++i3) {
			for (int i2 = 0; i2 < gr[0].length; ++i2) {
				for (int i1 = 0; i1 < gr[0][0].length; ++i1) {
					if (gr[i3][i2][i1] > 0.001)
						++counter;
				}
			}
		}
		System.out.println("counter= " + counter);
		float[] gr1D = new float[counter];
		float[] classData1D = new float[counter];
		int i = 0;
		for (int i3 = 0; i3 < gr.length; ++i3) {
			for (int i2 = 0; i2 < gr[0].length; ++i2) {
				for (int i1 = 0; i1 < gr[0][0].length; ++i1) {
					if (gr[i3][i2][i1] > 0.001) {
						gr1D[i] = gr[i3][i2][i1];
						classData1D[i] = classData[i3][i2][i1];
						++i;
					}
				}
			}
		}
		PointsView pv = new PointsView(classData1D, gr1D);
		pv.setLineStyle(PointsView.Line.NONE);
		pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
		pv.setMarkSize(4);

		PointsView pv2 = new PointsView(classMean);
		pv2.setLineStyle(PointsView.Line.NONE);


		PlotPanel pp = new PlotPanel();
		pp.addTiledView(pv);
		pp.addTiledView(pv2);
		pp.setHLabel("Class");
		pp.setVLabel(typeWellData+units);
		pp.setTitle(typeWellData + " vs Class");

		PlotFrame pf = new PlotFrame(pp);
		pf.setFontSizeForSlide(.9, .9);
		pf.setVisible(true);
		pf.paintToPng(1000, 5, typeWellData + "crossplot_69.png");
	}
	
	
	

	private static int numfilters, totalsamples, totaltraces, n1Trace, n2Trace;
	private static int sN1, eN1, n1, sN2, eN2, n2, sN3, eN3, n3, startSample,
			endSample;
	private static double dt, fmin, fmax;
	private static Sampling cdp, time, freq;
	private static float[][][][] conAmp, frtrtiamp;
}
