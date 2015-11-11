package somseis;

import java.io.File;
import java.util.Random;

import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.SimplePlot;

//import gui.UIOccurenceCounter;

public class SOMBenchMark1 {
	public static void main(String[] args) {
		/*
		 * File f0 = new File("C://Users/Chris/Documents/!gb.seismic.binary");
		 * File f1 = new
		 * File("C://Users/Chris/workspace/Harmonic_Attributes/sweep1.binary");
		 * File f2 = new
		 * File("C://Users/Chris/workspace/Harmonic_Attributes/impulseat150.binary"
		 * ); String pf0 = f0.getPath(); String pf1 = f1.getPath(); String pf2 =
		 * f2.getPath();
		 * 
		 * System.out.println(f2.getAbsolutePath()); UIOccurenceCounter ui = new
		 * UIOccurenceCounter(pf0); float[][][] freqAmp = ui.getFreqs();
		 */

		float[][][] randomData = new float[2500][1500][2];
		Random r = new Random();
		float x, y;
		for (int i = 0; i < 2500; ++i) {
			for (int j = 0; j < 1500; ++j) {
				/*
				 * x = r.nextFloat(); y = r.nextFloat(); if (x > 0.5f && y
				 * <=-2*x+2) { randomData[i][j][0] = x; randomData[i][j][1] = y;
				 * } else if (x <= 0.5f && y <= 2*x) { randomData[i][j][0] = x;
				 * randomData[i][j][1] = y; }
				 */
				randomData[i][j][0] = r.nextFloat();
				randomData[i][j][1] = r.nextFloat();
			}
		}

		int numIter = 100000;
		int numAttributes = 2;
		int hLength = 12;
		SOM1 som1 = new SOM1(numIter, hLength, numAttributes);
		som1.train(randomData);
		float[][] nodes = som1.getNodes();

		float[] weightr1 = new float[hLength];
		float[] weightr2 = new float[hLength];
		for (int i = 0; i < hLength; ++i) {
			weightr1[i] = nodes[i][0];
			weightr2[i] = nodes[i][1];
		}
		// String iterations = Integer.toString(t);
		PointsView pointsr = new PointsView(weightr1, weightr2);

		pointsr.setMarkColor(java.awt.Color.RED);
		pointsr.setMarkStyle(PointsView.Mark.CROSS);
		pointsr.setMarkSize(4);
		PlotPanel pp = new PlotPanel();
		pp.addTiledView(pointsr);
		pp.setLimits(0, 0, 1, 1);

		pp.addTitle("100000");
		// pp.setLimits(0, 0, 1, 1);
		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);
		/*
		 * String iterations = Integer.toString(numIter); PointsView points =
		 * new PointsView(weight1,weight2);
		 * points.setMarkColor(java.awt.Color.RED);
		 * points.setMarkStyle(PointsView.Mark.CROSS); points.setMarkSize(4);
		 * points.setLineStyle(PointsView.Line.NONE); PlotPanel pp = new
		 * PlotPanel(); pp.addTiledView(points); pp.addTitle(iterations);
		 * pp.setLimits(0, 0, 1, 1); PlotFrame pf = new PlotFrame(pp);
		 * pf.setVisible(true);
		 */
	}
}
