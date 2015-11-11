package somseis;

import java.io.File;
import java.util.Random;

import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.SimplePlot;
//import gui.UIOccurenceCounter;

public class SOMBenchMark {
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
		float x;
		float y;
		for (int i = 0; i < 2500; ++i) {
			for (int j = 0; j < 1500; ++j) {
				x = r.nextFloat();
				y = r.nextFloat();
				if (x <= .5 && y <= 2 * x) {
					randomData[i][j][0] = x;
					randomData[i][j][1] = y;
				} else if (x > .5 && y <= (-2 * x + 2)) {
					randomData[i][j][0] = x;
					randomData[i][j][1] = y;
				}

			}
		}

		int numIter = 100000;
		int numAttributes = 2;
		int hLength = 8;
		int vLength = 8;
		SOM2 som = new SOM2(numIter, hLength, vLength, numAttributes);
		som.train2D(randomData);
		float[][][] nodes = som.getNodes();

		float[][] weightr1 = new float[hLength][vLength];
		float[][] weightr2 = new float[hLength][vLength];
		float[][] weightc1 = new float[vLength][hLength];
		float[][] weightc2 = new float[vLength][hLength];
		for (int i = 0; i < vLength; ++i) {
			for (int j = 0; j < hLength; ++j) {
				weightc1[i][j] = nodes[i][j][0];
				weightc2[i][j] = nodes[i][j][1];
				weightr1[j][i] = nodes[i][j][0];
				weightr2[j][i] = nodes[i][j][1];
			}
		}

		String iterations = Integer.toString(numIter);
		PointsView pointsr = new PointsView(weightr1, weightr2);
		PointsView pointsc = new PointsView(weightc1, weightc2);

		pointsr.setMarkColor(java.awt.Color.RED);
		pointsr.setMarkStyle(PointsView.Mark.CROSS);
		pointsr.setMarkSize(4);
		pointsc.setMarkColor(java.awt.Color.RED);
		pointsc.setMarkStyle(PointsView.Mark.CROSS);
		pointsc.setMarkSize(4);
		PlotPanel pp = new PlotPanel();
		pp.addTiledView(pointsr);
		pp.addTiledView(pointsc);

		pp.addTitle(iterations);
		pp.setLimits(-4, -4, 4, 4);
		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);

	}
}
