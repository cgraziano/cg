package display;

import java.awt.image.ColorModel;
import java.awt.image.IndexColorModel;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.awt.ColorMapped;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.ColorBar;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;

/**
 * Plots multiple plots of 1D or 2D data using the same y axis.
 * 
 * @author Chris
 * 
 */
public class MultiplePlot {
	/**
	 * Constructs one PlotFrame with 3 TiledViews that share the same y axis,
	 * but different x axes. The default color bar is a grey scale, with white
	 * being the highest value.
	 * 
	 * @param ydata1
	 *            - Y-axis sampling for all plots
	 * @param xdata1
	 *            - X-axis sampling for data1.
	 * @param data1
	 *            - 1st set of 2D data. The first dimension: X - axis, the
	 *            second dimension: Y - axis.
	 * @param xdata2
	 *            - X-axis sampling for data2.
	 * @param data2
	 *            - 2nd set of 2D data. The first dimension: X - axis, the
	 *            second dimension: Y - axis.
	 * @param xdata3
	 *            - X-axis sampling for data3.
	 * @param data3
	 *            - 3rd set of 2D data. The first dimension: X - axis, the
	 *            second dimension: Y - axis.
	 */
	public MultiplePlot(Sampling ydata1, Sampling xdata1, float[][] data1,
			Sampling xdata2, float[][] data2, Sampling xdata3, float[][] data3) {
		pv = new PixelsView[3];
		// construct the 3 pixel views to be in one PlotPanel.
		pv[0] = new PixelsView(ydata1, xdata1, data1);
		pv[0].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv[1] = new PixelsView(ydata1, xdata2, data2);
		pv[1].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv[2] = new PixelsView(ydata1, xdata3, data3);
		pv[2].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);

		// constructs one plot panel to hold 3 plots.
		pp = new PlotPanel((int) 1, (int) 3);
		pp.addTiledView((int) 0, (int) 0, pv[0]);// The first parameter is the
													// row, the second is the
		pp.addTiledView((int) 0, (int) 1, pv[1]);// column, the third is the
													// specific PixelView
		pp.addTiledView((int) 0, (int) 2, pv[2]);// that will be in that
													// TiledView.

		pf = new PlotFrame(pp);
		pf.setVisible(true);
	}

	public MultiplePlot(Sampling ydata1, Sampling xdata1, float[][] data1,
			Sampling xdata2, float[][] data2, Sampling xdata3, float[][] data3,
			Sampling xdata4, float[][] data4) {
		pv = new PixelsView[4];
		// construct the 3 pixel views to be in one PlotPanel.
		pv[0] = new PixelsView(ydata1, xdata1, data1);
		pv[0].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv[1] = new PixelsView(ydata1, xdata2, data2);
		pv[1].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv[2] = new PixelsView(ydata1, xdata3, data3);
		pv[2].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv[3] = new PixelsView(ydata1, xdata4, data4);
		pv[3].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);

		pp = new PlotPanel((int) 1, (int) 4);
		pp.addTiledView((int) 0, (int) 0, pv[0]);// The first parameter is the
													// row, the second is the
		pp.addTiledView((int) 0, (int) 1, pv[1]);// column, the third is the
													// specific PixelView
		pp.addTiledView((int) 0, (int) 2, pv[2]);// that will be in that
													// TiledView.
		pp.addTiledView((int) 0, (int) 3, pv[3]);// that will be in that
													// TiledView.

		pf = new PlotFrame(pp);
		pf.setVisible(true);
	}

	/**
	 * Constructs one PlotFrame with 3 TiledViews that share the same y axis,
	 * but different x axes. The default color bar is a grey scale, with white
	 * being the highest value.
	 * 
	 * @param ydata1
	 *            - Y-axis sampling for all plots
	 * @param data1
	 *            - 1st set of 1D data. The first dimension: Y - axis.
	 * @param xdata2
	 *            - X-axis sampling for data2.
	 * @param data2
	 *            - 1st set of 2D data. The first dimension: X - axis, the
	 *            second dimension: Y - axis.
	 * @param xdata3
	 *            - X-axis sampling for data3.
	 * @param data3
	 *            - 2nd set of 2D data. The first dimension: X - axis, the
	 *            second dimension: Y - axis.
	 */
	public MultiplePlot(Sampling ydata1, float[] data1, Sampling xdata2,
			float[][] data2, Sampling xdata3, float[][] data3) {
		pv = new PixelsView[3];
		// Only 2 PixelsViews are created because a PixelView cannot deal with
		// 1D data
		pv[1] = new PixelsView(ydata1, xdata2, data2);
		pv[2] = new PixelsView(ydata1, xdata3, data3);

		// constructs one plot panel to hold 3 plots.
		pp = new PlotPanel((int) 1, (int) 3);
		pp.addPoints((int) 0, (int) 0, ydata1, data1);
		pp.addTiledView((int) 0, (int) 1, pv[1]);
		pp.addTiledView((int) 0, (int) 2, pv[2]);

		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);
	}

	/**
	 * Constructs one PlotFrame with 2 TiledViews that share the same y axis,
	 * but different x axes. The default color bar is a grey scale, with white
	 * being the highest value.
	 * 
	 * @param ydata1
	 *            - Y-axis sampling for all plots
	 * @param data1
	 *            - 1st set of 1D data. The first dimension: Y - axis.
	 * @param xdata2
	 *            - X-axis sampling for data2.
	 * @param data2
	 *            - 1st set of 2D data. The first dimension: X - axis, the
	 *            second dimension: Y - axis.
	 */
	public MultiplePlot(Sampling ydata1, float[] data1, Sampling xdata2,
			float[][] data2) {
		pv = new PixelsView[2];
		// construct the 1 pixel views to be in one PlotPanel.
		pv[1] = new PixelsView(ydata1, xdata2, data2);
		pv[1].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);

		// constructs one plot panel to hold 3 plots.
		pp = new PlotPanel((int) 1, (int) 2,
				PlotPanel.Orientation.X1DOWN_X2RIGHT);
		pp.addPoints((int) 0, (int) 0, ydata1, data1);
		pp.addTiledView((int) 0, (int) 1, pv[1]);

		pf = new PlotFrame(pp);
		pf.setVisible(true);
	}

	/**
	 * Constructs one PlotFrame with 2 TiledViews that share the same y axis,
	 * but different x axes. The default color bar is a grey scale, with white
	 * being the highest value.
	 * 
	 * @param ydata1
	 *            - Y-axis sampling for all plots
	 * @param xdata1
	 *            - X-axis sampling for data1.
	 * @param data1
	 *            - 1st set of 2D data. The first dimension: Y - axis, the
	 *            second dimension: Y - axis.
	 * @param xdata2
	 *            - X-axis sampling for data2.
	 * @param data2
	 *            - 2nd set of 2D data. The first dimension: X - axis, the
	 *            second dimension: Y - axis.
	 */
	public MultiplePlot(Sampling ydata1, Sampling xdata1, float[][] data1,
			Sampling xdata2, float[][] data2) {
		pv = new PixelsView[2];
		if (data1.length == 1) {
			// pv[0] = new PixelsView(ydata1, xdata1, data1);
			// pv[0].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
			pv[1] = new PixelsView(ydata1, xdata2, data2);
			pv[1].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);

			// constructs one plot panel to hold 3 plots.
			pp = new PlotPanel((int) 1, (int) 2,
					PlotPanel.Orientation.X1DOWN_X2RIGHT);
			pp.addPoints((int) 0, (int) 0, ydata1, data1[0]);
			// pp.addTiledView((int)0,(int)0, pv[0]);
			pp.addTiledView((int) 0, (int) 1, pv[1]);
			pf = new PlotFrame(pp);
			pf.setVisible(true);
			// pf.paintToPng(1000, 5,
			// "Classified Data 2 Frequencies, 100,000 iterations, 8x8 SOM, Color Scheme not dealt with");
		} else {
			// constructs the 2 pixel views to be in one PlotPanel.
			pv[0] = new PixelsView(ydata1, xdata1, data1);
			pv[0].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
			pv[1] = new PixelsView(ydata1, xdata2, data2);
			pv[1].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);

			// constructs one plot panel to hold 3 plots.
			pp = new PlotPanel((int) 1, (int) 2,
					PlotPanel.Orientation.X1DOWN_X2RIGHT);
			pp.addTiledView((int) 0, (int) 0, pv[0]);
			pp.addTiledView((int) 0, (int) 1, pv[1]);

			pf = new PlotFrame(pp);
			pf.setVisible(true);
		}
		/*
		 * // The menu bar. JMenu fileMenu = new JMenu("File");
		 * fileMenu.setMnemonic('F'); fileMenu.add(new
		 * SaveAsPngAction(_plotFrame)).setMnemonic('a'); fileMenu.add(new
		 * ExitAction()).setMnemonic('x'); JMenuBar menuBar = new JMenuBar();
		 * menuBar.add(fileMenu); _plotFrame.setJMenuBar(menuBar);
		 */

	}

	/**
	 * Constructs one PlotFrame with 2 TiledViews that share the same y axis,
	 * but different x axes. The default color bar is a grey scale, with white
	 * being the highest value.
	 * 
	 * @param ydata1
	 *            - Y-axis sampling for all plots
	 * @param xdata1
	 *            - X-axis sampling for data1.
	 * @param data1
	 *            - 1st set of 2D data. The first dimension: Y - axis, the
	 *            second dimension: Y - axis.
	 * @param xdata2
	 *            - X-axis sampling for data2.
	 * @param data2
	 *            - 2nd set of 2D data. The first dimension: X - axis, the
	 *            second dimension: Y - axis.
	 */
	public MultiplePlot(Sampling ydata1, Sampling xdata1, float[][] data1,
			Sampling xdata2, float[][] data2, ColorMap cm, float opacity) {
		pv = new PixelsView[2];
		// if (data1.length == 1){

		IndexColorModel icm = cm.getColorModel();
		// ColorMap cmT = new ColorMap(ColorMap.setAlpha(icm, opacity));
		pv[0] = new PixelsView(ydata1, xdata1, data1);
		pv[0].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv[1] = new PixelsView(ydata1, xdata2, data2);
		pv[1].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv[1].setColorModel(cm.getColorModel());

		// constructs one plot panel to hold 3 plots.
		pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		pp.addTiledView(pv[0]);
		pp.addTiledView(pv[1]);
		pf = new PlotFrame(pp);
		pf.setVisible(true);
		// pf.paintToPng(1000, 5,
		// "Classified Data 2 Frequencies, 100,000 iterations, 8x8 SOM, Color Scheme not dealt with");
		// }
		// else {
		/*
		 * //constructs the 2 pixel views to be in one PlotPanel. pv[0] = new
		 * PixelsView(ydata1, xdata1, data1);
		 * pv[0].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT); pv[1] =
		 * new PixelsView(ydata1, xdata2, data2);
		 * pv[1].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		 * 
		 * //constructs one plot panel to hold 3 plots. pp = new
		 * PlotPanel((int)1, (int)2, PlotPanel.Orientation.X1DOWN_X2RIGHT);
		 * pp.addTiledView((int)0,(int)0,pv[0]);
		 * pp.addTiledView((int)0,(int)1,pv[1]);
		 * 
		 * pf = new PlotFrame(pp); pf.setVisible(true); }
		 */
		/*
		 * // The menu bar. JMenu fileMenu = new JMenu("File");
		 * fileMenu.setMnemonic('F'); fileMenu.add(new
		 * SaveAsPngAction(_plotFrame)).setMnemonic('a'); fileMenu.add(new
		 * ExitAction()).setMnemonic('x'); JMenuBar menuBar = new JMenuBar();
		 * menuBar.add(fileMenu); _plotFrame.setJMenuBar(menuBar);
		 */

	}

	public MultiplePlot(float[][][] data2) {
		PixelsView pv1 = new PixelsView(data2);
		pv1.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		pp.addTiledView(pv1);
		pf = new PlotFrame(pp);
		pf.setVisible(true);
	}

	public void setColorModel(int n, IndexColorModel cm) {
		pv[n].setColorModel(cm);
	}

	public void setColorBar(String label) {
		pp.addColorBar(label);
	}

	public void setYLabel(String label) {
		pp.setVLabel(label);
	}

	public void setXLabel(int n, String label) {
		pp.setHLabel(n, label);
	}

	public void setTitle(String title) {
		pp.setTitle(title);
	}

	public void setFont(float size) {
		pf.setFontSize(size);
	}

	/**
	 * The data inbetween data points will have been linearly interpolated.
	 */
	public void setLinearInterpolation(int n) {
		pv[n].setInterpolation(PixelsView.Interpolation.LINEAR);
	}

	/**
	 * The data inbetween data points will have been nearest neighbor
	 * interpolated.
	 */
	public void setNearestInterpolation(int n) {
		pv[n].setInterpolation(PixelsView.Interpolation.NEAREST);
	}

	public void createPNG(double dpi, double width, String fileName) {
		pf.paintToPng(dpi, width, fileName);
	}

	private PixelsView[] pv;
	private PlotPanel pp;
	private PlotFrame pf;
	private ColorBar cb;
}
