package display;

import java.awt.image.IndexColorModel;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;

/**
 * Plots one plot of 1D or 2D data.
 * 
 * @author Chris
 * 
 */
public class SinglePlot {

	public SinglePlot(float[][] data) {
		pv = new PixelsView(data);
		pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		pp.addTiledView(pv);
		pp.addColorBar();
		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);
	}

	/**
	 * Plots 2D data with a y and x axis sampling. The default color bar is
	 * white represents high values, while black represents low values.
	 * 
	 * @param yaxis
	 *            - Y axis uniform sampling
	 * @param xaxis
	 *            - X axis uniform sampling
	 * @param data
	 *            - 2D data. The first dimension: X - axis, the second
	 *            dimension: Y - axis.
	 * @param percmin
	 *            The lowest percentile that is displayed.
	 * @param percmax
	 *            The highest percentile that is displayed.
	 */
	public SinglePlot(Sampling yaxis, Sampling xaxis, float[][] data,
			float percmin, float percmax) {
		pv = new PixelsView(yaxis, xaxis, data);
		pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
		pv.setPercentiles(percmin, percmax);

		pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		pp.addTiledView(pv);
		pp.addColorBar();

	}

	/**
	 * Plots 1D data with a y axis sampling. The default color bar is white
	 * represents high values, while black represents low values.
	 * 
	 * @param yaxis
	 *            - Y axis uniform sampling
	 * @param data
	 *            - 1D data, the first dimension is the y axis.
	 */
	public SinglePlot(Sampling yaxis, float[] data) {
		pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		pp.addPoints(data);
		pp.addColorBar();
	}

	/**
	 * Set the Title and Axis label and the overall font size. The title font
	 * will be 1.5 times larger than the set font size.
	 * 
	 * @param title
	 *            String
	 * @param yaxis
	 *            String
	 * @param xaxis
	 *            String
	 */
	public void setTitleAndAxis(String title, String yaxis, String xaxis,
			float fontsize) {
		pp.addTitle(title);
		pp.setVLabel(yaxis);
		pp.setHLabel(xaxis);

		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);
		pf.setFontSize((float) fontsize);
	}

	/**
	 * The data inbetween data points will have been linearly interpolated.
	 */
	public void setLinearInterpolation() {
		pv.setInterpolation(PixelsView.Interpolation.LINEAR);
	}

	/**
	 * The data inbetween data points will have been nearest neighbor
	 * interpolated.
	 */
	public void setNearestInterpolation() {
		pv.setInterpolation(PixelsView.Interpolation.NEAREST);
	}

	/**
	 * Set the color bar to display the same color for Pi and -Pi.
	 */
	public void setPhaseColorBar() {
		pv.setColorModel(ColorMap.getHue(0.0, 1.0));
	}

	public void setColorModel(IndexColorModel cm) {
		pv.setColorModel(cm);
	}

	public void setColorBar(String label) {
		pp.addColorBar(label);
	}

	public void setYLabel(String label) {
		pp.setVLabel(label);
	}

	public void setXLabel(String label) {
		pp.setHLabel(label);
	}

	public void setTitle(String title) {
		pp.setTitle(title);
	}

	/**
	 * Creates a PNG image of the corresponding plot.
	 * 
	 * @param dpi
	 *            the image resolution in dots per inch.
	 * @param width
	 *            the image width in inches.
	 * @param fileName
	 *            the name of the file to contain the PNG image.
	 */
	public void createPNG(double dpi, double width, String fileName) {
		pf.paintToPng(dpi, width, fileName);
	}

	private PixelsView pv;
	private PlotPanel pp;
	private PlotFrame pf;
}
