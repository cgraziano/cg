package color2D;

import java.awt.Color;
import java.awt.image.IndexColorModel;
import edu.mines.jtk.awt.ColorMap;

public class ColorMap2D {
	
	/**
	 * Creates a ColorMap that should be viewed in a square.
	 * -------->Red
	 * | 0  1  2  3 
	 * | 4  5  6  7  
	 * | 8  9  10 11
	 * | 12 13 14 15
	 * Green
	 * Each number is representing a different color that
	 * has some combination of the red and green color. No
	 * blue is involved in creating this color map.
	 * The numbers in the diagram that are next 
	 * to each other will have similar colors in the Color Map.
	 *
	 * @return	A 4x4 ColorMap.
	 */
	public static ColorMap get16BoxColorMap() {

		Color[] c = new Color[256];
		int count = 0;
		for (int i2 = 0; i2 < 4; ++i2) {
			for (int i1 = 0; i1 < 4; ++i1) {
				for (int i = 0; i < 16; ++i) {
					c[count] = new Color(63 + i1 * 63, 63 + i2 * 63, 0);
					++count;
				}

			}
		}
		return new ColorMap(0, 15, c);
	}

	/**
	 * Creates a ColorMap that should be viewed in a square.
	 * -------->Red
	 * | 0  1  2  3 
	 * | 4  5  6  7  
	 * | 8  9  10 11
	 * | 12 13 14 15
	 * Green
	 * Each number is representing a different color that
	 * has some combination of the red and green color. No
	 * blue is involved in creating this color map.
	 * The numbers in the diagram that are next 
	 * to each other will have similar colors in the Color Map.
	 * This color map will only let one of the numbers be visible, while
	 * the rest of the colors will be completely transparent.
	 * 
	 * 
	 * @param classNum One of the numbers displayed in the
	 * 					diagram above to be visible.
	 * @return A 4x4 ColorMap with one color visible and 
	 * 					the rest transparent.
	 */
	public static ColorMap get16BoxColorMapSingleTrans(int classNum) {

		Color[] c = new Color[256];
		float[] alpha = new float[256];
		int count = 0;
		for (int i2 = 0; i2 < 4; ++i2) {
			for (int i1 = 0; i1 < 4; ++i1) {
				for (int i = 0; i < 16; ++i) {
					c[count] = new Color(63 + i1 * 63, 63 + i2 * 63, 0);
					if (i2 * 4 + i1 == classNum)
						alpha[count] = .5f;
					else
						alpha[count] = 0;
					++count;
				}

			}
		}

		ColorMap cm = new ColorMap(0, 15, c);
		IndexColorModel icm = cm.getColorModel();
		IndexColorModel icm2 = cm.setAlpha(icm, alpha);
		ColorMap cm2 = new ColorMap(0, 15, icm2);
		return cm2;
	}
	
	/**
	 * Similar to the method get16BoxColorMap, except that the user
	 * can set how many rows and columns the ColorMap should have and 
	 * each row contains the same spectrum of colors, but the next row
	 * in the ColorMap will be slightly darker. The
	 * user can set how much darker the next row will become.
	 * 
	 * 
	 * @param n1 The number of rows.
	 * @param n2 The number of columns.
	 * @param darkBy The darkness added to the next row each time.
	 * @return A n1xn2 ColorMap
	 */
	public static ColorMap getColorMap(int n1, int n2, int darkBy) {
		Color[] color = new Color[256];

		int rowLength = 256 / n1;
		int colLength = 256 / n2;
		int clrWidth = 256 / (n1 * n2);
		int widthCount = 0;
		int[] r = new int[rowLength];
		int[] g = new int[rowLength];
		int[] b = new int[rowLength];

		for (int ir = 0; ir < n1; ++ir) {
			// Loop through one row of our 2D color map, which mimics a red to
			// blue hue color map
			for (int ic = 0; ic < rowLength; ++ic) {
				int darken = -ir * darkBy;
				int m = -(255) / (rowLength / 2);

				// This if statement controls the thickness of each color in the
				// colorbar/array.
				if (ic % clrWidth == 0) {
					// System.out.println("i = "+ic+" rowLength "+rowLength);
					if (ic <= rowLength / 2) {
						r[ic] = m * ic + 255 + darken;
						g[ic] = -m * ic + darken;

						// check to make sure none of the color values are below
						// 0 or above 255
						if (r[ic] > 255)
							r[ic] = 255;
						if (r[ic] < 0)
							r[ic] = 0;
						if (g[ic] > 255)
							g[ic] = 255;
						if (g[ic] < 0)
							g[ic] = 0;

					} else {
						g[ic] = m * ic + 510 + darken;
						b[ic] = -m * ic + -256 + darken;

						// check to make sure none of the color values are below
						// 0 or above 255
						if (g[ic] > 255)
							g[ic] = 255;
						if (g[ic] < 0)
							g[ic] = 0;
						if (b[ic] > 255)
							b[ic] = 255;
						if (b[ic] < 0)
							b[ic] = 0;
					}
				} else {
					r[ic] = r[ic - 1];
					g[ic] = g[ic - 1];
					b[ic] = b[ic - 1];
				}
				++widthCount;
				color[ir * rowLength + ic] = new Color(r[ic], g[ic], b[ic]);

				
			}

		}
		

		ColorMap cm = new ColorMap(0, n1*n2, color);
		return cm;
	}

	/**
	 * Similar to the method get16BoxColorMapSingleTrans, except that the user
	 * can set how many rows and columns the ColorMap should have and 
	 * each row contains the same spectrum of colors, but the next row
	 * in the ColorMap will be slightly darker. The
	 * user can set how much darker the next row will become.
	 * This color map will only let one of the numbers be visible, while
	 * the rest of the colors will be completely transparent.
	 * 
	 * 
	 * @param n1 The number of rows.
	 * @param n2 The number of columns.
	 * @param darkBy The darkness added to the next row each time.
	 * @param classNum One of the numbers displayed in the
	 * 					diagram above to be visible.
	 * @return A n1xn2 ColorMap with one color visible and the rest
	 * 					transparent.
	 */
	public static ColorMap getColorMapSingleTransparent(int n1, int n2, int darkBy, int classNum) {
		Color[] color = new Color[256];

		int rowLength = 256 / n1;
		int colLength = 256 / n2;
		int clrWidth = 256 / (n1 * n2);
		int widthCount = 0;
		int catCount = 0;

		int[] r = new int[rowLength];
		int[] g = new int[rowLength];
		int[] b = new int[rowLength];
		float[] alpha = new float[256];

		for (int ir = 0; ir < n1; ++ir) {
			// Loop through one row of our 2D color map, which mimics a red to
			// blue hue color map
			for (int ic = 0; ic < rowLength; ++ic) {
				int darken = -ir * darkBy;
				int m = -(255) / (rowLength / 2);

				// This if statement controls the thickness of each color in the
				// colorbar/array.
				if (ic % clrWidth == 0) {
					// System.out.println("i = "+ic+" rowLength "+rowLength);
					if (ic <= rowLength / 2) {
						r[ic] = m * ic + 255 + darken;
						g[ic] = -m * ic + darken;

						// check to make sure none of the color values are below
						// 0 or above 255
						if (r[ic] > 255)
							r[ic] = 255;
						if (r[ic] < 0)
							r[ic] = 0;
						if (g[ic] > 255)
							g[ic] = 255;
						if (g[ic] < 0)
							g[ic] = 0;

					} else {
						g[ic] = m * ic + 510 + darken;
						b[ic] = -m * ic + -256 + darken;

						// check to make sure none of the color values are below
						// 0 or above 255
						if (g[ic] > 255)
							g[ic] = 255;
						if (g[ic] < 0)
							g[ic] = 0;
						if (b[ic] > 255)
							b[ic] = 255;
						if (b[ic] < 0)
							b[ic] = 0;
					}
					++catCount;
				} else {
					r[ic] = r[ic - 1];
					g[ic] = g[ic - 1];
					b[ic] = b[ic - 1];
				}
				++widthCount;
				if (catCount - 1 == classNum) {
					alpha[ir * rowLength + ic] = .5f;
					color[ir * rowLength + ic] = new Color(r[ic], g[ic], b[ic]);
				} else {
					alpha[ir * rowLength + ic] = 0.0f;
					color[ir * rowLength + ic] = new Color(r[ic], g[ic], b[ic]);
				}

				
			}

		}
		
		double min = 0;
		double max = n1*n2;
		ColorMap cm = new ColorMap(min, max, color);
		IndexColorModel icm = cm.getColorModel();
		IndexColorModel icm2 = cm.setAlpha(icm, alpha);
		ColorMap cm2 = new ColorMap(min, max, icm2);

		return cm2;
	}

	
}
