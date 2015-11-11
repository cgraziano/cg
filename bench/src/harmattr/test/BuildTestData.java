package test;

import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.ByteOrder;

import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.sgl.ImagePanelGroup;
import edu.mines.jtk.sgl.SimpleFrame;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Writes a binary file that mimics a 2D seismic data set. The length and the
 * number of samples can be changed.
 * 
 * @author Chris
 * 
 */
public class BuildTestData {

	/**
	 * 
	 * @param i
	 *            Index used to increment the time the sweep is on, incremented
	 *            in buildSweepData
	 * @param dt
	 *            The sampling rate
	 * @param slope
	 *            =(maxf-minf)/(maxt-mint), used to know how the sweep
	 *            increments with time
	 * @param maxf
	 *            The maximum frequency of the sweep
	 * @param mint
	 *            The minimum frequency of the sweep
	 * @return The sweep value of a cosine that is dependent on the index i,
	 *         maxf, minf, slope, and mint, but mostly i.
	 */
	public static float sweep(int i, float dt, float slope, float minf,
			float mint) {

		double t = mint + i * dt;

		return (float) cos(2.0 * PI * (minf * t + slope * t * t / 2));

	}

	/**
	 * 
	 * @param i
	 *            Index used to increment the time the sweep is on, incremented
	 *            in buildSweepData
	 * @param dt
	 *            The sampling rate
	 * @param slope
	 *            =(maxf-minf)/(maxt-mint), used to know how the sweep
	 *            increments with time
	 * @param maxf
	 *            The maximum frequency of the sweep
	 * @param mint
	 *            The minimum frequency of the sweep
	 * @return The sweep value of a cosine that is dependent on the index i,
	 *         maxf, minf, slope, and mint, but mostly i.
	 */
	public static float constF(int i, float dt, float f) {
		float mint = 0;
		double t = mint + i * dt;

		return (float) cos(2.0 * PI * (f * t));

	}

	public static float[][] build2DSweepData(int nt, int nx, float dt,
			float maxf, float minf, float mint) {
		float[][] testdata = new float[nx][nt]; // Set ny to 1 to mimic 2D data
												// from a 3D format

		float slope = (maxf - minf) / (dt * nt - mint);

		// build one array of data
		for (int t = 0; t < nt; ++t) {
			testdata[0][t] = sweep(t, dt, slope, minf, mint);
		}
		SimplePlot.asSequence(testdata[0]);
		// instead of building the same array over and over again,
		// copy the first array to the other arrays.
		for (int x = 1; x < nx; ++x) {
			for (int t = 0; t < nt; ++t) {

				testdata[x][t] = testdata[0][t];

			}
		}

		return testdata;

	}

	public static float[][] build2DConstFData(int nt, int nx, float dt, float f) {
		float[][] testdata = new float[nx][nt]; // Set ny to 1 to mimic 2D data
												// from a 3D format

		// build one array of data
		for (int t = 0; t < nt; ++t) {
			testdata[0][t] = constF(t, dt, f);
		}
		SimplePlot.asSequence(testdata[0]);
		// instead of building the same array over and over again,
		// copy the first array to the other arrays.
		for (int x = 1; x < nx; ++x) {
			for (int t = 0; t < nt; ++t) {

				testdata[x][t] = testdata[0][t];

			}
		}

		return testdata;

	}

	public static float[][] build2DConstFAmpDecData(int nt, int nx, float dt,
			float f) {
		float[][] testdata = new float[nx][nt]; // Set ny to 1 to mimic 2D data
												// from a 3D format

		float slopet = 1 / (nt);
		// build one array of data
		for (int t = 0; t < nt; ++t) {
			testdata[0][t] = ((float) (nt - t) / (float) nt) * constF(t, dt, f);
		}
		SimplePlot.asSequence(testdata[0]);
		// instead of building the same array over and over again,
		// copy the first array to the other arrays.
		for (int x = 1; x < nx; ++x) {
			for (int t = 0; t < nt; ++t) {

				testdata[x][t] = testdata[0][t];

			}
		}

		return testdata;

	}

	public static float[][][] build3DSweepData(int nt, int nx, int ny,
			float dt, float maxf, float minf, float mint) {
		float[][][] testdata = new float[nx][ny][nt]; // Set ny to 1 to mimic 2D
														// data from a 3D format

		float slope = (maxf - minf) / (dt * nt - mint);

		// build one array of data
		for (int t = 0; t < nt; ++t) {
			testdata[0][0][t] = sweep(t, dt, slope, minf, mint);
		}
		// SimplePlot.asSequence(testdata[0]);
		// instead of building the same array over and over again,
		// copy the first array to the other arrays.
		for (int x = 1; x < nx; ++x) {
			for (int y = 1; y < ny; ++y) {
				for (int t = 0; t < nt; ++t) {

					testdata[x][y][t] = testdata[0][0][t];

				}
			}
		}
		System.out.println("Data built.");

		return testdata;

	}

	public static float[][] impulseData(int nt, int nx, float dt, int imptime) {
		float[][] testdata = new float[nx][nt]; // Set ny to 1 to mimic 2D data
												// from a 3D format

		// instead of building the same array over and over again,
		// copy the first array to the other arrays.
		for (int x = 1; x < nx; ++x) {
			for (int t = 0; t < nt; ++t) {

				testdata[x][t] = 0;

			}
			testdata[x][imptime] = (float) 1.0;
		}
		SimplePlot.asPixels(testdata);

		return testdata;

	}

	/**
	 * Writes a 2D float array to a file in binary. The fileName parameter is
	 * 
	 * @param testdata
	 * @param fileName
	 */
	public static void write2DData(float[][] testdata, String fileName) {

		try {
			FileOutputStream fos = new FileOutputStream(fileName);
			DataOutputStream dos = new DataOutputStream(fos);

			ArrayFile af = new ArrayFile(fileName, "rw",
					ByteOrder.LITTLE_ENDIAN, ByteOrder.LITTLE_ENDIAN);
			af.writeFloats(testdata);
			SimplePlot.asPixels(testdata);
			af.close();

		} catch (IOException e) {
			System.out.println("File cannot be found");
			throw new RuntimeException(e);// to stop program so as to not waste
											// time.
		}

	}

	/**
	 * Writes a 3D float array to a file in binary. The fileName parameter is
	 * 
	 * @param testdata
	 * @param fileName
	 */
	public static void write3DData(float[][][] testdata, String fileName) {

		try {
			FileOutputStream fos = new FileOutputStream(fileName);
			DataOutputStream dos = new DataOutputStream(fos);

			// int n3 = testdata.length;
			// int n2 = testdata[0].length;
			// int n1 = testdata[0][0].length;

			ArrayFile af = new ArrayFile(fileName, "rw",
					ByteOrder.LITTLE_ENDIAN, ByteOrder.LITTLE_ENDIAN);
			af.writeFloats(testdata);
			ImagePanelGroup img = new ImagePanelGroup(testdata);

			SimpleFrame sf = new SimpleFrame();
			sf.addImagePanels(img);
			af.close();
			System.out.println("Data Written");

		} catch (IOException e) {
			System.out.println("File cannot be found");
			throw new RuntimeException(e);// to stop program so as to not waste
											// time.
		}

	}
}
