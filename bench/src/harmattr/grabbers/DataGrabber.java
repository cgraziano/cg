package grabbers;

import java.io.DataInputStream;

import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.sgl.*;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.ByteOrder;

/**
 * A class devoted to getting binary data from a specified file and possibly
 * shortening it if necessary. The files this class needs to have any headers
 * stripped from the data file and the data file should be in the Big Endian
 * format.
 * 
 * @author Chris Graziano
 */
public class DataGrabber {

	/**
	 * Gets 2D data from a file and puts it into a 2D array of floats.
	 * 
	 * @param fileName
	 *            The file's path and name Example: File file = new
	 *            File("C://Users/Chris/Documents/Data.dat"); String fileName =
	 *            file.getPath(); fileName would then be the String used in this
	 *            method
	 * @param n2
	 *            The number of samples in the slow dimension of the data.
	 * @param n1
	 *            The number of samples in the fast dimension of the data.
	 * @return The 2D data in a 2D array of floats.
	 */
	public static float[][] grab2DDataFromFile(String fileName, int n2, int n1) {
		float[][] data = new float[n2][n1];

		try {
			ArrayFile af = new ArrayFile(fileName, "r", ByteOrder.BIG_ENDIAN,
					ByteOrder.LITTLE_ENDIAN);
			af.readFloats(data);

			System.out.println("Data extracted from " + fileName);

			af.close();

		} catch (IOException e) {
			System.out.println("File Cannot be found!");
			throw new RuntimeException(e);
		}

		return data;
	}

	/**
	 * Gets 3D data from a file and puts it into a 3D array of floats.
	 * 
	 * @param fileName
	 *            The file's path and name Example: File file = new
	 *            File("C://Users/Chris/Documents/Data.dat"); String fileName =
	 *            file.getPath(); fileName would then be the String used in this
	 *            method
	 * @param n3
	 *            The number of samples in the slowest dimension of the data.
	 * @param n2
	 *            The number of samples in the medium speed dimension of the
	 *            data.
	 * @param n1
	 *            The number of samples in the fast dimension of the data.
	 * @return The 3D data in a 3D array of floats.
	 */
	public static float[][][] grab3DDataFromFile(String fileName, int n3,
			int n2, int n1) {
		float[][][] data = new float[n3][n2][n1];

		try {
			ArrayFile af = new ArrayFile(fileName, "r", ByteOrder.BIG_ENDIAN,
					ByteOrder.LITTLE_ENDIAN);
			af.readFloats(data);

			System.out.println("Data extracted from file.");

			af.close();

		} catch (IOException e) {
			System.out.println("File Cannot be found!");
			throw new RuntimeException(e);
		}
		return data;
	}

	/**
	 * Shortens the 2D array to the specified index values
	 * 
	 * @param origArray
	 *            The array to shorten
	 * @param s2
	 *            The index of the origArray that the new array will start at in
	 *            the slow dimension.
	 * @param e2
	 *            The index of the origArray that the new array will end at in
	 *            the slow dimension.
	 * @param s1
	 *            The index of the origArray that the new array will start at in
	 *            the fast dimension.
	 * @param e1
	 *            The index of the origArray that the new array will end at in
	 *            the fast dimension.
	 * @return The shortened array.
	 */
	public static float[][] shorten2DData(float[][] origArray, int s2, int e2,
			int s1, int e1) {
		int n2 = e2 - s2 + 1;
		int n1 = e1 - s1 + 1;
		float[][] shortarray = new float[n2][n1];

		for (int i = s2, ii = 0; i < e2 && ii < n2; ++i, ++ii) {
			for (int j = s1, jj = 0; j < e1 && jj < n1; ++j, ++jj) {
				shortarray[ii][jj] = origArray[i][j];
			}
		}
		return shortarray;
	}

	/**
	 * Shortens the 3D array to the specified index values
	 * 
	 * @param origArray
	 *            The array to shorten
	 * @param S3
	 *            The index of the origArray that the new array will start at in
	 *            the slow dimension.
	 * @param e3
	 *            The index of the origArray that the new array will end at in
	 *            the slow dimension.
	 * @param s2
	 *            The index of the origArray that the new array will start at in
	 *            the second slowest dimension.
	 * @param e2
	 *            The index of the origArray that the new array will end at in
	 *            the second slowest dimension.
	 * @return The shortened array.
	 * @param s1
	 *            The index of the origArray that the new array will start at in
	 *            the fast dimension.
	 * @param e1
	 *            The index of the origArray that the new array will end at in
	 *            the fast dimension.
	 * @return The shortened array.
	 */
	public static float[][][] shorten3DData(float[][][] originalarray, int s3,
			int e3, int s2, int e2, int s1, int e1) {
		int n3 = e3 - s3 + 1;
		int n2 = e2 - s2 + 1;
		int n1 = e1 - s1 + 1;

		float[][][] shortarray = new float[n3][n2][n1];

		for (int i = s3, ii = 0; i < e3 && ii < n3; ++i, ++ii) {
			for (int j = s2, jj = 0; j < e2 && jj < n2; ++j, ++jj) {
				for (int k = s1, kk = 0; k < e1 && kk < n1; ++k, ++kk) {
					shortarray[ii][jj][kk] = originalarray[i][j][k];
				}
			}
		}

		return shortarray;
	}

}
