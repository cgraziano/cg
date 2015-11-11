package morletwavelettransform;

import static java.lang.Math.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteOrder;
import edu.mines.jtk.dsp.RecursiveGaussianFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayFile;
import edu.mines.jtk.io.ArrayOutputStream;
import static edu.mines.jtk.util.ArrayMath.cabs;
import static edu.mines.jtk.util.ArrayMath.log10;
import static edu.mines.jtk.util.ArrayMath.div;
import static edu.mines.jtk.util.ArrayMath.sqrt;

/**
 * Decomposes a signal into a number of complex sub-bands equal
 * a specified number of filters or Gabor-Morlet wavelets. The
 * amplitude and phase of these sub-bands can then be calculated.
 * 
 * @author Chris Graziano
 * 
 */
public class MorletTransform {

	/***
	 * Constructs a Gabor-Morlet Transform with a time sampling, the number of
	 * frequencies the signal will be broken up into(has to be atleast 1
	 * frequency), and the frequency band that is to be analyzed. Note: if the
	 * frequency band is 0-80 Hz and the number of frequencies is two, the
	 * frequencies that will be analyzed will be 0 and 80 Hz.
	 * 
	 * @param st
	 *            The time sampling of the signal 
	 * @param nf
	 *            The number of filters
	 * @param fmin
	 *            The minimum frequency (Hertz)
	 * @param fmax
	 *            The maximum frequency (Hertz)
	 */
	public MorletTransform(Sampling st, int nf, double fmin, double fmax) {
		centerfreqs = CenterFreqs.calcDoubleCenterFreqs(fmin, fmax,
				st.getDelta(), nf);
		this.nf = nf;

		for (int i = 0; i < nf; ++i) {
			System.out.println(centerfreqs[i] / st.getDelta() + " Hz ----- "
					+ centerfreqs[i] + " cycles/sample");
		}
		this.n = st.getCount();

		c = new float[nf][n];
		s = new float[nf][n];

		sigma = new double[nf];
		rgf = new RecursiveGaussianFilter[nf];
		double w0 = 5.336;
		double sigmaScale = w0 / (2.0 * PI);
		for (int i = 0; i < nf; ++i) {
			sigma[i] = sigmaScale / centerfreqs[i];
			rgf[i] = new RecursiveGaussianFilter(sigma[i]);
		}

		double inside = 0;
		for (int f = 0; f < nf; ++f) {

			inside = 2 * PI * centerfreqs[f];

			for (int j = 0; j < n; ++j) {
				c[f][j] = (float) cos(inside * j);
				s[f][j] = (float) sin(inside * j);
			}
		}
		xreal = new float[n];
		ximag = new float[n];
		yreal = new float[n];
		yimag = new float[n];
		scale = new float[n];

	}

	/**
	 * Applies the MorletTransform to one signal for a certain frequency.
	 * @param x The signal to be decomposed.
	 * @param f The frequency used to decompose the signal (cycles/sample)
	 * @return A complex sub-band.
	 */
	private float[] apply(float[] x, int f) {

		// create the real and imaginary parts of the input trace (step 1)
		for (int j = 0; j < n; ++j) {
			scale[j] = (float) sqrt(sqrt(PI) * 2.0 * sigma[f]);
			xreal[j] = c[f][j] * x[j] * scale[f];
			ximag[j] = -s[f][j] * x[j] * scale[f];
		}
		// convolve gaussian to x real and imaginary to create y real and
		// imaginary (step 2),
		// which gives the Gabor-Morlet wavelet its gaussian envelope.
		rgf[f].apply0(xreal, yreal);
		rgf[f].apply0(ximag, yimag);

		subbands1 = new float[2 * n];
		// Apply complex exponential that was taken out before convolution(step
		// 3)
		// This is done seperately to increase speed by taking this out of the
		// convolution, step2.
		for (int j = 0; j < n; ++j) {
			subbands1[2 * j] = (yreal[j] * c[f][j] - yimag[j] * s[f][j]);
			subbands1[2 * j + 1] = (yimag[j] * c[f][j] + yreal[j] * s[f][j]);

			

		}

		return subbands1;

	}

	/**
	 * Applies the MorletTransform to one signal for the 
	 * range of frequencies set in the constructor.
	 * @param x The signal to be decomposed.
	 * @return A set of complex sub-bands equal to 
	 * the number of filters set in the constructor.
	 */
	public float[][] apply(float[] x) {
		float[][] subbands = new float[nf][2 * n];
		for (int f = 0; f < subbands.length; ++f) {
			subbands[f] = apply(x, f);
			System.out.println("Morlet Frequeny " + f + " done");
		}
		return subbands;

	}

	/**
	 * Applies the MorletTransform to a set of signals for the 
	 * range of frequencies set in the constructor.
	 * @param x The signals to be decomposed.
	 * @return A set of complex sub-bands equal to 
	 * the number of filters set in the constructor for each signal.
	 */
	public float[][][] apply(float[][] x) {
		int n2 = x.length;
		float[][][] subbands = new float[nf][n2][n];
		for (int f = 0; f < nf; ++f) {
			for (int tr = 0; tr < n2; ++tr) {
				subbands[f][tr] = apply(x[tr], f);
			}
			System.out.println("Morlet Frequeny " + f + " done");

		}
		return subbands;

	}

	/**
	 * Applies the MorletTransform to a set of signals for the 
	 * range of frequencies set in the constructor.
	 * @param x The signals to be decomposed.
	 * @return A set of complex sub-bands equal to 
	 * the number of filters set in the constructor for each signal.
	 */
	public float[][][][] apply(float[][][] x) {
		int n3 = x.length;
		int n2 = x[0].length;
		int n1 = x[0][0].length;

		subbands4 = new float[nf][n3][n2][n1];
		for (int f = 0; f < nf; ++f) {
			for (int i3 = 0; i3 < n3; ++i3) {
				for (int tr = 0; tr < n2; ++tr) {
					subbands4[f][i3][tr] = apply(x[i3][tr], f);
				}
			}
			System.out.println("Morlet Frequeny " + f + " done");
		}
		return subbands4;
	}

	/**
	 * Applies the MorletTransform to a set of signals for the 
	 * range of frequencies set in the constructor.
	 * Sends the output of each frequency for an entire data set to a different
	 * file. This means that the first frequency of the entire data set will
	 * have its own file and that the second frequency of the entire data set
	 * will have its own file and so on and so forth up to the number of
	 * frequencies used. IN THIS APPLY METHOD, THE MAGNITUDE OF THE COMPLEX
	 * SUB-BANDS HAS ALREADY BEEN CALCULATED!
	 * 
	 * @param x
	 *            The signals to be decomposed.
	 * @param baseName
	 *            The file name
	 */
	public void applyThenWrite(float[][][] x, String baseName) {
		int n3 = x.length;
		int n2 = x[0].length;
		int n1 = x[0][0].length;

		// Create a writers to write each frequency data set out to a file
		ArrayOutputStream[] aos = new ArrayOutputStream[nf];
		for (int jf = 0; jf < nf; ++jf) {
			try {
				aos[jf] = new ArrayOutputStream(baseName + jf + ".dat",
						ByteOrder.BIG_ENDIAN);
			} catch (IOException e) {

			}
		}
		for (int jf = 0; jf < nf; ++jf) {
			try {

				for (int i3 = 0; i3 < n3; ++i3) {
					for (int tr = 0; tr < n2; ++tr) {
						subbands1 = apply(x[i3][tr], jf);
						aos[jf].writeFloats(cabs(subbands1));
					}
				}
				aos[jf].close();
			} catch (IOException e) {
				System.out.println("IOProblem");
			}
			System.out.println("Morlet Frequeny " + jf + " done");

		}

	}

	/**
	 * Constructs a new sampling of the log of the center frequencies. The
	 * frequencies are in cycles/sample. 
	 * @param dt The sampling interval of the signal(s)
	 * @return A sampling of the logarithmically sampled frequencies.
	 */
	public static Sampling getLogFrequencySampling(double dt) {
		double[] cfHz = div(centerfreqs, dt);
		// get log of centerfreqs
		int nf = cfHz.length;
		double ff = log10(cfHz[0]);
		double lf = log10(cfHz[nf - 1]);
		double df = (lf - ff) / (nf - 1);

		Sampling slogcf = new Sampling(nf, df, ff);
		return slogcf;

	}

	/**
	 * Writes the sampling of the frequencies to it's own file. 
	 * @param file The file name.
	 * @param s The sampling to write to the file.
	 */
	public static void writeSampling(File file, Sampling s) {
		ArrayFile af;
		float[] sa = { s.getCount(), (float) s.getDelta(), (float) s.getFirst() };
		try {
			af = new ArrayFile(file, "rw", ByteOrder.LITTLE_ENDIAN,
					ByteOrder.BIG_ENDIAN);
			try {
				af.writeFloats(sa);
				System.out.println("Morlet Sampling Written");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} catch (FileNotFoundException e) {
			System.out.println("FileNotFound");
		}
	}

	private float[] subbands1;
	private float[][][][] subbands4;
	private double[] sigma;
	private int n, nf;
	private double a;
	private float[][] c, s;
	private RecursiveGaussianFilter[] rgf;
	private static double[] centerfreqs;
	float[] xreal;
	float[] ximag;

	private float[] yreal;
	private float[] yimag;
	private float[] scale;

}
