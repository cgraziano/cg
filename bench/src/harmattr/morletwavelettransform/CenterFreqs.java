package morletwavelettransform;

import static java.lang.Math.*;

/**
 * Creates the center frequencies to be used in the Morlet Transform.
 * @author Chris Graziano
 *
 */
public class CenterFreqs {
	/**
	 * 
	 * @param freqmin
	 *            The minimum bandwidth set by user (Hertz)
	 * @param freqmax
	 *            The maximum bandwidth set by user (Hertz)
	 * @param dt
	 *            The sampling interval of the input data (seconds/sample)
	 * @param numfilters
	 *            The number of filters set by user.
	 * @return A 1D array of floats that has a length equal to the number of
	 *         filters and contains logarithmically sampled frequencies. THESE
	 *         FREQUENCIES ARE IN CYCLES/SAMPLE, NOT HERTZ (CYCLES/SECOND)
	 *         ANYMORE!
	 */
	public static float[] calcFloatCenterFreqs(double freqmin, double freqmax,
			double dt, int numfilters) {

		double fmin = freqmin * dt;// fmin in (cycles/sample)
		double fmax = freqmax * dt;// fmax in (cycles/sample)

		float[] cfreq = new float[numfilters];

		// calculate mew value for logarithmic sampling (mew =
		// (fmax/fmin)^(1/numfilters))
		// Used in the equation Fcenterj = (mew^j)*Fstart
		double mew = pow(fmax / fmin, (1.0 / (double) (numfilters - 1.0)));
		for (int i = 0; i < numfilters; ++i) {
			cfreq[i] = (float) (pow(mew, i) * fmin);
		}

		return cfreq;
	}

	/**
	 * 
	 * @param freqmin
	 *            The minimum bandwidth set by user (Hertz)
	 * @param freqmax
	 *            The maximum bandwidth set by user (Hertz)
	 * @param dt
	 *            The sampling interval of the input data (seconds/sample)
	 * @param numfilters
	 *            The number of filters set by user (default is 21). Example:
	 *            Filter number (j=0(minimum frequency), 1, 2, 3... 20, 21)
	 * @return A 1D array of floats that has a length equal to the number of
	 *         filters and contains logarithmically sampled frequencies. THESE
	 *         FREQUENCIES ARE IN CYCLES/SAMPLE, NOT HERTZ (CYCLES/SECOND)
	 *         ANYMORE!
	 */
	public static double[] calcDoubleCenterFreqs(double freqmin,
			double freqmax, double dt, int numfilters) {
		double fmin = freqmin * dt;// fmin in (cycles/sample)
		double fmax = freqmax * dt;// fmax in (cycles/sample)

		double[] cfreq = new double[numfilters];

		// calculate mew value for logarithmic sampling (mew =
		// (fmax/fmin)^(1/numfilters))
		// Used in the equation Fcenterj = (mew^j)*Fstart
		double mew = pow(fmax / fmin, (1.0 / (double) (numfilters - 1.0)));
		for (int i = 0; i < numfilters; ++i) {
			cfreq[i] = (float) (pow(mew, i) * fmin);
		}

		return cfreq;
	}
}
