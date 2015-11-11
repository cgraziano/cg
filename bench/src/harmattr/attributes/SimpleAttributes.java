package attributes;

import static edu.mines.jtk.util.ArrayMath.cabs;
import static edu.mines.jtk.util.ArrayMath.carg;
import static edu.mines.jtk.util.ArrayMath.mul;

import static java.lang.Math.PI;

/**
 * 
 * @author Chris Graziano Calculates the amplitude or phase spectrum of a
 *         signal's complex sub-bands
 */
public class SimpleAttributes {
	/**
	 * Finds the amplitude spectrum of each subband of a single signal.
	 * 
	 * @param subbands
	 *            - First dimension: number of sub-bands. Second dimension: the
	 *            number of samples in each sub-band(twice the number of samples
	 *            in the original trace to account for real and imaginary
	 *            parts).
	 * @return The magnitude of the real and complex parts of each subband.
	 */
	public static float[][] findAmplitude(float[][] subbands) {
		int n2 = subbands.length;// The number of complex subbands.
		int n1 = subbands[0].length / 2;// The number of samples (the real and
										// imaginary part are counted
										// individually).

		float[][] amplitude = new float[n2][n1];

		for (int f = 0; f < n2; ++f) {
			amplitude[f] = cabs(subbands[f]);// Find the magnitude of each
												// complex subband.
		}

		return amplitude;
	}

	/**
	 * Finds the amplitude spectrum of the each subband of a 1D array of
	 * signals.
	 * 
	 * @param subbands
	 *            First dimension: number of sub-bands. Second dimension: number
	 *            of signals. Third dimension: the number of samples in each
	 *            sub-band(twice the number of samples in the original trace to
	 *            account for real and imaginary parts).
	 * @return The magnitude of the real and complex parts of each subband for
	 *         each signal
	 */
	public static float[][][] findAmplitude(float[][][] subbands) {
		int n3 = subbands.length;// The number of complex subbands.
		int n2 = subbands[0].length;// The number of signals.
		int n1 = subbands[0][0].length / 2;// The number of samples (the real
											// and imaginary part are counted
											// individually).

		float[][][] amplitude = new float[n3][n2][n1];

		for (int f = 0; f < n3; ++f) {
			for (int i2 = 0; i2 < n2; ++i2) {
				amplitude[f][i2] = cabs(subbands[f][i2]);// Find the magnitude
															// of each complex
															// subband.
			}
		}

		return amplitude;
	}

	/**
	 * Finds the amplitude spectrum of the each subband of a 2D array of
	 * signals.
	 * 
	 * @param subbands
	 *            First dimension: number of sub-bands. Second dimension: number
	 *            of signals in the slow d dimension of a 2D array of signals.
	 *            Third dimension: number of signals in the slow dimension of a
	 *            2D array of signals. Fourth dimension: the number of samples
	 *            in each sub-band(twice the number of samples in the original
	 *            trace to account for real and imaginary parts).
	 * @return The magnitude of the real and complex parts of each subband for
	 *         each signal
	 */
	public static float[][][][] findAmplitude(float[][][][] subbands) {
		int n4 = subbands.length;// The number of subbands.
		int n3 = subbands[0].length;// The number of signals in the slow
									// dimension of a 2D array of signals.
		int n2 = subbands[0][0].length;// The number of signals in the fast
										// dimension of a 2D array of signals.
		int n1 = subbands[0][0][0].length / 2;// The number of samples (the real
												// and imaginary part are
												// counted individually).

		float[][][][] amplitude = new float[n4][n3][n2][n1];

		for (int f = 0; f < n4; ++f) {
			for (int i3 = 0; i3 < n3; ++i3) {
				for (int i2 = 0; i2 < n2; ++i2) {
					amplitude[f][i3][i2] = cabs(subbands[f][i3][i2]);// Find the
																		// magnitude
																		// of
																		// each
																		// complex
																		// array.
				}
			}
			System.out.println("Amplitude Frequeny " + f + " done");

		}

		return amplitude;

	}

	/**
	 * Calculates the phase spectrum of the complex subbands of a single signal.
	 * 
	 * @param subbands
	 *            - First dimension: number of sub-bands. Second dimension: the
	 *            number of samples in each sub-band or twice the number of
	 *            samples in the orginal trace to account for real and imaginary
	 *            parts.
	 * @return phase - The phase of the real and imaginary parts of each
	 *         sub-band. phase[f][t] = arcTan(subbands[f][2*t +
	 *         1]/subbands[f][2*t]),
	 * 
	 */
	public static float[][] findPhase(float[][] subbands) {
		int n2 = subbands.length;// The number of subbands
		int n1 = subbands[0].length;// The number of samples (the real and
									// imaginary part are counted individually).

		float[][] phase = new float[n2][n1];

		for (int cf = 0; cf < n2; ++cf) {
			phase[cf] = mul(carg(subbands[cf]), (float) (180 / PI));
		}

		return phase;
	}

}