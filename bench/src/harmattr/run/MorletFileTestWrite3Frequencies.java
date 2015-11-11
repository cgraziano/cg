package run;

import java.io.File;
import morletwavelettransform.MorletTransform;
import edu.mines.jtk.dsp.Sampling;
import grabbers.DataGrabber;

public class MorletFileTestWrite3Frequencies {
	public static void main(String[] args) {
		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");
		String sFile = file.getPath();

		float[][][] rawData = DataGrabber.grab3DDataFromFile(sFile, 161, 357,
				401);

		// File Information
		int sN3 = 0;
		int eN3 = 160;

		int sN2 = 0;
		int eN2 = 356;

		int sN1 = 0;
		int eN1 = 400;
		double dt = .004;

		int n1 = eN1 - sN1 + 1;
		int n2 = eN2 - sN2 + 1;
		int n3 = eN3 - sN3 + 1;
		
		// Frequencies to analyze and number of filters
		double fmin = 8;// Hz
		double fmax = 50.0;// Hz
		int numfilters = 21;

		//MorletTransform Creation.
		Sampling time = new Sampling(n1, dt, 0);
		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);

		gf.applyThenWrite(rawData, "TPMorlet");
		System.out.println("Morlet Complete");

		File mSampFile = new File(
				"C:/Users/Chris/Data/FrequencyData/MorletSamp.dat");
		gf.writeSampling(mSampFile, gf.getLogFrequencySampling(dt));

	}
}
