package run;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteOrder;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayFile;
import grabbers.DataGrabber;
import morletwavelettransform.MorletTransform;
import somseis.SOM2;

public class MorletFileTestWriteClassData {
	public static void main(String[] args) {

		File file = new File("C:/Users/Chris/Data/401_4_600TPData/tpsz.dat");

		String sFile = file.getPath();

		float[][][] rawData = DataGrabber.grab3DDataFromFile(sFile, 161, 357,
				401);

		File fileR = new File(
				"C:/Users/Chris/Data/ClassData/8to60Hz21F4x4som.dat");

		// File Information
		sN3 = 0;
		eN3 = 160;

		sN2 = 0;
		eN2 = 356;

		sN1 = 0;
		eN1 = 400;
		dt = .004;

		n1 = eN1 - sN1 + 1;
		n2 = eN2 - sN2 + 1;
		n3 = eN3 - sN3 + 1;

		

		// Frequencies to analyze and number of filters
		fmin = 8;// Hz
		fmax = 60.0;// Hz
		numfilters = 21;

		// SOM creation
		int numIter = 100000;
		int numAttributes = numfilters;
		int n2SOM = 4;
		int n1SOM = 4;
		SOM2 som = new SOM2(numIter, n2SOM, n1SOM, numAttributes,
				"C:/Users/Chris/workspace/Harmonic_Attributes/",
				"TP21F_4x4Morlet");

		//MorletTransform creation
		Sampling time = new Sampling(n1, dt, 0);
		MorletTransform gf = new MorletTransform(time, numfilters, fmin, fmax);

		gf.applyThenWrite(rawData, "TP21F_4x4Morlet");

		System.out.println("train");
		som.train3DFile(n3, n2, n1);
		float[][][] classData = som.classify3DDataFile();


		ArrayFile af;
		try {
			af = new ArrayFile(fileR, "rw", ByteOrder.LITTLE_ENDIAN,
					ByteOrder.BIG_ENDIAN);
			try {
				af.writeFloats(classData);
				System.out.println("ClassData Written");
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

		} catch (FileNotFoundException e) {
			System.out.println("FileNotFound");
		}
	}

	private static int numfilters;
	private static int sN1, eN1, n1, sN2, eN2, n2, sN3, eN3, n3;
	private static double dt, fmin, fmax;

}
