package test;

import attributes.SimpleAttributes;
import morletwavelettransform.MorletTransform;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.sgl.*;

public class Test {
	public static void main(String[] args) {
		int nx = 2142;
		int ny = 100;
		int nt = 300;
		float dt = (float) .004;
		float maxf = 50;
		float minf = 10;
		float f = 30;
		float tstart = 0;
		float[][] testData = new float[nx][nt];

		// sweep build
		testData = BuildTestData.build2DConstFData(nt, nx, dt, f);

		// BuildTestData.write2DData(testdata, "constSweepFDec=30");

		// impulse response
		// testdata = BuildTestData.impulseData(nt, nx, dt, 150);
		// BuildTestData.write2DData(testdata, "sweep1.binary");
		SimplePlot.asPixels(testData);
		Sampling time = new Sampling(nt, dt, 0);
		MorletTransform gf = new MorletTransform(time, 3, 30, 50);
		Sampling freq = MorletTransform.getLogFrequencySampling(dt);

		float[][][] subbands = gf.apply(testData);
		System.out.println("Morlet Complete");

		float[][][] amp = SimpleAttributes.findAmplitude(subbands);
		System.out.println("Amplitude Complete");

		SimplePlot.asPoints(amp[0][0]);

		System.out.println("test");

	}
}
