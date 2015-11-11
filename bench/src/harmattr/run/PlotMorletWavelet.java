package run;

import static java.lang.Math.*;

import java.awt.Color;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PointsView;

public class PlotMorletWavelet {
	public static void main(String[] args) {
		int start = 0;
		int end = 200;
		float shift = .4f;
		int n = end - start;
		
		float[] c = new float[n];
		float[] s = new float[n];

		float[] f = { 8, 20 };
		int nf = f.length;
		for (int i1 = 0; i1 < nf; ++i1) {
			float sig = (float) (5.336f / (2 * PI * f[i1]));
			float tInt = .004f;
			Sampling samp = new Sampling(n, tInt, start);
			for (int i = 0; i < n; ++i) {
		
				c[i] = (float) (exp(-.5 * ((float) (i * tInt - shift) / sig)
						* ((float) (i * tInt - shift) / sig)) * cos(2 * PI
						* f[i1] * (i * tInt - shift)));
				s[i] = (float) (exp(-.5 * ((float) (i * tInt - shift) / sig)
						* ((float) (i * tInt - shift) / sig)) * sin(2 * PI
						* f[i1] * (i * tInt - shift)));

				
				
			}

			PointsView pvc = new PointsView(samp, c);
			pvc.setLineColor(Color.RED);
			PointsView pvs = new PointsView(samp, s);
			pvs.setLineColor(Color.BLUE);

			PlotPanel pp = new PlotPanel();
			pp.setHLabel("Time (s)");
			pp.setVLabel("Amplitude");

			pp.addTiledView(pvc);
			pp.addTiledView(pvs);

			PlotFrame pf = new PlotFrame(pp);
			pf.setFontSizeForSlide(.9, .9);
			pf.paintToPng(1000, 5, "MorletWavelet" + f[i1] + ".png");
			pf.setVisible(true);
		}

	}
}
