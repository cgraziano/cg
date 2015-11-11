/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package wwarp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.mosaic.*;//For debugging
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Stopwatch;
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes and plots the amplitude and phase spectrums of the input signal.
 * 
 * @author Chris Graziano, Colorado School of Mines
 * @version 2014.09.08
 */
public class AmpSpectrum {

  /**
   * If true, the amplitude spectrum displayed is normalized, so that the maximum amplitude is 
   * 1. The default setting is false.
   */
  public void normalizeAmplitudeSpectrum(boolean norm) {
    _norm = norm;
  }

  public void setAmpMax(float amax) {
    _maxset = true;
    _amax = amax;
  }

  public void plotAmplitudeSpectrum(Sampling st, float[] p, int itmin, int itmax, 
    boolean simplePlot, boolean forSlide, boolean forPrint, String pngDir, String title,
    double fracWidth, double fracHeight, double aspectRatio)
  {
    //Time sampling for the time window specified.
    int nt = itmax-itmin;
    double dt = st.getDelta();
    double ft = st.getValue(itmin);

    float[] subp = zerofloat(nt);
    Sampling subst = new Sampling(nt,dt,ft);
    for (int it=0; it<nt; ++it) 
      subp[it] = p[itmin+it];

    //Frequency sampling
    int nfft = FftReal.nfftSmall(4*nt);//more time sample, the finer freq. samples
    int nf = nfft/2+1;
    double df = 1.0/(nfft*dt);
    double ff = 0.0;
    Sampling fs = new Sampling(nf,df,ff);
    float[] amp = computeAmplitudeSpectrum(subst,fs,nfft,subp);

    plotSpectrum(fs,amp,title,simplePlot,forSlide,forPrint,pngDir,fracWidth,fracHeight,aspectRatio);
  }

  public void plotAmplitudeSpectrum(Sampling st, float[][] p, int itmin, int itmax, 
    boolean simplePlot, boolean forSlide, boolean forPrint, String pngDir, String title,
    double fracWidth, double fracHeight, double aspectRatio)
  {
    //Time sampling for the time window specified.
    int nx = p.length;
    int nt = itmax-itmin;
    double dt = st.getDelta();
    double ft = st.getValue(itmin);

    float[][] subp = zerofloat(nt,nx);
    Sampling subst = new Sampling(nt,dt,ft);
    for (int ix=0; ix<nx; ++ix) 
      for (int it=0; it<nt; ++it) 
        subp[ix][it] = p[ix][itmin+it];

    //Frequency sampling
    int nfft = FftReal.nfftSmall(4*nt);//more time sample, the finer freq. samples
    int nf = nfft/2+1;
    double df = 1.0/(nfft*dt);
    double ff = 0.0;
    Sampling fs = new Sampling(nf,df,ff);
    float[] amptest = computeAmplitudeSpectrum(subst,fs,nfft,subp[0]);
    float[] amp = new float[amptest.length];
    float[] subp1 = new float[nt];
    for (int ix=0; ix<nx; ++ix) 
      amp = add(amp,computeAmplitudeSpectrum(subst,fs,nfft,subp[ix]));

    plotSpectrum(fs,amp,title,simplePlot,forSlide,forPrint,pngDir,fracWidth,fracHeight,aspectRatio);
  }


//private
////////////////////////////////////////////////////////////////////////////////////////////
  private float _amax;
  private boolean _norm = false;
  private boolean _maxset= false;

  private float[] computeAmplitudeSpectrum(Sampling st, Sampling fs, int nfft, float[] p) {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nf = fs.getCount();
    double df = fs.getDelta();
    double ff = fs.getFirst();

    //Real-to-complex fast Fourier transform.
    FftReal fft = new FftReal(nfft);
    float[] cf = zerofloat(2*nf);
    copy(nt,p,cf);
    fft.realToComplex(-1,cf,cf);

    //Adjust phase for possibly non-zero time of first sample.
    float[] wft = rampfloat(0.0f,(float) (-2.0*FLT_PI*df*ft),nf);
    cf = cmul(cf,cmplx(cos(wft),sin(wft)));

    float[] af = cabs(cf);
    return af;
  }

  private void plotSpectrum(Sampling sf, float[] spec, String title, 
    boolean simplePlot, boolean forSlide, boolean forPrint, String pngDir,
    double fracWidth, double fracHeight, double aspectRatio) 
  {
    //Amplitude spectrum normalized
    System.out.println("dekunut");
    dump(spec);
    if (_norm) {
      float amax = max(max(spec),FLT_EPSILON);
      spec = mul(1.0f/amax,spec);
    }
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.LOWER_LEFT);
    if (_norm)
      sp.setVLabel("Amplitude (normalized)");
    else
      sp.setVLabel("Amplitude)");
    sp.setHLabel("Frequency (Hz)");
    sp.setVLimits(0,max(spec)+0.01);
    sp.setSize(500,170);

    if (_maxset)
      sp.setVLimits(0,_amax);
    PointsView pv = sp.addPoints(sf,spec);
    pv.setLineWidth(2.0f);
    if (simplePlot == true) {
      sp.addTitle(title);
    }
    if (forPrint == true) {
      sp.setFontSizeForPrint(10.0,469.0);
      pngDir = pngDir+title+"paper"+"onecol.png";
      sp.paintToPng(720.0,6.51,pngDir);

    }
    if (forSlide == true) {
      sp.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio);
      pngDir = pngDir+title+"w"+fracWidth+"h"+fracHeight+"slide.png";
      sp.paintToPng(720.0,3.0,pngDir);
    }
  }




}

