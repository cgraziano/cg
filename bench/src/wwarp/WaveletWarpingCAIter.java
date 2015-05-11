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
 * Estimates one wavelet from alignment by warping sequences or images.
 * The two sequences or images are assumed to have been convolved with 
 * the same wavelet. Warping of one sequence or image to align with 
 * the other will cause the wavelet to be stretched or squeezed, and this
 * distortion enables us to estimate the wavelet.
 *
 * For images, convolution with the wavelet is assumed to be in only 
 * the 1st dimension. For definiteness, this 1st dimension is assumed
 * to be time.
 *
 * @author Chris Graziano, Colorado School of Mines
 * @version 2014.09.08
 */

 public class WaveletWarpingCAIter {
   /**
   * Sets the min-max range of times used to estimate wavelet.
   * @param itmin minimum time, in samples.
   * @param itmax maximum time, in samples.
   */
  public void setTimeRange(int itmin, int itmax) {
    _itmin = itmin;
    _itmax = itmax;
  }

  /**
   * Sets the stability factor by which to scale zero-lag of correlations.
   * A factor slightly greater than one may stabilize estimates of
   * inverse wavelets A.
   * @param sfac stability factor.
   */
  public void setStabilityFactor(double sfac) {
    _sfac = sfac;
  }

  /**
   * Sets weights that determine how much each samples in f and g 
   * influences the esimtimate of the inverse wavelets. The weighting
   * system created is a boxcar.
   * @param boxVal the value of the boxcar.
   * @param itmin the first sample to have a weight of boxVal.
   * @param itmax the last sample to have a weight of boxVal.
   * @param ntf the length of the weighting vector (should be the length of f).
   */
   public void setBoxWeights(float boxVal, int itmin, int itmax, int ntf) {
     float[] w = zerofloat(ntf);
     for (int i=itmin; i<=itmax; ++i) {
       w[i] = boxVal;
     }
     _w1D = w;
     _weights = true;
   }

  /**
   * Sets weights that determine how much each samples in f and g 
   * influences the esimtimate of the inverse wavelets. The weighting
   * system created is a boxcar.
   * @param boxVal the value of the boxcar.
   * @param itmin the first time sample to have a weight of boxVal.
   * @param itmax the last time sample to have a weight of boxVal.
   * @param ixmin the first space sample to have a weight of boxVal.
   * @param ixmax the last space sample to have a weight of boxVal.
   * @param ntf the number of time samples in the weighting vector (should be nf)
   * @param nx the number of space samples in the weighting vector (should be nx)
   */
   public void setBoxWeights(float boxVal, int itmin, int itmax, int ixmin, int ixmax, 
     int ntf, int nx) {
     float[][] w = zerofloat(ntf,nx);
     for (int ix = ixmin; ix<=ixmax; ++ix) {
       for (int it=itmin; it<=itmax; ++it) { 
         w[ix][it] = boxVal;
       }
     }
     _w2D = w;
     _weights = true;
   }


  /**
   * Sets weights that determine how much each sample in f and g 
   * influences the esimtimate of the inverse wavelets. The weighting
   * system created is a gaussian.
   * @param peak the maximum height of the gaussian curve.
   * @param low the value the gaussian falls off towards.
   * @param itmin the first sample to have a weight of boxVal.
   * @param itmax the last sample to have a weight of boxVal.
   * @param ntf the length of the weighting vector (should be the length of f).
   * @param nx the number of space samples in the weighting vector (should be nx)
   */
   public void setGaussWeights(float peak, float low, int itmin, int itmax, int ntf) {
     float[] w = zerofloat(ntf);
     float center = (itmax+itmin)/2.0f;
     float width = (itmax-itmin)/2.0f;
     for (int i=0; i<ntf; ++i) {
       w[i] = peak*exp(-(i-center)*(i-center)/(2.0f*width*width))+low;
     }
     _w1D = w;
     _weights = true;
   }

  /**
   * Sets weights that determine how much each sample in f and g 
   * influences the esimtimate of the inverse wavelets. The weighting
   * system created is a 2D gaussian.
   * @param peak the maximum height of the gaussian curve.
   * @param low the value the gaussian falls off towards.
   * @param itmin the first sample to have a weight of boxVal.
   * @param itmax the last sample to have a weight of boxVal.
   * @param ixmin the first space sample to have a weight of boxVal.
   * @param ixmax the last space sample to have a weight of boxVal.
   * @param ntf the length of the weighting vector (should be the length of f).
   * @param nx the number of space samples in the weighting vector (should be nx)
   */
   public void setGaussWeights(float peak, float low, int itmin, int itmax, 
     int ixmin, int ixmax, int ntf, int nx) {
     float[][] w = zerofloat(ntf,nx);
     float centert = (itmax+itmin)/2.0f;
     float centerx = (ixmax+ixmin)/2.0f;
     float widtht = (itmax-itmin)/2.0f;
     float widthx = (ixmax-ixmin)/2.0f;
     float x = 0.0f;
     float t = 0.0f;
     for (int ix=0; ix<nx; ++ix) {
       for (int it=0; it<ntf; ++it) {
         x = (ix-centerx)*(ix-centerx)/(2.0f*widthx*widthx); 
         t = (it-centert)*(it-centert)/(2.0f*widtht*widtht); 
         w[ix][it] = peak*exp(-(x+t))+low;
       }
     }
     _w2D = w;
     _weights = true;
   }

  public float[] getWeights1D() {
    return _w1D;
  }

  public float[][] getWeights2D() {
    return _w2D;
  }

  /**
   * Sets weights that determine how much each samples in f and g 
   * influences the esimtimate of the inverse wavelets.
   * @param weights an array of weights corresponding to each sample in f and sg.
   */
   public void setWeights(float[] w) {
     _w1D = w;
     _weights = true;
   }

  /**
   * Sets weights that determine how much each samples in f and g 
   * influences the esimtimate of the inverse wavelets.
   * @param weights an array of weights corresponding to each sample in f and sg.
   */
   public void setWeights(float[][] w) {
     _w2D = w;
     _weights = true;
   }

  public float[][] getWaveletC(
    int niter, int nb, int kb, int nc, int kc, float[] u, float[] f, float[] g)
  {
    Check.argument(-nb<kb,"-nb<kb");
    Check.argument(kb<=0,"kb<=0");

    int nf = f.length;
    int ng = g.length;
    //Set default weights (equal weight to all samples)
    if (_weights==false) {
      _w1D = fillfloat(1.0f,nf);
    }
    trace("niter = "+niter);

    float sumsqdiff = 0.0f;
    float sumsqdiffpre = 0.0f;
    float[] b = new float[nb];
    float[] c = new float[nc];
    float[] b0 = new float[nb];
    float[] c0 = new float[nc];
    float[] d = new float[nf];
    float[] bg = new float[nf];
    float[] sbg = new float[nf];
    float[] csbg = new float[nf];
    _sumsqdiffpre = new float[niter];
    _sumsqdiff = new float[niter];
    Warper warp = new Warper();

    b0[-kb] = 1.0f;
    c0[-kc] = 1.0f;
    _sumsqdiffpre[0] = rmsDiff(nb,kb,b0,nc,kc,c0,u,f,g);
    if (_cset==true) {
      c = _givenc;
    }
    else 
      c = getWaveletH(nc,kc,nb,kb,b0,u,f,g);
    trace("rms diff. = "+rmsDiff(nb,kb,b0,nc,kc,c,u,f,g)+" (shaping filter)");
    _sumsqdiff[0] = rmsDiff(nb,kb,b,nc,kc,c,u,f,g);

    for (int iter=1; iter<niter; ++iter) {
      //trace("****iter = "+iter+"**********");
      b = getWaveletH(nc,kc,c,nb,kb);
      //SimplePlot.asPoints(applyA(nb,kb,b,c));
      b = normalizeSSD1(b);
      //trace("rms diff before = "+rmsDiff(nb,kb,b,nc,kc,c0,u,f,g));
      _sumsqdiffpre[iter] = rmsDiff(nb,kb,b,nc,kc,c0,u,f,g);
      c = getWaveletH(nc,kc,nb,kb,b,u,f,g);
      //c = getWaveletH(nb,kb,b,nc,kc);
      //trace("rms diff after = "+rmsDiff(nb,kb,b,nc,kc,c,u,f,g));
      //SimplePlot.asPoints(c);
      _sumsqdiff[iter] = rmsDiff(nb,kb,b,nc,kc,c,u,f,g);
      /*if (iter>niter-2) {
        SimplePlot.asPoints(b);
        SimplePlot.asPoints(c);
      }*/
      //trace("******************************");

    }
    trace("rms diff. = "+rmsDiff(nb,kb,b,nc,kc,c,u,f,g)+" (last iteration)");
    return new float[][]{c,b};
  }

  public float rmsDiff(int nb, int kb, float[] b, int nc, int kc, float[] c, 
    float[] u, float[] f, float[] g) {
    Warper warp = new Warper();
    float[] bg = applyA(nb,kb,b,g);
    float[] sbg = warp.applyS(u,bg);
    float[] csbg = applyH(nc,kc,c,sbg);
    float[] d = sub(f,csbg);
    return (float)dot(d,d);
  }

  public float[] getInverseA(int na, int ka, int nh, int kh, float[] h,
    float[] u, float[] f, float[] g)
  {
    int nt = u.length;
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");
    // Matrix P = HSLG.
    Warper warp = new Warper();
    float[][] p = new float[na][];
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[] dg = delay(lag,g);
      float[] sdg = warp.applyS(u,dg);
      p[ia] = applyH(nh,kh,h,sdg);
    }
    // Auto-correlation H'H and cross-correlation H'delta.
    float[] one = {1.0f};
    float[] ch1 = new float[na];
    float[] chh = new float[na];
    xcor(nh,kh,h,1,0,one,na,ka,ch1);
    xcor(nh,kh,h,nh,kh,h,na, 0,chh);
    // Weight applied to HA = I term.
    double wrms = _wha*rms(f);
    double w = wrms*wrms;
    // Matrix C = P'P+wH'H and vector b = P'f+wH'delta.
    DMatrix c = new DMatrix(na,na);
    DMatrix b = new DMatrix(na,1);
    for (int ia=0; ia<na; ++ia) {
      for (int ja=0; ja<na; ++ja) {
        double cij = dot(p[ia],p[ja])+w*chh[abs(ia-ja)];
        c.set(ia,ja,cij);
      }
      c.set(ia,ia,c.get(ia,ia)*_sfac);
      double bi = dot(p[ia],f)+w*ch1[ia];
      b.set(ia,0,bi);
    }
    //System.out.println("c=\n"+c);
    //System.out.println("b=\n"+b);
    // Solve for inverse filter a using Cholesky decomposition of C.
    // Normalize a such that rms(a) = 1.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    float amax = 0.0f;
    for (int ia=0; ia<na; ++ia) {
    aa[ia] = (float)a.get(ia,0);
    amax = max(amax,abs(aa[ia]));
    }
    return aa;
    //return mul(aa,1.0f/amax);
    }


  /**
   * Estimates the wavelet h from the inverse wavelet a.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   */
  public float[] getWaveletH(int na, int ka, float[] a, int nh, int kh) {
    float[] one = new float[na];
    one[-ka] = 1.0f;
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,na,ka,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);
    caa[0] *= _sfac;
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }
  public float[] getWaveletH(
    int nh, int kh, int na, int ka, float[] a,
    float[] u, float[] f, float[] g)
  {
    Warper warp = new Warper();
    int nt = u.length;
    // Sequence q = SAg.
    float[] ag = applyA(na,ka,a,g);
    float[] q = warp.applyS(u,ag);
    //SimplePlot.asPoints(ag);
    //SimplePlot.asPoints(q);
    // Autocorrelation Q'Q and crosscorrelation Q'f.
    float[] cqf = new float[nh];
    float[] cqq = new float[nh];
    int mt = (0<=_itmin && _itmin<_itmax && _itmax<nt)?1+_itmax-_itmin:nt;
    f = copy(mt,_itmin,f);
    q = copy(mt,_itmin,q);
    xcor(mt,0,q,mt,0,f,nh,kh,cqf);
    xcor(mt,0,q,mt,0,q,nh,0,cqq);
    //cqq[0] *= _sfac;
    // Solve for wavelet h.
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cqq);
    return stm.solve(cqf);
  }
/**
   * Applies the specified inverse wavelet A.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyA(int na, int ka, float[] a, float[] f) {
    return convolve(na,ka,a,f);
  }
  public float[][] applyA(int na, int ka, float[] a, float[][] f) {
    return convolve(na,ka,a,f);
  }

   /**
   * Applies the specified wavelet H.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   * @param h array of coefficients for the wavelet h.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyH(int nh, int kh, float[] h, float[] f) {
    return convolve(nh,kh,h,f);
  }
  public float[][] applyH(int nh, int kh, float[] h, float[][] f) {
    return convolve(nh,kh,h,f);
  }

  /**
   * Normalizes an array, so that the 
   * sum of square differences is 1.
   */
  public float[] normalizeSSD1(float[] x) {
    int nx = x.length;
    float sum = 0.0f;
    float[] y = new float[nx];
    for (int i=0; i<nx; ++i)
      sum += x[i]*x[i];
    sum = sqrt(sum);
    y = div(x,sum);
    return y;
  }


  public float[] getSumSqDiff() {
    return _sumsqdiff;
  }
  public float[] getSumSqDiffPre() {
    return _sumsqdiffpre;
  }

  public void setFirstC(float[] c) {
    _cset = true;
    _givenc = c;
    //SimplePlot.asPoints(_givenc);

  }
  
  
  /**
   * Returns the rms value of the image/trace.
   */
  public float rms(float[] x) {
    return (float)sqrt(dotK(x,x)/x.length);
  }
  public float rms(float[][] x) {
    return (float)sqrt(dot(x,x)/x.length/x[0].length);
  }
  /**
   * Returns the rms value of the image specified between 
   * the two time indices.
   */
  public float rms(int itmin, int itmax, float[][] x) {
    int nt = itmax-itmin+1;
    return (float)sqrt(dot(itmin,itmax,x,x)/x.length/nt);
  }
  /**
   * Returns the rms value of the PS image specified between 
   * the two time indices 228 and 700.
   */
  public float rmsPS(float[][] x) {
    int nt = 700-228+1;
    return (float)sqrt(dotPS(x,x)/x.length/nt);
  }


///////////////////////////////////////////////////////////////////////////
  // private
  private float _wha = 0.0f;
  private double _sfac = 1.0;
  private int _itmin = -1;
  private int _itmax = -1;
  private int _ng0 = 0;//staring size of array/image.
  private float[] _sumsqdiff;
  private float[] _sumsqdiffpre;
  private float[] _a;
  private float[] _w1D;
  private float[] _givenc;
  private float[][] _w2D;
  private double[][] _eigvect;
  private double[][] _eigval;
  private boolean _weights;
  private boolean _cset = false;
  private DMatrix _z;


  private double dotK(float[] x, float[] y) {
    int nt = x.length;
    double sum = 0.0;
    for (int it=0; it<nt; ++it) 
      sum += x[it]*y[it];
    return sum;
  }
  private double dot(float[] x, float[] y) {
    int nt = x.length;
    int itlo = (_itmin>0)?_itmin:0;
    int ithi = (_itmax>0)?_itmax:nt-1;
    float[] w = fillfloat(1.0f,nt);
    if (_weights) {
      w = _w1D;
    }
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) 
      sum += x[it]*y[it]*w[it];
    return sum;
  }

  private double dot(float[] x, float[] y, float[] w) {
    int nt = x.length;
    int itlo = (_itmin>0)?_itmin:0;
    int ithi = (_itmax>0)?_itmax:nt-1;
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) 
      sum += x[it]*y[it]*w[it];
    return sum;
  }

  private double dot(float[][] x, float[][] y) {
    int nx = x.length;
    int nt = x[0].length;
    float[][] w = fillfloat(1.0f,nt,nx);
    if (_weights) {
      w = _w2D;
    }
    double sum = 0.0;
    for (int ix=0; ix<nx; ++ix) 
      sum += dot(x[ix],y[ix],w[ix]);
    return sum;
  }
  private double dot(int itmin, int itmax, float[] x, float[] y) {
    int nt = x.length;
    int itlo = itmin;
    int ithi = itmax; 
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) 
      sum += x[it]*y[it];
    return sum;
  }
  private double dot(int itmin, int itmax, float[][] x, float[][] y) {
    int n = x.length;
    double sum = 0.0;
    for (int i=0; i<n; ++i) 
      sum += dot(itmin,itmax,x[i],y[i]);
    return sum;
  }

  private double dotPS(float[] x, float[] y) {
    int nt = x.length;
    int itlo = 228;//(_itmin>0)?_itmin:0;
    int ithi = 700;//(_itmax>0)?_itmax:nt-1;
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) 
      sum += x[it]*y[it];
    return sum;
  }
  private double dotPS(float[][] x, float[][] y) {
    int n = x.length;
    double sum = 0.0;
    for (int i=0; i<n; ++i) 
      sum += dotPS(x[i],y[i]);
    return sum;
  }
    
  /**
   * Given an amount of squeezing (r>1) or stretching (r<1), constructs
   * the appropriate band-pass filter.
   */
  private BandPassFilter constructBPF(float r) {
    float width = 0.20f/r;
    BandPassFilter bpf = new BandPassFilter(0.0,0.5/r,width,0.01);
    return bpf;
  }

  /**
   * Returns y(t) = x(t-lag).
   */
  private static float[] delay(int lag, float[] x) {
    int nt = x.length;
    int itlo = max(0,lag);   // 0 <= it-lag
    int ithi = min(nt,nt+lag); // it-lag < nt
    float[] y = new float[nt];
    for (int it=0; it<itlo; ++it)
      y[it] = 0.0f;
    for (int it=itlo; it<ithi; ++it)
      y[it] = x[it-lag];
    for (int it=ithi; it<nt; ++it)
      y[it] = 0.0f;
    return y;
  }
  private static float[][] delay(int lag, float[][] x) {
    int n = x.length;
    float[][] y = new float[n][];
    for (int i=0; i<n; ++i)
      y[i] = delay(lag,x[i]);
    return y;
  }

  

  /**
   * Returns y(t) = h(t)*x(t), where * denotes convolution.
   */
  private static float[] convolve(int nh, int kh, float[] h, float[] x) {
    int nt = x.length;
    float[] y = new float[nt];
    convolve(nh,kh,h,x,y);
    return y;
  }
  private static void convolve(
    int nh, int kh, float[] h, float[] f,  float[] g)
  {
    int nt = f.length;
    conv(nh,kh,h,nt,0,f,nt,0,g);
  }
  private static float[][] convolve(int nh, int kh, float[] h, float[][] x) {
    int n = x.length;
    int nt = x[0].length;
    float[][] y = new float[n][nt];
    for (int i=0; i<n; ++i)
      convolve(nh,kh,h,x[i],y[i]);
    return y;
  }

  private double[] convertFToD(float[] f) {
    int nf = f.length;
    double[] d = new double[nf];
    for (int i=0; i<nf; ++i)
      d[i] = f[i];
    return d;
  }

  private float[][] convertDToF(double[][] d) {
    int nd2 = d.length;
    int nd1 = d[0].length;
    float[][] f = new float[nd2][nd1];
    for (int i2=0; i2<nd2; ++i2)
      for (int i1=0; i1<nd1; ++i1)
        f[i2][i1] = (float) d[i2][i1];
    return f;
  }

  private void plotAmplitudeSpectrum(Sampling st, float[] f, 
    int itmin, int itmax, String title) {
    //Time sampling for the specified time window.
    int nt = itmax-itmin;
    double dt = st.getDelta();
    double ft = st.getValue(itmin);
    float[] subf = zerofloat(nt);
    Sampling subst = new Sampling(nt,dt,ft);
    for (int i=0; i<nt; ++i) 
      subf[i] = f[itmin+i];

    //Frequency sampling
    int nfft = FftReal.nfftSmall(4*nt);//more time sample, the finer freq. samples
    int nf = nfft/2+1;
    double df = 1.0/(nfft*dt);
    double ff = 0.0;
    Sampling sf = new Sampling(nf,df,ff);
    float[] amp = computeAmplitudeSpectrum(subst,sf,nfft,subf);
    plotSpectrum(sf,amp,title);
  }

  private float[] computeAmplitudeSpectrum(Sampling st, Sampling sf, int nfft, float[] f) {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nf = sf.getCount();
    double df = sf.getDelta();
    double ff = sf.getFirst();

    //Real-to-complex fast Fourier transform.
    FftReal  fft = new FftReal(nfft);
    float[] cf = zerofloat(2*nf);
    copy(nt,f,cf);
    fft.realToComplex(-1,cf,cf);

    //Adjust phase for possibly non-zero time of first sample.
    float[] wft = rampfloat(0.0f,(float)(-2.0f*FLT_PI*df*ft),nf);
    cf = cmul(cf,cmplx(cos(wft),sin(wft)));

    float[] af = cabs(cf);
    //Amplitude spectrum normalized
    //float amax = max(max(af),FLT_EPSILON);
    //af = mul(1.0f/amax,af);
    return af;
  }

  private void plotSpectrum(Sampling sf,float[] f,String title) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.LOWER_LEFT);
    sp.setVLabel("Amplitude");
    sp.setHLabel("Frequency (Hz)");
    sp.setSize(750,400);
    sp.addTitle(title);
    //sp.setVLimits(0.0,1.0);
    PointsView pv = sp.addPoints(sf,f);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
 }



