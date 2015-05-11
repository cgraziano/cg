/****************************************************************************
Copyright (c) 2014, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package wwarp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates a wavelet from alignment by warping of sequences or images.
 * The two sequences or images are assumed to have been convolved with the
 * same wavelet. Warping of one sequence or image to align with the other will
 * cause the wavelet to be stretched or squeezed, and this distortion enables
 * us to estimate the wavelet.
 * <p>
 * For images, convolution with the wavelet is assumed to be in only the 1st
 * dimension. For definiteness, this 1st dimension is assumed to be time in
 * the documentation below.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2014.01.27
 */
public class WaveletWarpingCA {

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
   * Sets the min-max range of frequencies used to estimate the wavelet.
   * If the specified min-max bounds on frequency are not a subset of the
   * zero-Nyquist range [0,0.5], then no bandpass filter is used. The default
   * is to use no bandpass filter.
   * @param fmin minimum frequency, in cycles/sample.
   * @param fmax maximum frequency, in cycles/sample.
   */
  public void setFrequencyRange(double fmin, double fmax) {
    if (fmin<fmax && (0.0<fmin || fmax<0.5)) {
      _bpf = new BandPassFilter(fmin,fmax,0.05,0.01);
    } else {
      _bpf = null;
    }
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
   * Returns inverse wavelet a estimated by warping one sequence to another.
   * The sequences are related by warping such that f[t] ~ g[u[t]].
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param u array of samples for warping u[t].
   * @param f array of samples for sequence f[t].
   * @param g array of samples for sequence g[t]
   * @return array of coefficients for the inverse wavelet a.
   */
  public float[] getInverseA(
    int na, int ka, float[] u, float[] f, float[] g)
  {
    int nt = u.length;
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    // Differences d for all lags of inverse wavelet a.
    float[][] d = computeDifferences(na,ka,u,f,g);

    // The matrix C and right-hand-side vector b, for Ca = b. For zero lag, we
    // have a0 = a[-ka] = 1, so that only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags.
    int ma = na-1;
    DMatrix c = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jc=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double cij = dot(d[ia],d[ja]);
        c.set(ic,jc,cij);
        ++jc;
      }
      c.set(ic,ic,c.get(ic,ic)*(1.0+_sfac));
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,1.0*bi);
      ++ic;
    }
    //System.out.println("c=\n"+c);
    //System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of C.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // lag 0, so a0 = 1
      } else {
        aa[ia] = (float)a.get(ic,0);
        ++ic;
      }
    }
    return aa;
  }

  public float[] getInverseA(
    int na, int ka, float[][] u, float[][] f, float[][] g)
  {
    int nt = u.length;
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    // Differences d for all lags of inverse wavelet a.
    float[][][] d = computeDifferences(na,ka,u,f,g);

    // The matrix C and right-hand-side vector b, for Ca = b. For zero lag, we
    // have a0 = a[-ka] = 1, so that only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags.
    int ma = na-1;
    DMatrix c = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jc=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double cij = dot(d[ia],d[ja]);
        c.set(ic,jc,cij);
        ++jc;
      }
      c.set(ic,ic,c.get(ic,ic)*(1.0+_sfac));
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    //System.out.println("c=\n"+c);
    //System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of C.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // lag 0, so a0 = 1
      } else {
        aa[ia] = (float)a.get(ic,0);
        ++ic;
      }
    }
    return aa;
  }

  /**
   * Estimates the wavelet h from the inverse wavelet a.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   */
  public float[] getWaveletC(int na, int ka, float[] a, int nh, int kh) {
    float[] one = {1.0f};
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }

  /**
   * Applies the specified wavelet H.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   * @param h array of coefficients for the wavelet h.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyC(int nh, int kh, float[] h, float[] f) {
    return convolve(nh,kh,h,f);
  }
  public float[][] applyC(int nh, int kh, float[] h, float[][] f) {
    return convolve(nh,kh,h,f);
  }

  /**
   * Applies the bandpass filter B, if any was specified.
   * If no bandpass filter has been specified, then this method simply returns
   * a copy of the specified input sequence.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyB(float[] f) {
    float[] g = new float[f.length];
    if (_bpf!=null) {
      _bpf.apply(f,g);
    } else {
      copy(f,g);
    }
    return g;
  }
  public float[][] applyB(float[][] f) {
    int n = f.length;
    float[][] g = new float[n][];
    for (int i=0; i<n; ++i)
      g[i] = applyB(f[i]);
    return g;
  }

  /**
   * Estimates that shaping filter that will shape SBg to f.
   * @param nc number of samples in wavelet c.
   * @param kc the sample index for a[0].
   * @param nb number of samples in the wavelet h.
   * @param kb the sample index for h[0].
   * @param b array of coefficients for the inverse wavelet a.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   */
  public float[] getWaveletC(
    int nc, int kc, int nb, int kb, float[] b, float stabFact,
    float[] u, float[] f, float[] g)
  {
    int nt = u.length;
    Warper warp = new Warper();

    // Sequence q = SBg.
    float[] bg = applyC(nb,kb,b,g);
    float[] q0 = warp.applyS(u,bg);

    //Q'Q
    DMatrix qq = new DMatrix(nc,nc);
    for (int ic=0,lagi=kc; ic<nc; ++ic,++lagi) {
      float[] qi = delay(lagi,q0);
      for (int jc=0,lagj=kc; jc<nc; ++jc,++lagj) {
        float[] qj = delay(lagj,q0);
        double qiqj = dot(qi,qj);
        qq.set(ic,jc,qiqj);
        if (ic==jc)
          qq.set(ic,jc,((1.0+(double) stabFact)*qq.get(ic,jc)));
      }
    }
    //trace(qq.toString());
    //Q'f
    DMatrix qf = new DMatrix(nc,1);
    for (int ic=0,lagi=kc; ic<nc; ++ic,++lagi) {
      float[] qi = delay(lagi,q0);
      double qif = dot(qi,f);
      qf.set(ic,0,qif);
    }

    // Solve for wavelet C.
    DMatrixChd chd = new DMatrixChd(qq);
    DMatrix h = chd.solve(qf);
    return convertDToF(h.getArray());
  }
  public float[] getWaveletC(
    int nc, int kc, int nb, int kb, float[] b, float stabFact,
    float[][] u, float[][] f, float[][] g)
  {
    int nt = u[0].length;
    Warper warp = new Warper();

    // Sequence q = SBg.
    float[][] bg = applyC(nb,kb,b,g);
    float[][] q0 = warp.applyS(u,bg);

    //Q'Q
    DMatrix qq = new DMatrix(nc,nc);
    for (int ic=0,lagi=kc; ic<nc; ++ic,++lagi) {
      float[][] qi = delay(lagi,q0);
      for (int jc=0,lagj=kc; jc<nc; ++jc,++lagj) {
        float[][] qj = delay(lagj,q0);
        double qiqj = dot(qi,qj);
        qq.set(ic,jc,qiqj);
        if (ic==jc)
          qq.set(ic,jc,((1.0+(double) stabFact)*qq.get(ic,jc)));
      }
    }
    
    //Q'f
    DMatrix qf = new DMatrix(nc,1);
    for (int ic=0,lagi=kc; ic<nc; ++ic,++lagi) {
      float[][] qi = delay(lagi,q0);
      double qif = dot(qi,f);
      qf.set(ic,0,qif);
    }

    // Solve for wavelet C.
    DMatrixChd chd = new DMatrixChd(qq);
    DMatrix h = chd.solve(qf);
    trace("h:");
    trace(h.toString());
    return convertDToF(h.getArray());
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private double _sfac = 0.0;
  private int _itmin = -1;
  private int _itmax = -1;
  private BandPassFilter _bpf;

  /**
   * Returns the array of differences D = B(SLG-F).
   */
  private float[][] computeDifferences(
    int na, int ka, float[] u, float[] f, float[] g)
  {
    int ntf = u.length;
    float[][] d = new float[na][];
    float[] dfi = new float[ntf];
    float[] sdgi = new float[ntf];
    float[] di = new float[ntf];
    Warper warp = new Warper();
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      dfi = delay(lagi,f);
      sdgi = warp.applyS(u,delay(lagi,g));
      di = sub(sdgi,dfi);
      d[ia] = applyB(di);
    }
    return d;
  }
  private float[][][] computeDifferences(
    int na, int ka, float[][] u, float[][] f, float[][] g)
  {
    int ntf = u[0].length;
    int nx = u.length;
    float[][][] d = new float[na][][];
    float[][] dfi = new float[nx][ntf];
    float[][] sdgi = new float[nx][ntf];
    float[][] di = new float[nx][ntf];
    Warper warp = new Warper();
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      dfi = delay(lagi,f);
      sdgi = warp.applyS(u,delay(lagi,g));
      di = sub(sdgi,dfi);
      d[ia] = applyB(di);
    }
    return d;
  }

  private double dot(float[] x, float[] y) {
    int nt = x.length;
    int itlo = (_itmin>0)?_itmin:0;
    int ithi = (_itmax>0)?_itmax:nt-1;
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) {
      sum += x[it]*y[it];
    }
    return sum;
  }

  private double dot(float[][] x, float[][] y) {
    int nx = x.length;
    int nt = x[0].length;
    double sum = 0.0;
    for (int ix=0; ix<nx; ++ix) 
      sum += dot(x[ix],y[ix]);
    return sum;
  }
  private double dotAll(float[] x, float[] y) {
    int nt = x.length;
    double sum = 0.0;
    for (int it=0; it<nt; ++it) {
      sum += x[it]*y[it];
    }
    return sum;
  }
  private double dotAll(float[][] x, float[][] y) {
    int nx = x.length;
    int nt = x[0].length;
    double sum = 0.0;
    for (int ix=0; ix<nx; ++ix) 
      sum += dotAll(x[ix],y[ix]);
    return sum;
  }
  private double dot(int itmin, int itmax, float[] x, float[] y) {
    int nt = x.length;
    int itlo = (itmin>0)?itmin:0;
    int ithi = (itmax>0)?itmax:nt-1;
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) {
      sum += x[it]*y[it];
    }
    return sum;
  }

  private double dot(int itmin, int itmax, float[][] x, float[][] y) {
    int nx = x.length;
    int nt = x[0].length;
    double sum = 0.0;
    for (int ix=0; ix<nx; ++ix) 
      sum += dot(itmin,itmax,x[ix],y[ix]);
    return sum;
  }

  /**
   * Returns y(t) = x(t-lag).
   */
  public static float[] delay(int lag, float[] x) {
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

  private float[] convertDToF(double[] d) {
    int nf = d.length;
    float[] f = new float[nf];
    for (int i=0; i<nf; ++i)
      f[i] = (float) d[i];
    return f;
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
  private static void trace(String s) {
    System.out.println(s);
  }


}













