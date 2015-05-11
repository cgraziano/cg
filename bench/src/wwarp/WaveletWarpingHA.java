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
 * @version 2014.02.19
 */
public class WaveletWarpingHA {

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
   * Sets the relative weight of the HA = I term. If this weight is zero, then
   * H is effectively decoupled from A. As this weight is increased, A will be
   * forced to be the inverse of H, such that HA = I. The default weight is
   * 0.1.
   */
  public void setWeightHA(double wha) {
    _wha = wha;
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

  public void setMaxPercentChange(float minRmsPercentChange) {
    _minRmsPercentChange = minRmsPercentChange;
  }

  public float[] getAllResRmsAllS() {
    return _allResRmsAllS;
  }

  public int getLastIter() {
    return _lastIter;
  }

  public float[][] getWaveletHInverseA(
    int na, int ka, float[] aGuess, 
    int nh, int kh, float[] hGuess,
    float[][] u, float[][] f, float[][] g, int niter) 
  {
    _allResRmsAllS = new float[niter];
    float allResRmsInit = 0.0f;
    float allResRmsFina = 0.0f;
    float rmsPercentChange = 0.0f;
    float[] a = copy(aGuess);
    float[] h = copy(hGuess);
    for (int iter=0; iter<niter; ++iter) {
      allResRmsInit = rms(computeDataResidual(nh,kh,h,na,ka,a,u,f,g)); 
      
      //Solve for a and h
      a = getInverseA(na,ka,nh,kh,h,u,f,g);
      h = getWaveletH(nh,kh,na,ka,a,u,f,g);

      allResRmsFina = rms(computeDataResidual(nh,kh,h,na,ka,a,u,f,g)); 
      rmsPercentChange = percentChange(allResRmsFina,allResRmsInit);
      _lastIter = iter;

      if (iter == 0)
        _allResRmsAllS[0] = allResRmsInit;
      else
        _allResRmsAllS[iter] = allResRmsInit;
      _allResRmsAllS[iter] = allResRmsInit;
      _allResRmsAllS[iter+1] = allResRmsFina;
      trace("iter = "+iter);
      trace("allResRmsInit = "+allResRmsInit);
      trace("allResRmsFina = "+allResRmsFina);
      trace("rmsPercentChange = "+rmsPercentChange);
      if (abs(rmsPercentChange)<_minRmsPercentChange) {
        return new float[][]{h,a};
      }

    }
    return new float[][]{h,a};
  }


  /**
   * Returns inverse wavelet a estimated by warping one sequence to another.
   * The sequences are related by warping such that f[t] ~ g[u[t]].
   * The specified wavelet h is just the current best estimate of h.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param nh number of samples in the inverse wavelet a.
   * @param kh the sample index for a[0].
   * @param h array of coefficients for the wavelet h.
   * @param u array of samples for warping u[t].
   * @param f array of samples for sequence f[t].
   * @param g array of samples for sequence g[t]
   * @return array of coefficients for the inverse wavelet a.
   */
  public float[] getInverseA(
    int na, int ka, int nh, int kh, float[] h, 
    float[] u, float[] f, float[] g)
  {
    int nt = u.length;
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    // Matrix P = HSLG.
    float[][] p = new float[na][];
    Warper warp = new Warper();
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
  public float[] getInverseA(
    int na, int ka, int nh, int kh, float[] h, 
    float[][] u, float[][] f, float[][] g)
  {
    int nt = u.length;
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    // Matrix P = HSLG.
    Warper warp = new Warper();
    float[][][] p = new float[na][][];
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[][] dg = delay(lag,g);
      float[][] sdg = warp.applyS(u,dg);
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
    // Normalize a such that max(abs(a)) = 1.
    System.out.println(c.toString());
    System.out.println(b.toString());
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
    float[] one = {1.0f};
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);
    caa[0] *= _sfac;
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }
  public float[] getWaveletH(
    int nh, int kh, int na, int ka, float[] a,
    float[] u, float[] f, float[] g)
  {
    int nt = u.length;

    // Sequence q = SLAg.
    Warper warp = new Warper();
    float[] ag = applyA(na,ka,a,g);
    float[] q = warp.applyS(u,ag);

    // Autocorrelation Q'Q and crosscorrelation Q'f.
    float[] cqf = new float[nh];
    float[] cqq = new float[nh];
    int mt = (0<=_itmin && _itmin<_itmax && _itmax<nt)?1+_itmax-_itmin:nt;
    f = copy(mt,_itmin,f);
    q = copy(mt,_itmin,q);
    xcor(mt,0,q,mt,0,f,nh,kh,cqf);
    xcor(mt,0,q,mt,0,q,nh, 0,cqq);

    // Autocorrelation A'A and crosscorrelation A'delta.
    float[] one = {1.0f};
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);

    // Weight applied to HA = I term.
    float wrms = (float)(_wha*rms(f));
    float w = wrms*wrms;

    // Matrix Q'Q+wA'A and vector Q'f+wA'delta.
    for (int ih=0; ih<nh; ++ih) {
      cqq[ih] += w*caa[ih];
      cqf[ih] += w*ca1[ih];
    }
    cqq[0] *= _sfac;

    // Solve for wavelet h.
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cqq);
    return stm.solve(cqf);
  }
  public float[] getWaveletH(
    int nh, int kh, int na, int ka, float[] a,
    float[][] u, float[][] f, float[][] g)
  {
    int nt = u[0].length;
    int nx = u.length;

    // Autocorrelation A'A and crosscorrelation A'delta.
    float[] one = {1.0f};
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);

    // Weight applied to HA = I term.
    float wrms = (float)(_wha*rms(f));
    float w = wrms*wrms;

    // Sequence q = SLAg.
    Warper warp = new Warper();
    float[][] ag = applyA(na,ka,a,g);
    float[][] q = warp.applyS(u,ag);

    // Autocorrelation Q'Q and crosscorrelation Q'f.
    float[] cqf = new float[nh];
    float[] cqq = new float[nh];
    float[] tqf = new float[nh];
    float[] tqq = new float[nh];
    int mt = (0<_itmin && _itmin<_itmax && _itmax<nt)?1+_itmax-_itmin:nt;
    f = copy(mt,nx,_itmin,0,f);
    q = copy(mt,nx,_itmin,0,q);
    for (int ix=0; ix<nx; ++ix) {
      xcor(mt,0,q[ix],mt,0,f[ix],nh,kh,tqf);
      xcor(mt,0,q[ix],mt,0,q[ix],nh, 0,tqq);
      for (int ih=0; ih<nh; ++ih) {
        cqf[ih] += tqf[ih];
        cqq[ih] += tqq[ih];
      }
    }

    // Matrix Q'Q+wA'A and vector Q'f+wA'delta.
    for (int ih=0; ih<nh; ++ih) {
      cqq[ih] += w*caa[ih];
      cqf[ih] += w*ca1[ih];
    }
    cqq[0] *= _sfac;

    // Solve for wavelet h.
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cqq);
    return stm.solve(cqf);
  }
  public float[] getWaveletHShifts(
    int nh, int kh, int na, int ka, float[] a,
    float[] u, float[] f, float[] g)
  {
    int nt = u.length;

    // Sequence q = SAg.
    float[] ag = applyA(na,ka,a,g);
    float[] q = applyS(u,ag);

    // Autocorrelation Q'Q and crosscorrelation Q'f.
    float[] cqf = new float[nh];
    float[] cqq = new float[nh];
    int mt = (0<=_itmin && _itmin<_itmax && _itmax<nt)?1+_itmax-_itmin:nt;
    f = copy(mt,_itmin,f);
    q = copy(mt,_itmin,q);
    xcor(mt,0,q,mt,0,f,nh,kh,cqf);
    xcor(mt,0,q,mt,0,q,nh, 0,cqq);

    // Autocorrelation A'A and crosscorrelation A'delta.
    float[] one = {1.0f};
    float[] ca1 = new float[nh];
    float[] caa = new float[nh];
    xcor(na,ka,a,1,0,one,nh,kh,ca1);
    xcor(na,ka,a,na,ka,a,nh, 0,caa);

    // Weight applied to HA = I term.
    float wrms = (float)(_wha*rms(f));
    float w = wrms*wrms;

    // Matrix Q'Q+wA'A and vector Q'f+wA'delta.
    for (int ih=0; ih<nh; ++ih) {
      cqq[ih] += w*caa[ih];
      cqf[ih] += w*ca1[ih];
    }
    cqq[0] *= _sfac;

    // Solve for wavelet h.
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cqq);
    return stm.solve(cqf);
  }

  /**
   * Returns the rms value of the image/trace.
   */
  public float rms(float[] x) {
    int nt = _itmax-_itmin+1;
    return (float)sqrt(dot(x,x)/nt);
  }
  public float rms(float[][] x) {
    int nt = _itmax-_itmin+1;
    return (float)sqrt(dot(x,x)/(nt*x.length));
  }
  public float rms(int itmin, int itmax, float[] x) {
    int nt = itmax-itmin+1;
    return (float)sqrt(dot(itmin,itmax,x,x)/nt);
  }
  public float rms(int itmin, int itmax, float[][] x) {
    int nt = itmax-itmin+1;
    return (float)sqrt(dot(itmin,itmax,x,x)/(nt*x.length));
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
   * Applies the low-pass anti-alias filter L.
   * If the specified warping includes squeezing, then this method attenuates
   * high frequencies that could be aliased during warping.
   * @param u array of warping times u(t).
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyL(float[] u, float[] f) {
    return aaf(RMAX,u,f);
  }
  public float[][] applyL(float[][] u, float[][] f) {
    return aaf(RMAX,u,f);
  }

  /**
   * Applies the warping operator S.
   * Does not apply an anti-alias low-pass filter.
   * @param u array of warping times u(t).
   * @param f array with input sequence f(t).
   * @return array with warped output sequence.
   */
  public float[] applyS(float[] u, float[] f) {
    return warp(u,f);
  }
  public float[][] applyS(float[][] u, float[][] f) {
    return warp(u,f);
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
    float[] bg = applyH(nb,kb,b,g);
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
    float[][] bg = applyH(nb,kb,b,g);
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
    return convertDToF(h.getArray());
  }


  /**
   * Applies the composite linear operator HSLA.
   * The sequence of operations is (1) convolution with the inverse wavelet a,
   * (2) anti-alias filtering (if squeezing), (3) warping, and (4) convolution
   * with the wavelet h.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param nh number of samples in the wavelet h.
   * @param kh the sample index for h[0].
   * @param h array of coefficients for the wavelet h.
   * @param u array[nt] of warping times u(t).
   * @param f array[nt] with input sequence.
   * @return array[nt] with output sequence.
   */
  public float[] applyHSLA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    float[] u, float[] f) 
  {
    int nt = f.length;
    float[] af = applyA(na,ka,a,f);
    float[] laf = applyL(u,af);
    float[] saf = applyS(u,laf);
    float[] hsaf = applyH(nh,kh,h,saf);
    return hsaf;
  }
  public float[][] applyHSLA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    float[][] u, float[][] f) 
  {
    int nt = f.length;
    float[][] af = applyA(na,ka,a,f);
    float[][] laf = applyL(u,af);
    float[][] saf = applyS(u,laf);
    float[][] hsaf = applyH(nh,kh,h,saf);
    return hsaf;
  }
  public float[] applyHSA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    float[] u, float[] f) 
  {
    int nt = f.length;
    float[] af = applyA(na,ka,a,f);
    float[] saf = applyS(u,af);
    float[] hsaf = applyH(nh,kh,h,saf);
    return hsaf;
  }
  /**
   * Displays the contents of P
   */
  public float[][] displayPShifts(int na, int ka, int nh, int kh, float[] h, 
    float[] u, float[] g) {
    // Matrix P = HSG.
    float[][] p = new float[na][];
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[] dg = delay(lag,g);
      float[] sdg = applyS(u,dg);
      p[ia] = applyH(nh,kh,h,sdg);
      p[ia] = delay(-lag,p[ia]);
    }
    return p;
  }

  /**
   * Makes the rms equal to 1 within a specified time range.
   */
   public float[] makeRms1(float[] x) {
     float rmsx = rms(x);
     float[] rms1x = div(x,rmsx);
     //Check
     trace("rms after normalization is "+rms(rms1x));
     return rms1x;
   }
   public float[][] makeRms1(float[][] x) {
     float rmsx = rms(x);
     float[][] rms1x = div(x,rmsx);
     //Check
     trace("rms after normalization is "+rms(rms1x));
     return rms1x;
   }
   public float[] makeRms1(int itmin, int itmax, float[] x) {
     float rmsx = rms(itmin,itmax,x);
     float[] rms1x = div(x,rmsx);
     //Check
     trace("rms after normalization is "+rms(itmin,itmax,rms1x));
     return rms1x;
   }
   public float[][] makeRms1(int itmin, int itmax, float[][] x) {
     float rmsx = rms(itmin,itmax,x);
     float[][] rms1x = div(x,rmsx);
     //Check
     trace("rms after normalization is "+rms(itmin,itmax,rms1x));
     return rms1x;
   }



  ///////////////////////////////////////////////////////////////////////////
  // private

  private static final float RMAX = 10.0f; // limits anti-alias filter
  private static final SincInterp _si = 
    SincInterp.fromErrorAndFrequency(0.01,0.40);

  private double _wha = 0.1;
  private double _sfac = 1.0;
  private int _itmin = -1;
  private int _itmax = -1;
  private float[] _allResRmsAllS;
  private float _minRmsPercentChange = 0.000f;
  private int _lastIter = 0;
  

  private double dot(float[] x, float[] y) {
    int nt = x.length;
    int itlo = (_itmin>0)?_itmin:0;
    int ithi = (_itmax>0)?_itmax:nt-1;
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) 
      sum += x[it]*y[it];
    return sum;
  }
  private double dot(float[][] x, float[][] y) {
    int n = x.length;
    double sum = 0.0;
    for (int i=0; i<n; ++i) 
      sum += dot(x[i],y[i]);
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
   * Returns the largest squeezing r(t) = u'(t) not greater than rmax.
   * If less than or equal to one, then no squeezing is implied by u(t).
   */
  private float squeezing(float rmax, float[] u) {
    int nt = u.length;
    int itlo = max(1,_itmin);
    int ithi = min(_itmax,nt-1);
    float r = 0.0f;
    for (int it=itlo; it<=ithi; ++it) {
      float du = u[it]-u[it-1];
      if (r<du)
        r = du;
    }
    return min(r,rmax);
  }
  private float squeezing(float rmax, float[][] u) {
    int n = u.length;
    float r = 0.0f;
    for (int i=0; i<n; ++i)
      r = max(r,squeezing(rmax,u[i]));
    return r;
  }

  /**
   * If necessary, applies an anti-alias filter to the sequence x(t).
   * An anti-alias filter is necessary if the warping includes squeezing.
   */
  private float[] aaf(float rmax, float[] u, float[] x) {
    int nt = x.length;
    float r = squeezing(RMAX,u);
    if (r>1.0) {
      float[] y = new float[nt];
      BandPassFilter aaf = new BandPassFilter(0.0,0.5/r,0.10/r,0.01);
      aaf.apply(x,y);
      return y;
    } else {
      return copy(x);
    }
  }
  private float[][] aaf(float rmax, float[][] u, float[][] x) {
    float r = squeezing(RMAX,u);
    if (r>1.0) {
      int nx = x.length;
      int nt = x[0].length;
      float[][] y = new float[nx][nt];
      BandPassFilter aaf = new BandPassFilter(0.0,0.5/r,0.10/r,0.01);
      for (int ix=0; ix<nx; ++ix)
        aaf.apply(x[ix],y[ix]);
      return y;
    } else {
      return copy(x);
    }
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
   * Returns y(t) = x(u(t)).
   */
  private static float[] warp(float[] u, float[] x) {
    int nt = u.length;
    float[] y = new float[nt];
    _si.interpolate(x.length,1.0,0.0,x,nt,u,y);
    y[0] *= u[1]-u[0];
    for (int it=1; it<nt; ++it)
      y[it] *= u[it]-u[it-1];
    return y;
  }
  private static float[][] warp(float[][] u, float[][] x) {
    int n = u.length;
    float[][] y = new float[n][];
    for (int i=0; i<n; ++i)
      y[i] = warp(u[i],x[i]);
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

  private float percentChange(float xf, float xi) {
    return (xf-xi)/xi*100.0f;
  }

  private float[][] computeDataResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[][] u, float[][] f, float[][] g) 
  {
    Warper warp = new Warper();
    float[][] csbg = applyH(nc,kc,c,warp.applyS(u,applyH(nb,kb,b,g)));
    return sub(csbg,f);
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

