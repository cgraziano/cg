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
 * different wavelets. Warping one sequence or image to align with 
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

public class WaveletWarpingCBCyclic {
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
   * When the percent change of the RMS of the residuals is below this value,
   * the iterations stop.
   */
  public void setMinPercentChange(float minRmsPercentChange) {
    _minRmsPercentChange = minRmsPercentChange;
  }
  

  /**
   * Sets the stability factor by which to scale zero-lag of correlations.
   * A factor slightly greater than one may stabilize estimates of
   * wavelets c and inverse wavelets b.
   * @param sfac stability factor.
   */
  public void setStabilityFactor(double sfac) {
    _sfac = sfac;
  }

  /**
   * Estimates the wavelet c and the inverse wavelet b using a cyclic search.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index for b[0].
   * @param bguess array of coefficients for the first guess of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index for c[0].
   * @param cguess array of coefficients for the first guess of the wavelet c.
   * @param niter number of iterations the Gauss-Newton method will have.
   * @param u array of samples for warping u[t].
   * @param f array of samples for sequence f[t].
   * @param g array of samples for sequence g[t].
   * @param niter the maximum number of iterations.
   * @return array of coefficients for the inverse wavelet a.   
   */
  public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bGuess, int nc, int kc, float[] cGuess, 
    float[] u, float[] f, float[] g, int niter)
  {
    float[] c = copy(cGuess);
    trace("cGuess");
    dump(cGuess);
    float[] b = copy(bGuess);
    SimplePlot.asPoints(cGuess);
    _allResRmsAllS = new float[niter];
    float[] dataResInit = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
    float[] bPenaltyResInit = new float[nb];
    float[] cPenaltyResInit = new float[nc];
    _allResRmsAllS[0] = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
    trace("allResRmsAllS0 = "+_allResRmsAllS[0]);

    for (int iter=1; iter<niter; ++iter) {
      trace("iteration = "+iter);
      b = getInverseB(nb,kb,nc,kc,c,u,f,g);
      trace("in between rms residuals with b and previous c");
      trace("rms = "+rmsOfObjectiveFunction(computeDataResidual(nc,kc,c,nb,kb,b,u,f,g),bPenaltyResInit,cPenaltyResInit));
      c = getWaveletC(nc,kc,nb,kb,b,_sfac,u,f,g);
      if (iter%5==0) {
        SimplePlot sp = new SimplePlot();
        sp.addTitle("b iter = "+iter);
        sp.addPoints(b);
        SimplePlot sp1 = new SimplePlot();
        sp1.addTitle("c iter = "+iter);
        sp1.addPoints(c);
      }
      
      dataResInit = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
      _allResRmsAllS[iter] = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
      trace("allResRmsAllSiter = "+_allResRmsAllS[iter]);
      float rmsPercentChange = percentChange(_allResRmsAllS[iter],_allResRmsAllS[iter-1]);
      trace("rmsPercentChange = "+rmsPercentChange);
      _lastIter = iter;
      if (abs(rmsPercentChange)<_minRmsPercentChange) {
        return new float[][]{c,b};
      }
    }
    return new float[][]{c,b};
  }


  /**
   * Estimates the wavelet c and the inverse wavelet b using a cyclic search.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index for b[0].
   * @param bguess array of coefficients for the first guess of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index for c[0].
   * @param cguess array of coefficients for the first guess of the wavelet c.
   * @param niter number of iterations the Gauss-Newton method will have.
   * @param u array of samples for warping u[t].
   * @param f array of samples for sequence f[t].
   * @param g array of samples for sequence g[t].
   * @param niter the maximum number of iterations.
   * @return array of coefficients for the inverse wavelet a.   
   */
  public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bGuess, int nc, int kc, float[] cGuess, 
    float[][] u, float[][] f, float[][] g, int niter)
  {
    float[] c = copy(cGuess);
    trace("cGuess");
    dump(cGuess);
    float[] b = copy(bGuess);
    SimplePlot.asPoints(cGuess);
    _allResRmsAllS = new float[niter];
    float[][] dataResInit = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
    float[] bPenaltyResInit = new float[nb];
    float[] cPenaltyResInit = new float[nc];
    _allResRmsAllS[0] = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
    trace("allResRmsAllS0 = "+_allResRmsAllS[0]);

    for (int iter=1; iter<niter; ++iter) {
      trace("iteration = "+iter);
      b = getInverseB(nb,kb,nc,kc,c,u,f,g);
      trace("in between rms residuals with b and previous c");
      trace("rms = "+rmsOfObjectiveFunction(computeDataResidual(nc,kc,c,nb,kb,b,u,f,g),bPenaltyResInit,cPenaltyResInit));
      c = getWaveletC(nc,kc,nb,kb,b,_sfac,u,f,g);
      SimplePlot sp = new SimplePlot();
      sp.addTitle("b iter = "+iter);
      sp.addPoints(b);
      SimplePlot sp1 = new SimplePlot();
      sp1.addTitle("c iter = "+iter);
      sp1.addPoints(c);
      dataResInit = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
      _allResRmsAllS[iter] = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
      trace("allResRmsAllSiter = "+_allResRmsAllS[iter]);
      float rmsPercentChange = percentChange(_allResRmsAllS[iter],_allResRmsAllS[iter-1]);
      trace("rmsPercentChange = "+rmsPercentChange);
      _lastIter = iter;
      if (abs(rmsPercentChange)<_minRmsPercentChange) {
        return new float[][]{c,b};
      }
    }
    return new float[][]{c,b};
  }
 
  /**
   * Given an estimate of c, can solve for the inverse wavelet b.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index of b[0].
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index of c[0].
   * @param c estimate of the wavelet c.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   */
  public float[] getInverseB(
    int nb, int kb, int nc, int kc, float[] c, 
    float[] u, float[] f, float[] g)
  {
    int nt = u.length;
    Check.argument(-nb<kb,"-nb<kb");
    Check.argument(kb<=0,"kb<=0");

    // Matrix P = CSG.
    float[][] p = new float[nb][];
    Warper warp = new Warper();
    for (int ib=0,lag=kb; ib<nb; ++ib,++lag) {
      float[] dgi = delay(lag,g);
      float[] sdgi = warp.applyS(u,dgi);
      p[ib] = applyC(nc,kc,c,sdgi);
    }

    // Matrix PP = P'P and vector pf = P'f.
    DMatrix pp = new DMatrix(nb,nb);
    DMatrix pf = new DMatrix(nb,1);
    for (int ib=0; ib<nb; ++ib) {
      for (int jb=0; jb<nb; ++jb) {
        double ppij = dot(p[ib],p[jb]);
        pp.set(ib,jb,ppij);
      }
      pp.set(ib,ib,pp.get(ib,ib)*(1.0+_sfac));
      double pfi = dot(p[ib],f);
      pf.set(ib,0,pfi);
    }
    //System.out.println("c=\n"+c);
    //System.out.println("b=\n"+b);

    // Solve for inverse filter b using Cholesky decomposition of PP.
    // Normalize a such that rms(a) = 1.
    DMatrixChd chd = new DMatrixChd(pp);
    DMatrix b = chd.solve(pf);
    float[] bb = new float[nb];
    float bmax = 0.0f;
    for (int ib=0; ib<nb; ++ib) {
      bb[ib] = (float)b.get(ib,0);
      bmax = max(bmax,abs(bb[ib]));
    }
    return bb;
    //return mul(bb,1.0f/bmax);
  }
  public float[] getInverseB(
    int nb, int kb, int nc, int kc, float[] c, 
    float[][] u, float[][] f, float[][] g)
  {
    int nt = u[0].length;
    Check.argument(-nb<kb,"-nb<kb");
    Check.argument(kb<=0,"kb<=0");

    // Matrix P = CSG.
    float[][][] p = new float[nb][][];
    Warper warp = new Warper();
    for (int ib=0,lag=kb; ib<nb; ++ib,++lag) {
      float[][] dgi = delay(lag,g);
      float[][] sdgi = warp.applyS(u,dgi);
      p[ib] = applyC(nc,kc,c,sdgi);
    }

    // Matrix PP = P'P and vector pf = P'f.
    DMatrix pp = new DMatrix(nb,nb);
    DMatrix pf = new DMatrix(nb,1);
    for (int ib=0; ib<nb; ++ib) {
      for (int jb=0; jb<nb; ++jb) {
        double ppij = dot(p[ib],p[jb]);
        pp.set(ib,jb,ppij);
      }
      pp.set(ib,ib,pp.get(ib,ib)*(1.0+_sfac));
      double pfi = dot(p[ib],f);
      pf.set(ib,0,pfi);
    }
    //trace(pp.toString());
    //trace(pf.toString());
    //System.out.println("c=\n"+c);
    //System.out.println("b=\n"+b);

    // Solve for inverse filter b using Cholesky decomposition of PP.
    // Normalize a such that rms(a) = 1.
    DMatrixChd chd = new DMatrixChd(pp);
    DMatrix b = chd.solve(pf);
    float[] bb = new float[nb];
    float bmax = 0.0f;
    for (int ib=0; ib<nb; ++ib) {
      bb[ib] = (float)b.get(ib,0);
      bmax = max(bmax,abs(bb[ib]));
    }
    return bb;
    //return mul(bb,1.0f/bmax);
  }


  /**
   * Estimates the wavelet c from the inverse wavelet a.
   * @param na number of samples in the inverse wavelet a.
   * @param ka the sample index for a[0].
   * @param a array of coefficients for the inverse wavelet a.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index for c[0].
   */
  public float[] getWaveletC(int na, int ka, float[] a, int nc, int kc) {
    float[] one = {1.0f};
    float[] ca1 = new float[nc];
    float[] caa = new float[nc];
    xcor(na,ka,a,1,0,one,nc,kc,ca1);
    xcor(na,ka,a,na,ka,a,nc, 0,caa);
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
  }

  /**
   * Estimates that shaping filter that will shape SBg to f.
   * @param nc number of samples in wavelet c.
   * @param kc the sample index for a[0].
   * @param nb number of samples in the wavelet b.
   * @param kb the sample index for b[0].
   * @param b array of coefficients for the inverse wavelet a.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   */
  public float[] getWaveletC(
    int nc, int kc, int nb, int kb, float[] b, double stabFact,
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
          qq.set(ic,jc,((1.0+stabFact)*qq.get(ic,jc)));
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
    int nc, int kc, int nb, int kb, float[] b, double stabFact,
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
          qq.set(ic,jc,((1.0+stabFact)*qq.get(ic,jc)));
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

  /*public float[] getWaveletCOld(
    int nh, int kh, int na, int ka, float[] a,
    float[][] u, float[][] f, float[][] g)
  {
    int nt = u[0].length;
    int nx = u.length;


    // Sequence q = SLAg.
    Warper warp = new Warper();
    float[][] ag = applyC(na,ka,a,g);
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

    // Solve for wavelet h.
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cqq);
    trace("new");
    trace("cqq");
    dump(cqq);
    trace("cqf");
    dump(cqf);
    trace("h:");
    dump(stm.solve(cqf));
    return stm.solve(cqf);
  }
  */



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
  public float rmsAll(float[] x) {
    int nt = x.length;
    return (float)sqrt(dotAll(x,x)/nt);
  }
  public float rmsAll(float[][] x) {
    int nt = x[0].length;
    return (float)sqrt(dotAll(x,x)/(nt*x.length));
  }
  public float rms(int itmin, int itmax, float[] x) {
    int nt = itmax-itmin+1;
    return (float)sqrt(dot(itmin,itmax,x,x)/nt);
  }
  public float rms(int itmin, int itmax, float[][] x) {
    int nt = itmax-itmin+1;
    return (float)sqrt(dot(itmin,itmax,x,x)/(nt*x.length));
  }

  public float rmsOfObjectiveFunction(
    float[] dataResiduals, float[] bPenaltyResiduals, float[] cPenaltyResiduals) 
  {
    int npenb = bPenaltyResiduals.length;
    int npenc = cPenaltyResiduals.length;
    int nt = _itmax-_itmin+1;
    int n = nt+npenb+npenc;
    double data2NormSq = dot(dataResiduals,dataResiduals);
    double bPenalty2NormSq = dotAll(bPenaltyResiduals,bPenaltyResiduals);
    double cPenalty2NormSq = dotAll(cPenaltyResiduals,cPenaltyResiduals);
    return (float) sqrt((data2NormSq+bPenalty2NormSq+cPenalty2NormSq)/n);
  }

  public float rmsOfObjectiveFunction(
    float[][] dataResiduals, float[] bPenaltyResiduals, float[] cPenaltyResiduals) 
  {
    int npenb = bPenaltyResiduals.length;
    int npenc = cPenaltyResiduals.length;
    int nt = _itmax-_itmin+1;
    int nx = dataResiduals.length;
    int n = nx*nt+npenb+npenc;
    double data2NormSq = dot(dataResiduals,dataResiduals);
    double bPenalty2NormSq = dotAll(bPenaltyResiduals,bPenaltyResiduals);
    double cPenalty2NormSq = dotAll(cPenaltyResiduals,cPenaltyResiduals);
    return (float) sqrt((data2NormSq+bPenalty2NormSq+cPenalty2NormSq)/n);
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

   public float[] getAllResRmsAllS() {
    return _allResRmsAllS;
  }

  public int getLastIter() {
    return _lastIter;
  }




  //Private
  private double _sfac = 0.0;
  private int _itmin = -1;
  private int _itmax = -1;
  private int _lastIter = 0;
  private float[] _allResRmsAllS;
  private float _minRmsPercentChange;


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
  private float[] computeDataResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
   float[] u, float[] f, float[] g) 
  {
    Warper warp = new Warper();
    float[] csbg = applyC(nc,kc,c,warp.applyS(u,applyC(nb,kb,b,g)));
    return sub(csbg,f);
  }
  private float[][] computeDataResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[][] u, float[][] f, float[][] g) 
  {
    Warper warp = new Warper();
    float[][] csbg = applyC(nc,kc,c,warp.applyS(u,applyC(nb,kb,b,g)));
    return sub(csbg,f);
  }
  private float percentChange(float xf, float xi) {
    return (xf-xi)/xi*100.0f;
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

