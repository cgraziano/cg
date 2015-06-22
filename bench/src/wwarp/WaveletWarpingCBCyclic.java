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

public class WaveletWarpingCBCyclic{
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
   * The 2D and 3D methods can only be parallized.
   */
  public void setParallel(boolean parallel) {
    Parallel.setParallel(parallel);
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
    checkArguments(nb,kb,nc,kc);
    initializeAllRequiredFields(f,niter,nb,nc);
    setCAndB(cGuess,bGuess);

    for (int iter=0; iter<niter; ++iter) 
    {
      Stopwatch sw = new Stopwatch();
      sw.start();
      trace("*****************************************");
      trace("iteration = "+iter);
      computeAndStoreAllInitialResiduals(nc,kc,_c,nb,kb,_b,u,f,g);
      computeAndStoreInitialMeasuresOfResiduals(_dataResInit1D,_bPenaltyResInit,_cPenaltyResInit);
      _b = getInverseB(nb,kb,nc,kc,_c,u,f,g);
      _c = getWaveletC(nc,kc,nb,kb,_b,_sfac,u,f,g);
      computeAndStoreAllFinalResiduals(nc,kc,_c,nb,kb,_b,u,f,g);
      computeAndStoreFinalMeasuresOfResiduals(_dataResFina1D,_bPenaltyResInit,_cPenaltyResFina);
      calcAndStoreAllDiagnosticMeasures(iter);
      printBAndC();
      printIterationCompletionTime(sw);
      _lastIter = iter;
      trace("*****************************************");
      if (isRmsPercentChangeSmallerThanMinimumPercentChange())
      {
        return new float[][]{_c,_b};
      }
    }
    return new float[][]{_c,_b};
  }
    public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bGuess, int nc, int kc, float[] cGuess, 
    float[][] u, float[][] f, float[][] g, int niter)
  {
    checkArguments(nb,kb,nc,kc);
    initializeAllRequiredFields(f,niter,nb,nc);
    setCAndB(cGuess,bGuess);

    for (int iter=0; iter<niter; ++iter) 
    {
      Stopwatch sw = new Stopwatch();
      sw.start();
      trace("*****************************************");
      trace("iteration = "+iter);
      computeAndStoreAllInitialResiduals(nc,kc,_c,nb,kb,_b,u,f,g);
      computeAndStoreInitialMeasuresOfResiduals(_dataResInit2D,_bPenaltyResInit,_cPenaltyResInit);
      _b = getInverseB(nb,kb,nc,kc,_c,u,f,g);
      _c = getWaveletC(nc,kc,nb,kb,_b,_sfac,u,f,g);
      trace("b");
      dump(_b);
      trace("c");
      dump(_c);
      computeAndStoreAllFinalResiduals(nc,kc,_c,nb,kb,_b,u,f,g);
      computeAndStoreFinalMeasuresOfResiduals(_dataResFina2D,_bPenaltyResInit,_cPenaltyResFina);
      calcAndStoreAllDiagnosticMeasures(iter);
      printBAndC();
      printIterationCompletionTime(sw);
      _lastIter = iter;
      trace("*****************************************");
      if (isRmsPercentChangeSmallerThanMinimumPercentChange())
      {
        return new float[][]{_c,_b};
      }
    }
    return new float[][]{_c,_b};
  }
  public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bGuess, int nc, int kc, float[] cGuess, 
    float[][][] u, float[][][] f, float[][][] g, int niter)
  {
    checkArguments(nb,kb,nc,kc);
    initializeAllRequiredFields(f,niter,nb,nc);
    setCAndB(cGuess,bGuess);

    for (int iter=0; iter<niter; ++iter) 
    {
      Stopwatch sw = new Stopwatch();
      sw.start();
      Stopwatch swb = new Stopwatch();
      Stopwatch swc = new Stopwatch();
      trace("*****************************************");
      trace("iteration = "+iter);
      computeAndStoreAllInitialResiduals(nc,kc,_c,nb,kb,_b,u,f,g);
      computeAndStoreInitialMeasuresOfResiduals(_dataResInit3D,_bPenaltyResInit,_cPenaltyResInit);
      swb.start();
      _b = getInverseB(nb,kb,nc,kc,_c,u,f,g);
      swb.stop();
      trace("b time = "+swb.time());
      swc.start();
      _c = getWaveletC(nc,kc,nb,kb,_b,_sfac,u,f,g);
      trace("c time = "+swc.time());
      trace("b");
      dump(_b);
      trace("c");
      dump(_c);
      computeAndStoreAllFinalResiduals(nc,kc,_c,nb,kb,_b,u,f,g);
      computeAndStoreFinalMeasuresOfResiduals(_dataResFina3D,_bPenaltyResInit,_cPenaltyResFina);
      calcAndStoreAllDiagnosticMeasures(iter);
      printBAndC();
      printIterationCompletionTime(sw);
      _lastIter = iter;
      trace("*****************************************");
      if (isRmsPercentChangeSmallerThanMinimumPercentChange())
      {
        return new float[][]{_c,_b};
      }
    }
    return new float[][]{_c,_b};
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

    // Matrix P = CSG.
    float[][] p = new float[nb][];
    Warper warp = new Warper();
    for (int ib=0,lag=kb; ib<nb; ++ib,++lag) {
      p[ib] = applyC(nc,kc,c,warp.applyS(u,delay(lag,g)));
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
    return mul(bb,1.0f/bmax);
  }
  public float[] getInverseB(
    final int nb, final int kb, final int nc, final int kc, final float[] c, 
    final float[][] u, final float[][] f, final float[][] g)
  {
    int nx = u.length;
    int nt = u[0].length;

    // Matrix P = CSG.
    final float[][][] p = new float[nb][][];
    Parallel.loop(nb,new Parallel.LoopInt() 
    {
      public void compute(int ib) 
      {
        int lag = kb+ib;
        Warper warp = new Warper();
        p[ib] = applyC(nc,kc,c,warp.applyS(u,delay(lag,g)));
      }
    });

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

    /*
    // Matrix PP = P'P and vector pf = P'f.
    final double[][] ppTemp = new double[nb][nb];
    final double[] pfTemp = new double[nb];
    Parallel.loop(nb,new Parallel.LoopInt() 
    {
      public void compute(int ib) 
      {
        Warper warp = new Warper();
        int lagi = kb+ib;
        for (int jb=0; jb<nb; ++jb) {
          int lagj = kb+jb;
          double ppij = dot(applyC(nc,kc,c,warp.applyS(u,delay(lagi,g))),applyC(nc,kc,c,warp.applyS(u,delay(lagj,g))));
          ppTemp[ib][jb] = ppij;
        }
        ppTemp[ib][ib] = ppTemp[ib][ib]*(1.0+_sfac);
        double pfi = dot(applyC(nc,kc,c,warp.applyS(u,delay(lagi,g))),f);
        pfTemp[ib] = pfi;
      }
    });
    DMatrix pp = new DMatrix(ppTemp);
    DMatrix pf = new DMatrix(nb,1,pfTemp);
    */

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
    return mul(bb,1.0f/bmax);
  }
  public float[] getInverseB(
    final int nb, final int kb, final int nc, final int kc, final float[] c, 
    final float[][][] u, final float[][][] f, final float[][][] g)
  {
    int nx = u.length;
    int nt = u[0].length;

    // Matrix PP = P'P and vector pf = P'f.
    final double[][] ppTemp = new double[nb][nb];
    final double[] pfTemp = new double[nb];
    Parallel.loop(nb,new Parallel.LoopInt() 
    {
      public void compute(int ib) 
      {
        Warper warp = new Warper();
        int lagi = kb+ib;
        float[][][] pi = applyC(nc,kc,c,warp.applyS(u,delay(lagi,g)));
        for (int jb=0; jb<nb; ++jb) {
          int lagj = kb+jb;
          double ppij = dot(pi,applyC(nc,kc,c,warp.applyS(u,delay(lagj,g))));
          ppTemp[ib][jb] = ppij;
        }
        ppTemp[ib][ib] = ppTemp[ib][ib]*(1.0+_sfac);
        double pfi = dot(pi,f);
        pfTemp[ib] = pfi;
      }
    });
    DMatrix pp = new DMatrix(ppTemp);
    DMatrix pf = new DMatrix(nb,1,pfTemp);

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
    return mul(bb,1.0f/bmax);
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
    float[] q = warp.applyS(u,bg);

    // Autocorrelation Q'Q and crosscorrelation Q'f.
    float[] cqf = new float[nc];
    float[] cqq = new float[nc];
    int mt = (0<=_itmin && _itmin<_itmax && _itmax<nt)?1+_itmax-_itmin:nt;
    f = copy(mt,_itmin,f);
    q = copy(mt,_itmin,q);
    xcor(mt,0,q,mt,0,f,nc,kc,cqf);
    xcor(mt,0,q,mt,0,q,nc, 0,cqq);
    float cqfMax = max(cqf);
    float cqqMax = max(cqq);
    float max = max(cqqMax,cqfMax);
    cqq = div(cqq,max);
    cqf = div(cqf,max);
    trace("cqqSum");
    dump(cqq);
    trace("cqfSum");
    dump(cqf);


    // Solve for wavelet C.
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cqq);
    return stm.solve(cqf);
  }
  public float[] getWaveletC(
    int nc, int kc, int nb, int kb, float[] b, double stabFact,
    float[][] u, float[][] f, float[][] g)
  {
    int nx = u.length;
    int nt = u[0].length;
    Warper warp = new Warper();

    // Sequence q = SBg.
    float[][] bg = applyC(nb,kb,b,g);
    float[][] q = warp.applyS(u,bg);

    // Autocorrelation Q'Q and crosscorrelation Q'f.
    float[] cqf = new float[nc];
    float[] cqq = new float[nc];
    float[] cqfSum = new float[nc];
    float[] cqqSum = new float[nc];
    int mt = (0<=_itmin && _itmin<_itmax && _itmax<nt)?1+_itmax-_itmin:nt;
    f = copy(mt,nx,_itmin,0,f);
    q = copy(mt,nx,_itmin,0,q);
    for (int ix=0; ix<nx; ++ix) {
      xcor(mt,0,q[ix],mt,0,f[ix],nc,kc,cqf);
      xcor(mt,0,q[ix],mt,0,q[ix],nc, 0,cqq);
      cqfSum = add(cqfSum,cqf);
      cqqSum = add(cqqSum,cqq);
    }
    float cqfMax = max(cqfSum);
    float cqqMax = max(cqqSum);
    float max = max(cqqMax,cqfMax);
    cqqSum = div(cqqSum,max);
    cqfSum = div(cqfSum,max);

    // Solve for wavelet C.
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cqqSum);
    return stm.solve(cqfSum);
  }
  public float[] getWaveletC(
    int nc, int kc, int nb, int kb, float[] b, double stabFact,
    float[][][] u, float[][][] f, float[][][] g)
  {
    int nx3 = u.length;
    int nx2 = u[0].length;
    int nt = u[0][0].length;
    Warper warp = new Warper();

    // Sequence q = SBg.
    float[][][] bg = applyC(nb,kb,b,g);
    float[][][] q = warp.applyS(u,bg);

    // Autocorrelation Q'Q and crosscorrelation Q'f.
    float[] cqf = new float[nc];
    float[] cqq = new float[nc];
    float[] cqfSum = new float[nc];
    float[] cqqSum = new float[nc];
    int mt = (0<=_itmin && _itmin<_itmax && _itmax<nt)?1+_itmax-_itmin:nt;
    f = copy(mt,nx2,nx3,_itmin,0,0,f);
    q = copy(mt,nx2,nx3,_itmin,0,0,q);
    for (int ix3=0; ix3<nx3; ++ix3) {
      for (int ix2=0; ix2<nx2; ++ix2) {
        xcor(mt,0,q[ix3][ix2],mt,0,f[ix3][ix2],nc,kc,cqf);
        xcor(mt,0,q[ix3][ix2],mt,0,q[ix3][ix2],nc, 0,cqq);
        cqfSum = add(cqfSum,cqf);
        cqqSum = add(cqqSum,cqq);
      }
    }
    float cqfMax = max(cqfSum);
    float cqqMax = max(cqqSum);
    float max = max(cqqMax,cqfMax);
    cqqSum = div(cqqSum,max);
    cqfSum = div(cqfSum,max);

    // Solve for wavelet C.
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cqqSum);
    return stm.solve(cqfSum);
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
  public float[][][] applyC(int nh, int kh, float[] h, float[][][] f) {
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
  public float rmsOfObjectiveFunction(
    float[][][] dataResiduals, float[] bPenaltyResiduals, float[] cPenaltyResiduals) 
  {
    int npenb = bPenaltyResiduals.length;
    int npenc = cPenaltyResiduals.length;
    int nt = _itmax-_itmin+1;
    int nx3 = dataResiduals.length;
    int nx2 = dataResiduals[0].length;
    int n = nx3*nx2*nt+npenb+npenc;
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

   public float[] getAllResRmsS() {
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
  private float _allResRmsInit;
  private float _allResRmsFina;
  private float _minRmsPercentChange;
  private float _rmsPercentChange;
  private float[] _rmsPercentChangeS;
  private float[] _allResRmsAllS;
  private float[] _c;
  private float[] _b;
  private float[] _bPenaltyResInit;
  private float[] _bPenaltyResFina;
  private float[] _cPenaltyResInit;
  private float[] _cPenaltyResFina;
  private float[] _dataResInit1D;
  private float[] _dataResFina1D;
  private float[][] _dataResInit2D;
  private float[][] _dataResFina2D;
  private float[][][] _dataResInit3D;
  private float[][][] _dataResFina3D;


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
  private double dot(float[][][] x, float[][][] y) {
    int nx3 = x.length;
    int nx2 = x[0].length;
    int nt = x[0][0].length;
    double sum = 0.0;
    for (int ix3=0; ix3<nx3; ++ix3) {
      for (int ix2=0; ix2<nx2; ++ix2) {
        sum += dot(x[ix3][ix2],y[ix3][ix2]);
      }
    }
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
  private double dotAll(float[][][] x, float[][][] y) {
    int nx3 = x.length;
    int nx2 = x[0].length;
    int nt = x[0][0].length;
    double sum = 0.0;
    for (int ix3=0; ix3<nx3; ++ix3) {
      for (int ix2=0; ix2<nx2; ++ix2) {
        sum += dotAll(x[ix3][ix2],y[ix3][ix2]);
      }
    }
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
  private double dot(int itmin, int itmax, float[][][] x, float[][][] y) {
    int nx3 = x.length;
    int nx2 = x[0].length;
    int nt = x[0][0].length;
    double sum = 0.0;
    for (int ix3=0; ix3<nx3; ++ix3) {
      for (int ix2=0; ix2<nx2; ++ix2) {
        sum += dot(itmin,itmax,x[ix3][ix2],y[ix3][ix2]);
      }
    }
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
  private static float[][][] delay(int lag, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    float[][][] y = new float[n3][n2][];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        y[i3][i2] = delay(lag,x[i3][i2]);
      }
    }
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
  private static float[][][] convolve(int nh, int kh, float[] h, float[][][] x) {
    int n3 = x.length;
    int n2 = x[0].length;
    int nt = x[0][0].length;
    float[][][] y = new float[n3][n2][nt];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        convolve(nh,kh,h,x[i3][i2],y[i3][i2]);
      }
    }
    return y;
  }
  private float[] computeDataResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
   float[] u, float[] f, float[] g) 
  {
    Warper warp = new Warper();
    return sub(applyC(nc,kc,c,warp.applyS(u,applyC(nb,kb,b,g))),f);
  }
  private float[][] computeDataResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[][] u, float[][] f, float[][] g) 
  {
    Warper warp = new Warper();
    return sub(applyC(nc,kc,c,warp.applyS(u,applyC(nb,kb,b,g))),f);
  }
  private float[][][] computeDataResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[][][] u, float[][][] f, float[][][] g) 
  {
    Warper warp = new Warper();
    return sub(applyC(nc,kc,c,warp.applyS(u,applyC(nb,kb,b,g))),f);
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

  private boolean isRmsPercentChangeSmallerThanMinimumPercentChange() {
    return abs(_rmsPercentChange)<_minRmsPercentChange;
  }

  private void calcAndStoreAllDiagnosticMeasures(int iter) {
    trace("_allResRmsInit = "+_allResRmsInit);
    trace("_allResRmsFina = "+_allResRmsFina);
    _allResRmsAllS[iter] = _allResRmsInit;
    _allResRmsAllS[iter+1] = _allResRmsFina;
    if (iter==0) 
    {
      _rmsPercentChangeS[0] = 0.0f;//RMS percent change only exists with each iteration.
                                   //Set this to zero because initially there is not a RMS percent change, but
                                   //there is a starting RMS value. This starting value is in _allResRmsAllS[0].
                                   //The overall goal is to make _rmsPercentChangeS and _allResRmsAllS have the
                                   //same length for plotting purposes.
      _rmsPercentChange = percentChange(_allResRmsAllS[iter+1],_allResRmsAllS[iter]);
      _rmsPercentChangeS[iter+1] = _rmsPercentChange;
    }
    else
    {
      _rmsPercentChange = percentChange(_allResRmsAllS[iter+1],_allResRmsAllS[iter]);
      _rmsPercentChangeS[iter+1] = _rmsPercentChange;
    }
  }


  private void computeAndStoreInitialMeasuresOfResiduals(float[] dataResInit, float[] bPenaltyResInit, float[] cPenaltyResInit) 
  {
    _allResRmsInit = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
    trace("allResRmsInit = "+_allResRmsInit);
  }
  private void computeAndStoreInitialMeasuresOfResiduals(float[][] dataResInit, float[] bPenaltyResInit, float[] cPenaltyResInit) 
  {
    _allResRmsInit = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
    trace("allResRmsInit = "+_allResRmsInit);
  }
  private void computeAndStoreInitialMeasuresOfResiduals(float[][][] dataResInit, float[] bPenaltyResInit, float[] cPenaltyResInit) 
  {
    _allResRmsInit = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
    trace("allResRmsInit = "+_allResRmsInit);
  }

  private void computeAndStoreFinalMeasuresOfResiduals(float[] dataResFina, float[] bPenaltyResFina, float[] cPenaltyResFina) 
  {
    _allResRmsFina = rmsOfObjectiveFunction(dataResFina,bPenaltyResFina,cPenaltyResFina);
    trace("allResRmsFina= "+_allResRmsFina);
  }
  private void computeAndStoreFinalMeasuresOfResiduals(float[][] dataResFina, float[] bPenaltyResFina, float[] cPenaltyResFina) 
  {
    _allResRmsFina = rmsOfObjectiveFunction(dataResFina,bPenaltyResFina,cPenaltyResFina);
    trace("allResRmsFina = "+_allResRmsFina);
  }
  private void computeAndStoreFinalMeasuresOfResiduals(float[][][] dataResFina, float[] bPenaltyResFina, float[] cPenaltyResFina) 
  {
    _allResRmsFina = rmsOfObjectiveFunction(dataResFina,bPenaltyResFina,cPenaltyResFina);
    trace("allResRmsFina = "+_allResRmsFina);
  }

  private void checkArguments(int nb, int kb, int nc, int kc) {
    Check.argument(-nb<kb,"-nb<kb");
    Check.argument(kb<=0,"kb<=0");
    Check.argument(-nc<kc,"-nc<kc");
    Check.argument(kc<=0,"kc<=0");
  }

  private void initializeAllRequiredFields(float[] f, int niter, int nb, int nc) {
    int nt = f.length;
    _allResRmsAllS = new float[niter+1];//we use a length of niter+1 because we want to include the starting 
                                        //RMS value in this array, which would mean we have niter+1 measurements.
    _rmsPercentChangeS = new float[niter+1];//we use a length of niter+1 because we want _allResRmsAllS and 
                                            //_rmsPercentChangeS to have the same length.
    _b = new float[nb];
    _c = new float[nc];
    _bPenaltyResInit = new float[nb];
    _cPenaltyResInit = new float[nc];
    _bPenaltyResFina = new float[nb];
    _cPenaltyResFina = new float[nc];
    _dataResInit1D = new float[nt];
    _dataResFina1D = new float[nt];
  }
  private void initializeAllRequiredFields(float[][] f, int niter, int nb, int nc) {
    int nx2 = f.length;
    int nt = f[0].length;
    _allResRmsAllS = new float[niter+1];
    _rmsPercentChangeS = new float[niter+1];
    _b = new float[nb];
    _c = new float[nc];
    _bPenaltyResInit = new float[nb];
    _cPenaltyResInit = new float[nc];
    _bPenaltyResFina = new float[nb];
    _cPenaltyResFina = new float[nc];
    _dataResInit2D = new float[nx2][nt];
    _dataResFina2D = new float[nx2][nt];
  }
  private void initializeAllRequiredFields(float[][][] f, int niter, int nb, int nc) {
    int nx3 = f.length;
    int nx2 = f[0].length;
    int nt = f[0][0].length;
    _allResRmsAllS = new float[niter+1];
    _rmsPercentChangeS = new float[niter+1];
    _b = new float[nb];
    _c = new float[nc];
    _bPenaltyResInit = new float[nb];
    _cPenaltyResInit = new float[nc];
    _bPenaltyResFina = new float[nb];
    _cPenaltyResFina = new float[nc];
    _dataResInit3D = new float[nx3][nx2][nt];
    _dataResFina3D = new float[nx3][nx2][nt];
  }


  private void setCAndB(float[] c, float[] b){
    _c = copy(c);
    _b = copy(b);
  }

  private void computeAndStoreAllInitialResiduals(
      int nc, int kc, float[] c, int nb, int kb, float[] b, 
      float[] u, float[] f, float[] g)
  {
    _dataResInit1D = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
  }
  private void computeAndStoreAllInitialResiduals(
      int nc, int kc, float[] c, int nb, int kb, float[] b, 
      float[][] u, float[][] f, float[][] g)
  {
    _dataResInit2D = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
  }
  private void computeAndStoreAllInitialResiduals(
      int nc, int kc, float[] c, int nb, int kb, float[] b, 
      float[][][] u, float[][][] f, float[][][] g)
  {
    _dataResInit3D = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
  }

  private void computeAndStoreAllFinalResiduals(
      int nc, int kc, float[] c, int nb, int kb, float[] b, 
      float[] u, float[] f, float[] g)
  {
    _dataResFina1D = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
  }
  private void computeAndStoreAllFinalResiduals(
      int nc, int kc, float[] c, int nb, int kb, float[] b, 
      float[][] u, float[][] f, float[][] g)
  {
    _dataResFina2D = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
  }
  private void computeAndStoreAllFinalResiduals(
      int nc, int kc, float[] c, int nb, int kb, float[] b, 
      float[][][] u, float[][][] f, float[][][] g)
  {
    _dataResFina3D = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
  }

  private void printIterationCompletionTime(Stopwatch sw) {
    sw.stop();
    trace("Iteration Completion Time: "+sw.time());
  }
  private void printBAndC() {
    trace("b:");
    dump(_b);
    trace("c:");
    dump(_c);
  }




  
      
  private static void trace(String s) {
    System.out.println(s);
  }


    
}

