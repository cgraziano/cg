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
import edu.mines.jtk.opt.BrentMinFinder;
import java.awt.Color;
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates two wavelets from alignment by warping sequences or images.
 * The two sequences or images are assumed to have been convolved with 
 * different wavelets. Warping of one sequence or image to align with 
 * the other will cause the wavelet to be stretched or squeezed, and this
 * distortion enables this algorithm to estimate the wavelets.
 *
 * For images, convolution with the wavelet is assumed to be in only 
 * the 1st dimension. For definiteness, this 1st dimension is assumed
 * to be time.
 *
 * @author Chris Graziano, Colorado School of Mines
 * @version 2013.03.04
 */
 public class Start {
  /**
   * Sets the min-max range of times used to estimate wavelets.
   * @param itmin minimum time, in samples.
   * @param itmax maximum time, in samples.
   */
  public void setTimeRange(int itmin, int itmax) {
    _itmin = itmin;
    _itmax = itmax;
  }

  /**
   * Enables parallel processing.
   */
  public void setParallel(boolean parallel) {
    _parallel = parallel;
  }

  /**
   * The second part of the Gauss-Newton method is a line search.
   * This method defines the lower bound of values to be tested in the line
   * search. 
   * This lower bound also serves as the tolerance value for the line search 
   * algorithm. A value will be returned by the line search algorithm when
   * the distance between the current value of the line search and the true minimum is 
   * smaller than the tolerance value.
   * The upper bound of the values to be searched by the line search is set to be 1.0.
   * The lower bound and tolerance is 0.0001.
   */
  public void setLineSearchMinScale(float minsteplength) {
    _minsteplength = minsteplength;
  }

  /**
   * Sets the largest percent change of the rms of the residuals before the
   * Gauss-Newton method will stop iterating. The default percent change required before 
   * the Gauss-Newton method will stop iterating is 0.1 percent.
   * Note that you should NOT include the negative in the percent change, we are assuming 
   * that you wish to have a decrease in the residuals.
   * @param maxrmspercentchange The maximum percent change needed to stop this iterative method.
   */
  public void setMaxPercentChange(float maxrmspercentchange) {
    _maxrmspercentchange = maxrmspercentchange;
  }

  public void setPenalize(float alpha, float sfac,
    boolean penb, boolean penc, boolean simpleStab)
  {
    _penb = penb;
    _penc = penc;
    _alpha = alpha;
    _sfac = sfac;
    _simpleStab = simpleStab;
  }



  /**
   * Sets the stability factor by which to multiply the diagonal of the matrix X by.
   * A factor slightly greater than one may stabilize estimates of
   * inverse wavelet b and wavelet c.
   * @param sfac stability factor.
   */
  public void setStabilityFactor(double sfac) {
    _sfac = sfac;
  }

  /**
   * Estimates the wavelet c and the inverse wavelet b using the Gauss-Newton algorithm.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index for b[0].
   * @param bguess array of coefficients for the first guess of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index for c[0].
   * @param cguess array of coefficients for the first guess of the wavelet c.
   * @param niter number of iterations the Gauss-Newton method will have.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   */
  public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bguess, int nc, int kc, float[] cguess,
    float[] u, float[] f, float[] g, int niter)
  {
    //trace("bguess");//db
    //dump(bguess);//db
    //trace("cguess");//db
    //dump(cguess);//db
    _steplength = new float[niter];//db
    _rmsri = new float[niter];//db
    _rmsrf = new float[niter];//db
    _twonormgrad = new float[niter];//db
    _twonormdeltac = new float[niter];//db
    _twonormdeltab = new float[niter];//db
    _condnum = new float[niter];//db
    _twonormb = new float[niter];//db
    _twonormc = new float[niter];//db
    _rmspercentchange = new float[niter];//db

    int nbc = nc+nb;
    float[] b = copy(bguess);
    float[] c = copy(cguess);
    float[] bnew = new float[nb];
    float[] cnew = new float[nc];
    float[] deltabdeltac = new float[nbc];
    float[] deltab = new float[nb];
    float[] deltac = new float[nc];
    float[] ri = new float[f.length];
    float[] rf = new float[f.length];
    float rmsrf = 0.0f;
    float rmsri = 0.0f;
    float steplength = 0.0f;
    float pc = 0.0f;
    float rmspercentchange = 0.0f;
    DMatrix[] xv = new DMatrix[2];
    DMatrix v = new DMatrix(nbc,1);
    DMatrix x = new DMatrix(nbc,nbc);
    DMatrix deltabc = new DMatrix(nbc,1);
    //Gauss-Newton iterations
    for (int iter=0; iter<niter; ++iter) {
      trace("iter = "+iter);//db

      ri = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);//compute initial residuals
      rmsri = rms(ri);
      xv = buildApproxHessianAndGradient(nc,kc,c,nb,kb,b,ri,u,g);
      //xv = buildScaledIandV(nc,kc,c,nb,kb,b,ri,u,g);
      x = xv[0];
      v = xv[1];
      //printXandV(x,v,"Original");
      //printDiagonal(x,"Original X Diagonal");

      //////////Stabilize//////////////
      x = stabilizeX(nb,nc,x);
      printDiagonal(x,"Stabilized X Diagonal");
      /////////////////////////////////
  
      //////////Penalization////////////
      xv = penalize(nb,kb,b,nc,kc,c,_alpha,_sfac,_simpleStab,_penb,_penc,x,v);
      x = xv[0];
      v = xv[1];
      //printXandV(x,v,"Penalized");
      //printDiagonal(x,"Penalized X Diagonal");
      //////////////////////////////////
    
      //////////Constrain deltab0 to be 0////////////
      xv = constrainDeltaB0ToBe0(nc,nb,kb,x,v);
      x = xv[0];
      v = xv[1];
      //printXandV(x,v,"Constrained");
      //printDiagonal(x,"Constrained X Diagonal");
      //////////////////////////////////

      deltabdeltac = solveForSearchDirection(nb,nc,x,v);
      deltab = copy(nb,0,deltabdeltac);
      deltac = copy(nc,nb,deltabdeltac);
      //trace("deltab");
      //dump(deltab);
      //trace("deltac");
      //dump(deltac);

      steplength = solveForStepLength(nb,kb,b,deltab,nc,kc,c,deltac,u,f,g);
      deltab = mul(steplength,deltab);
      deltac = mul(steplength,deltac);
      bnew = add(b,deltab);
      cnew = add(c,deltac);

      rf = computeDataResidual(nc,kc,cnew,nb,kb,bnew,u,f,g);
      rmsrf = rms(rf);
      rmspercentchange = percentChange(rmsrf,rmsri);

      //Debugging Information
      _twonormb[iter] = twoNorm(b);
      _twonormc[iter] = twoNorm(c);
      _lastiter = iter;
      _rmsri[iter] = rmsri;
      _rmsrf[iter] = rmsrf;
      _twonormdeltab[iter] = sqrt(sum(pow(deltab,2)));
      _twonormdeltac[iter] = sqrt(sum(pow(deltac,2)));
      _rmspercentchange[iter] = rmspercentchange;
      _condnum[iter]= (float) x.cond();
      _twonormgrad[iter] = (float) sqrt(sum(pow(v.getArray(),2.0)));
      _steplength[iter] = steplength;
      trace("initial rmsr: "+rmsri);
      trace("final rmsr: "+rmsrf);
      //trace("rmsr change: "+(rmsrf-rmsri));
      //trace("Step Length: "+steplength);
      //trace("2Norm of Negative Gradient: "+_twonormgrad[iter]);
      //trace("2Norm of deltab: "+_twonormdeltab[iter]);
      //trace("2Norm of deltac: "+_twonormdeltac[iter]);
      //trace("Smallest Eigenvalue: "+getFirstEigenvalue(x));
      trace("X Condition Number: "+_condnum[iter]);
      printB(b);
      printC(c);
      //printDeltaB(deltab);
      //printDeltaC(deltac);
      //plotResiduals(rf,iter,rmspercentchange);

      if (rmspercentchange>0.0f) {
        trace("Stopped Iterating: Percent change is positive");
        trace("Maximum percent change before stopping = "+_maxrmspercentchange);
        trace("Percent change = "+pc);
        _lastiter = iter-1;
        if (_lastiter<0)
          _lastiter=0;
        return new float[][]{c,b};
      }
      
      //Update c and b to newest c and b for the next iteration.
      b = copy(bnew);
      c = copy(cnew);
      if (abs(rmspercentchange)<_maxrmspercentchange) {
        trace("Stopped Iterating: Percent change is below"+
        " maximum percent change required for stopping");
        trace("Maximum percent change before stopping = "+_maxrmspercentchange);
        trace("Percent change = "+pc);
        _lastiter = iter;
        return new float[][]{c,b};
      }
      
    }
    return new float[][]{c,b};
  }
  
  public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bguess, int nc, int kc, float[] cguess,
    float[][] u, float[][] f, float[][] g, int niter)
  {
    trace("bguess");//db
    dump(bguess);//db
    trace("cguess");//db
    dump(cguess);//db
    _steplength = new float[niter];//db
    _rmsri = new float[niter];//db
    _rmsrf = new float[niter];//db
    _twonormgrad = new float[niter];//db
    _twonormdeltac = new float[niter];//db
    _twonormdeltab = new float[niter];//db
    _condnum = new float[niter];//db
    _twonormb = new float[niter];//db
    _twonormc = new float[niter];//db
    _rmspercentchange = new float[niter];//db

    int nf = f[0].length;
    int nx = f.length;
    int nbc = nc+nb;
    float[] b = copy(bguess);
    float[] c = copy(cguess);
    float[] bnew = new float[nb];
    float[] cnew = new float[nc];
    float[] deltabdeltac = new float[nbc];
    float[] deltab = new float[nb];
    float[] deltac = new float[nc];
    float[][] ri = new float[nx][nf];
    float[][] rf = new float[nx][nf];
    float rmsrf = 0.0f;
    float rmsri = 0.0f;
    float steplength = 0.0f;
    float pc = 0.0f;
    double alpha = 0.1;
    float rmspercentchange = 0.0f;
    boolean pendeltab = false;
    boolean pendeltac = false;
    boolean penb = false;
    boolean penc = false;
    DMatrix[] xv = new DMatrix[2];
    DMatrix v = new DMatrix(nbc,1);
    DMatrix x = new DMatrix(nbc,nbc);
    DMatrix deltabc = new DMatrix(nbc,1);

    //Gauss-Newton iterations
    for (int iter=0; iter<niter; ++iter) {
      trace("iter = "+iter);

      ri = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);//compute initial residuals
      rmsri = rms(ri);
      xv = buildApproxHessianAndGradient(nc,kc,c,nb,kb,b,ri,u,g);
      //xv = buildScaledIandV(nc,kc,c,nb,kb,b,ri,u,g);
      x = xv[0];
      v = xv[1];
      //printXandV(x,v,"Original");

      //////////Stabilize//////////////
      x = stabilizeX(nb,nc,x);
      /////////////////////////////////

      //////////Penalization////////////
      xv = penalize(nb,kb,b,nc,kc,c,_alpha,_sfac,_simpleStab,_penb,_penc,x,v);
      x = xv[0];
      v = xv[1];
      //printXandV(x,v,"Penalized");
      //////////////////////////////////
    
      //////////Constrain deltab0 to be 0////////////
      xv = constrainDeltaB0ToBe0(nc,nb,kb,x,v);
      x = xv[0];
      v = xv[1];
      //printXandV(x,v,"Constrained");
      //////////////////////////////////

      deltabdeltac = solveForSearchDirection(nb,nc,x,v);
      deltab = copy(nb,0,deltabdeltac);
      deltac = copy(nc,nb,deltabdeltac);

      steplength = solveForStepLength(nb,kb,b,deltab,nc,kc,c,deltac,u,f,g);
      deltab = mul(steplength,deltab);
      deltac = mul(steplength,deltac);
      bnew = add(b,deltab);
      cnew = add(c,deltac);

      rf = computeDataResidual(nc,kc,cnew,nb,kb,bnew,u,f,g);
      rmsrf = rms(rf);
      rmspercentchange = percentChange(rmsrf,rmsri);

      //Debugging Information
      _twonormb[iter] = twoNorm(b);
      _twonormc[iter] = twoNorm(c);
      _lastiter = iter;
      _rmsri[iter] = rmsri;
      _rmsrf[iter] = rmsrf;
      _twonormdeltab[iter] = sqrt(sum(pow(deltab,2)));
      _twonormdeltac[iter] = sqrt(sum(pow(deltac,2)));
      _rmspercentchange[iter] = rmspercentchange;
      _condnum[iter]= (float) x.cond();
      _twonormgrad[iter] = (float) sqrt(sum(pow(v.getArray(),2.0)));
      _steplength[iter] = steplength;
      trace("initial rmsr: "+rmsri);
      trace("final rmsr: "+rmsrf);
      trace("rmsr change: "+(rmsrf-rmsri));
      trace("Step Length: "+steplength);
      trace("2Norm of Negative Gradient: "+_twonormgrad[iter]);
      trace("2Norm of deltab: "+_twonormdeltab[iter]);
      trace("2Norm of deltac: "+_twonormdeltac[iter]);
      trace("Smallest Eigenvalue: "+getFirstEigenvalue(x));
      trace("X Condition Number: "+_condnum[iter]);
      printB(b);
      printC(c);
      printDeltaB(deltab);
      printDeltaC(deltac);
      //plotResiduals(rf,iter,rmspercentchange);

      if (rmspercentchange>0.0f) {
        trace("Stopped Iterating: Percent change is positive");
        trace("Maximum percent change before stopping = "+_maxrmspercentchange);
        trace("Percent change = "+pc);
        _lastiter = iter-1;
        if (_lastiter<0)
          _lastiter=0;
        return new float[][]{c,b};
      }
      
      //Update c and b to newest c and b for the next iteration.
      b = copy(bnew);
      c = copy(cnew);
      if (abs(rmspercentchange)<_maxrmspercentchange) {
        trace("Stopped Iterating: Percent change is below"+
        " maximum percent change required for stopping");
        trace("Maximum percent change before stopping = "+_maxrmspercentchange);
        trace("Percent change = "+pc);
        _lastiter = iter;
        return new float[][]{c,b};
      }
      
    }
    return new float[][]{c,b};
  }

  public int getLastIter() {
    return _lastiter;

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
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(caa);
    return stm.solve(ca1);
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
  public float[] getWaveletH(
    int nc, int kc, int nb, int kb, float[] b, float stabfact,
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
          qq.set(ic,jc,((1.0+(double) stabfact)*qq.get(ic,jc)));
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
  public float[] getWaveletH(
    int nc, int kc, int nb, int kb, float[] b, float stabfact,
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
    x = div(x,sum);
    return x;
  }


  public float[] getRMSRf() {
    return _rmsrf;
  }
  public float[] getTwoNormGrad() {
    return _twonormgrad;
  }
  public float[] getStepLength() {
    return _steplength;
  }
  public float[] getRMSRi() {
    return _rmsri;
  }
  public float[] getTwoNormDeltaB() {
    return _twonormdeltab;
  }
  public float[] getTwoNormDeltaC() {
    return _twonormdeltac;
  }
  public float[] getCondNum() {
    return _condnum;
  }
  public float[] getTwoNormB() {
    return _twonormb;
  }
  public float[] getTwoNormC() {
    return _twonormc;
  }
  public float[] getRMSPercentChange() {
    return _rmspercentchange;
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
   * Makes the rms equal to 1 within a specified time range.
   */
   public float[] makerms1(float[] x) {
     float rmsx = rms(x);
     float[] rms1x = div(x,rmsx);
     //Check
     trace("rms after normalization is "+rms(rms1x));
     return rms1x;
   }
   public float[][] makerms1(float[][] x) {
     float rmsx = rms(x);
     float[][] rms1x = div(x,rmsx);
     //Check
     trace("rms after normalization is "+rms(rms1x));
     return rms1x;
   }
   public float[] makerms1(int itmin, int itmax, float[] x) {
     float rmsx = rms(itmin,itmax,x);
     float[] rms1x = div(x,rmsx);
     //Check
     trace("rms after normalization is "+rms(itmin,itmax,rms1x));
     return rms1x;
   }
   public float[][] makerms1(int itmin, int itmax, float[][] x) {
     float rmsx = rms(itmin,itmax,x);
     float[][] rms1x = div(x,rmsx);
     //Check
     trace("rms after normalization is "+rms(itmin,itmax,rms1x));
     return rms1x;
   }
   public float[] computeDataResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[] u, float[] f, float[] g) 
  {
    Warper warp = new Warper();
    float[] csbg = applyH(nc,kc,c,warp.applyS(u,applyH(nb,kb,b,g)));
    return sub(csbg,f);
  }
  public float[][] computeDataResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[][] u, float[][] f, float[][] g) 
  {
    Warper warp = new Warper();
    float[][] csbg = applyH(nc,kc,c,warp.applyS(u,applyH(nb,kb,b,g)));
    return sub(csbg,f);
  }

  public float[] computeResidualTest(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[] u, float[] f, float[] g) 
  {
    SimplePlot sp0 = new SimplePlot();
    SimplePlot sp1 = new SimplePlot();
    SimplePlot sp2 = new SimplePlot();
    SimplePlot sp3 = new SimplePlot();
    SimplePlot sp4 = new SimplePlot();
    Warper warp = new Warper();
    float[] bg = applyH(nb,kb,b,g);
    float[] sbg = warp.applyS(u,bg);
    float[] csbg = applyH(nc,kc,c,sbg);
    sp0.addPoints(g);
    sp1.addPoints(bg);
    sp2.addPoints(sbg);
    sp3.addPoints(csbg);
    sp4.addPoints(sub(csbg,f));
    sp0.addTitle("g");
    sp1.addTitle("bg");
    sp2.addTitle("sbg");
    sp3.addTitle("csbg");
    sp4.addTitle("csbg-f");
    return sub(csbg,f);
  }

///////////////////////////////////////////////////////////////////////////
  // private
  private double _sfac=0.00;
  private float _maxrmspercentchange = 0.000f;
  private float _minsteplength = 0.0001f;
  private int _itmin = -1;
  private int _itmax = -1;
  private int _ng0 = 0;//staring size of array/image.
  private int _lastiter = 0;
  private double _alpha;
  private boolean _penb = false;
  private boolean _penc = false;
  private boolean _simpleStab = false;
  private float[] _rmsri;
  private float[] _rmsrf;
  private float[] _twonormgrad;
  private float[] _twonormdeltac;
  private float[] _twonormdeltab;
  private float[] _condnum;
  private float[] _steplength;
  private float[] _twonormb;
  private float[] _twonormc;
  private float[] _rmspercentchange;
  private boolean _parallel = false;

  private DMatrix stabilizeX(int nb, int nc, DMatrix x) {
    trace("sfac = "+_sfac);
    for (int ib=0; ib<nb; ++ib)
      x.set(ib,ib,x.get(ib,ib)*(1.0+_sfac));
    for (int ic=0; ic<nc; ++ic) {
      int nbpic = nb+ic;
      x.set(nbpic,nbpic,x.get(nbpic,nbpic)*(1.0+_sfac));
    }
    return x;
  }

  private float[] solveForSearchDirection(int nb, int nc, DMatrix x, DMatrix v) {
    //solve SPD linear system of equations for deltab and deltac.
    DMatrixChd chd = new DMatrixChd(x);
    DMatrix deltabdeltac = chd.solve(v);
    return convertDToF(deltabdeltac.getArray());

  }

  private float solveForStepLength(int nb, int kb, float[] b, float[] deltab,
    int nc, int kc, float[] c, float[] deltac,
    float[] u, float[] f, float[] g)
  {
    ErrorFunction ef = new ErrorFunction();
    ef.setParameters(nc,kc,c,deltac,nb,kb,b,deltab,u,f,g);
    BrentMinFinder bf = new BrentMinFinder(ef);
    return (float) bf.findMin(_minsteplength,1.0,_minsteplength);
  }
  private float solveForStepLength(int nb, int kb, float[] b, float[] deltab,
    int nc, int kc, float[] c, float[] deltac,
    float[][] u, float[][] f, float[][] g)
  {
    ErrorFunction2D ef = new ErrorFunction2D();
    ef.setParameters(nc,kc,c,deltac,nb,kb,b,deltab,u,f,g);
    BrentMinFinder bf = new BrentMinFinder(ef);
    return (float) bf.findMin(_minsteplength,1.0,_minsteplength);
  }

  private float[][] separateDeltaBDeltaC(int nb, int nc, DMatrix deltabc) {
    double[] tempbc = deltabc.getArray();
    float[] deltab = new float[nb];
    float[] deltac = new float[nc];
    for (int ib=0; ib<nb; ++ib) {
      deltab[ib] = (float) tempbc[ib];
    }
    for (int ic=0; ic<nc; ++ic) {
      deltac[ic] = (float) tempbc[nb+ic];
    }
    return new float[][]{deltab,deltac};
  }


  private DMatrix[] constrainDeltaB0ToBe0(int nc, int nb, int kb, DMatrix x, DMatrix v) {
    int nbc = nb+nc;
    double xkbxkb = x.get(-kb,-kb);
    for (int icb=0; icb<nbc; ++icb) {
      x.set(icb,-kb,0.0);
      x.set(-kb,icb,0.0);
    }
    x.set(-kb,-kb,xkbxkb);
    v.set(-kb,0,0.0);
    return new DMatrix[]{x,v};
  }
  
  private DMatrix[] penalize(int nb, int kb, float[] b, int nc, int kc, float[] c, 
    double alpha, double sfac, boolean simpleStab, 
    boolean penb, boolean penc, 
    DMatrix x, DMatrix v) 
  {
    if (penb == true) {
      DMatrix[] xv = penalizeB(nc,nb,kb,b,x,v,alpha);
      x = xv[0];
      v = xv[1];
    }
    if (penc == true) {
      DMatrix[] xv = penalizeC(nb,nc,kc,c,x,v,alpha);
      x = xv[0];
      v = xv[1];
    }
    return new DMatrix[]{x,v};
  }

  private DMatrix penalizeDeltaB(int nc, int nb, int kb, DMatrix x, double sfac) {
    int nbc = nc+nb;
    double[] ww = new double[nbc];
    for (int ib=0; ib<nb; ++ib) {
      int lag = abs(-kb-ib);
      double stabfactor = sfac*x.get(ib,ib); 
      ww[ib] = stabfactor;
    }
    return x.plus(DMatrix.diagonal(ww));
  }
  private DMatrix penalizeDeltaC(int nb, int nc, int kc, DMatrix x, double sfac) {
    int nbc = nc+nb;
    double[] zz = new double[nbc];
    for (int ic=0; ic<nc; ++ic) {
      int lag = abs(-kc-ic);
      double stabfactor= sfac*x.get(nb+ic,nb+ic); 
      zz[nb+ic] = stabfactor;
    }
    return x.plus(DMatrix.diagonal(zz));
  }

  private DMatrix[] penalizeB(int nc, int nb, int kb, float[] b, 
    DMatrix x, DMatrix v, double alpha) 
  {
    int nbc = nc+nb;
    double[] ww = new double[nbc];
    DMatrix newb = new DMatrix(nbc,1);
    for (int ib=0; ib<nb; ++ib) {
      int lag = abs(-kb-ib);
      double normalize = alpha*x.get(ib,ib); 
      ww[ib] = (double) (lag*lag)*normalize;
    }
    for (int ib=0; ib<nb; ++ib) {
      newb.set(ib,0,(double) b[ib]*ww[ib]);
    }
    DMatrix newx = x.plus(DMatrix.diagonal(ww));
    DMatrix newv = v.minus(newb);
    return new DMatrix[]{newx,newv};
  }

  private DMatrix[] penalizeC(int nb, int nc, int kc, float[] c, 
    DMatrix x, DMatrix v, double alpha) 
  {
    int nbc = nc+nb;
    double[] zz = new double[nbc];
    DMatrix newc = new DMatrix(nbc,1);
    for (int ic=0; ic<nc; ++ic) {
      int lag = abs(-kc-ic);
      double normalize = alpha*x.get(nb+ic,nb+ic); 
      zz[nb+ic] = (double) (lag*lag)*normalize;
    }
    for (int ic=0; ic<nc; ++ic) {
      newc.set(nb+ic,0,(double) c[ic]*zz[nb+ic]);
    }
    DMatrix newx = x.plus(DMatrix.diagonal(zz));
    DMatrix newv = v.minus(newc);
    return new DMatrix[]{newx,newv};
  }

  private class ErrorFunction implements BrentMinFinder.Function {
        public void setParameters(int nc, int kc, float[] c, float[] deltac, 
                             int nb, int kb, float[] b, float[] deltab,
                             float[] u, float[] f, float[] g){
          _nc = nc;
          _nb = nb;
          _kc = kc;
          _kb = kb;
          _c = copy(c);
          _b = copy(b);
          _deltac = copy(deltac);
          _deltab = copy(deltab);
          _u = u;
          _f = f;
          _g = g;
        }
        public double evaluate(double steplength) {
          float[] b = add(_b,mul((float) steplength,_deltab));
          float[] c = add(_c,mul((float) steplength,_deltac));
          float[] r = computeDataResidual(_nc,_kc,c,_nb,_kb,b,_u,_f,_g);
          //trace("steplength = "+steplength);
          //trace("dotrr = "+dot(r,r));
          return dot(r,r);
        }
        private float[] _c,_b,_deltac,_deltab,_u,_f,_g;
        private int _nc,_nb,_kc,_kb;
  }

  private class ErrorFunction2D implements BrentMinFinder.Function {
        public void setParameters(int nc, int kc, float[] c, float[] deltac, 
                             int nb, int kb, float[] b, float[] deltab,
                             float[][] u, float[][] f, float[][] g){
          _nc = nc;
          _nb = nb;
          _kc = kc;
          _kb = kb;
          _c = copy(c);
          _b = copy(b);
          _deltac = copy(deltac);
          _deltab = copy(deltab);
          _u = u;
          _f = f;
          _g = g;
        }
        public double evaluate(double steplength) {
          float[] b = add(_b,mul((float) steplength,_deltab));
          float[] c = add(_c,mul((float) steplength,_deltac));
          float[][] r = computeDataResidual(_nc,_kc,c,_nb,_kb,b,_u,_f,_g);
          return dot(r,r);
        }
        private float[] _c,_b,_deltac,_deltab;
        private float[][] _u,_f,_g;
        private int _nc,_nb,_kc,_kb;
  }



  
  private float[] computeApproxResidual(int nc, int kc, float[] c, float[] deltac, 
                                        int nb, int kb, float[] b, float[] deltab,
                                        float[] u, float[] f, float[] g) 
  {
    int nt = u.length;
    Warper warp = new Warper();
    float[] bg = applyH(nb,kb,b,g);
    float[] sbg = warp.applyS(u,bg);
    float[] csbg = applyH(nc,kc,c,sbg);
    float[] r = sub(csbg,f);
    float[] q0 = warp.applyS(u,bg);
    float[] pdeltab = new float[nt];
    float[] qdeltac = new float[nt];

    for (int ib=0,lagi=kb; ib<nb; ++ib,++lagi) {
      float[] dig = delay(lagi,g);
      float[] sdig = warp.applyS(u,dig);
      float[] pi = applyH(nc,kc,c,sdig);
      pdeltab = add(pdeltab,mul(pi,deltab[ib]));
    }

    for (int ic=0,lagi=kc; ic<nc; ++ic,++lagi) {
      float[] qi = delay(lagi,q0);
      qdeltac = add(qdeltac,mul(qi,deltac[ic]));
    }
    return add(add(r,pdeltab),qdeltac);
  }

  private float[][] computeApproxResidual(int nc, int kc, float[] c, float[] deltac, 
                                        int nb, int kb, float[] b, float[] deltab,
                                        float[][] u, float[][] f, float[][] g) 
  {
    int nt = u[0].length;
    int nx = u.length;
    Warper warp = new Warper();
    float[][] bg = applyH(nb,kb,b,g);
    float[][] sbg = warp.applyS(u,bg);
    float[][] csbg = applyH(nc,kc,c,sbg);
    float[][] r = sub(csbg,f);
    float[][] q0 = warp.applyS(u,bg);
    float[][] pdeltab = new float[nx][nt];
    float[][] qdeltac = new float[nx][nt];

    for (int ib=0,lagi=kb; ib<nb; ++ib,++lagi) {
      float[][] dig = delay(lagi,g);
      float[][] sdig = warp.applyS(u,dig);
      float[][] pi = applyH(nc,kc,c,sdig);
      pdeltab = add(pdeltab,mul(pi,deltab[ib]));
    }

    for (int ic=0,lagi=kc; ic<nc; ++ic,++lagi) {
      float[][] qi = delay(lagi,q0);
      qdeltac = add(qdeltac,mul(qi,deltac[ic]));
    }
    return add(add(r,pdeltab),qdeltac);
  }



  /**
   * v = [P'] 
   *     |--|r
   *     [Q']
   * r = CSBg-f
   *
   *
   */
  private DMatrix[] buildApproxHessianAndGradient(int nc, int kc, float[] c, 
                               int nb, int kb, float[] b, 
                               float[] r, float[] u, float[] g) 
  {
    int nt = u.length;
    int nbc = nb+nc;
    Warper warp = new Warper();
    DMatrix v = new DMatrix(nbc,1);
    DMatrix x = new DMatrix(nbc,nbc);
    float[] q0 = warp.applyS(u,applyH(nb,kb,b,g));
    float[] qi = new float[nt];
    float[] qj = new float[nt];
    float[] pi = new float[nt];
    float[] pj = new float[nt];
    double vi = 0.0;
    double xij = 0.0;
    //P'P and Q'P and P'Q
    for (int ib=0,lagi=kb; ib<nb; ++ib,++lagi) {
      pi = applyH(nc,kc,c,warp.applyS(u,delay(lagi,g)));
      vi = dot(pi,r);
      v.set(ib,0,-vi);
      for (int jb=0,lagj=kb; jb<=ib; ++jb,++lagj) {
        pj = applyH(nc,kc,c,warp.applyS(u,delay(lagj,g)));
        xij = dot(pi,pj);
        x.set(ib,jb,xij);
        x.set(jb,ib,xij);
      }
      for (int jc=0,lagj=kc; jc<nc; ++jc,++lagj) {
        qi = delay(lagj,q0);
        vi = dot(qi,r);
        v.set(nb+jc,0,-vi);
        xij = dot(qi,pi);
        x.set(nb+jc,ib,xij);
        x.set(ib,nb+jc,xij);
      }
    }

    //Q'Q
    for (int ic=nb,lagi=kc; ic<nbc; ++ic,++lagi) {
      qi = delay(lagi,q0);
      for (int jc=nb,lagj=kc; jc<=ic; ++jc,++lagj) {
        qj = delay(lagj,q0);
        if (ic==nb && jc==nb) {
          //SimplePlot sp = new SimplePlot();
          //sp.addPoints(qi);
          //sp.addTitle("Last Iter = "+_lastiter);
        }
        xij = dot(qi,qj);
        x.set(ic,jc,xij);
        x.set(jc,ic,xij);
      }
    }
    return new DMatrix[]{x,v};
  }
  private DMatrix[] buildApproxHessianAndGradient(int nc, int kc, float[] c, 
                               int nb, int kb, float[] b, 
                               float[][] r, float[][] u, float[][] g) 
  {
    int nt = u[0].length;
    int nx = u.length;
    int nbc = nb+nc;
    Warper warp = new Warper();
    DMatrix v = new DMatrix(nbc,1);
    DMatrix x = new DMatrix(nbc,nbc);
    float[][] q0 = warp.applyS(u,applyH(nb,kb,b,g));
    float[][] qi = new float[nx][nt];
    float[][] qj = new float[nx][nt];
    float[][] pi = new float[nx][nt];
    float[][] pj = new float[nx][nt];
    double vi = 0.0;
    double xij = 0.0;
    //P'P and Q'P and P'Q
    for (int ib=0,lagi=kb; ib<nb; ++ib,++lagi) {
      pi = applyH(nc,kc,c,warp.applyS(u,delay(lagi,g)));
      vi = dot(pi,r);
      v.set(ib,0,-vi);
      for (int jb=0,lagj=kb; jb<=ib; ++jb,++lagj) {
        pj = applyH(nc,kc,c,warp.applyS(u,delay(lagj,g)));
        xij = dot(pi,pj);
        x.set(ib,jb,xij);
        x.set(jb,ib,xij);
      }
      for (int jc=0,lagj=kc; jc<nc; ++jc,++lagj) {
        qi = delay(lagj,q0);
        vi = dot(qi,r);
        v.set(nb+jc,0,-vi);
        xij = dot(qi,pi);
        x.set(nb+jc,ib,xij);
        x.set(ib,nb+jc,xij);
      }
    }

    //Q'Q
    for (int ic=nb,lagi=kc; ic<nbc; ++ic,++lagi) {
      qi = delay(lagi,q0);
      for (int jc=nb,lagj=kc; jc<=ic; ++jc,++lagj) {
        qj = delay(lagj,q0);
        xij = dot(qi,qj);
        x.set(ic,jc,xij);
        x.set(jc,ic,xij);
      }
    }
    return new DMatrix[]{x,v};
  }

  private DMatrix[] buildScaledIandV(int nc, int kc, float[] c, 
                               int nb, int kb, float[] b, 
                               float[] r, float[] u, float[] g) 
  {
    int nbc = nb+nc;
    DMatrix x = new DMatrix(nbc,nbc);
    DMatrix v = new DMatrix(nbc,1);
    for (int ib=0,lagib=kb; ib<nb; ++ib,++lagib) {
      Warper warp = new Warper();
      float[] dig = delay(lagib,g);
      float[] sdig = warp.applyS(u,dig);
      float[] pi = applyH(nc,kc,c,sdig);
      double vi = dot(pi,r);
      double xij = dot(pi,pi);
      v.set(ib,0,-vi);
      x.set(ib,ib,xij);
    }
    Warper warp = new Warper();
    float[] ag = applyH(nb,kb,b,g);
    float[] q0 = warp.applyS(u,ag);
    for (int ic=0,lagic=kc; ic<nc; ++ic,++lagic) {
      float[] qi = delay(lagic,q0);
      double vi = dot(qi,r);
      double xij = dot(qi,qi);
      v.set(nb+ic,0,-vi);
      x.set(nb+ic,nb+ic,xij);
    }
    return new DMatrix[]{x,v};
  }
  private DMatrix[] buildScaledIandV(int nc, int kc, float[] c, 
                               int nb, int kb, float[] b, 
                               float[][] r, float[][] u, float[][] g) 
  {
    int nbc = nb+nc;
    DMatrix x = new DMatrix(nbc,nbc);
    DMatrix v = new DMatrix(nbc,1);
    for (int ib=0,lagib=kb; ib<nb; ++ib,++lagib) {
      Warper warp = new Warper();
      float[][] dig = delay(lagib,g);
      float[][] sdig = warp.applyS(u,dig);
      float[][] pi = applyH(nc,kc,c,sdig);
      double vi = dot(pi,r);
      double xij = dot(pi,pi);
      v.set(ib,0,-vi);
      x.set(ib,ib,xij);
    }
    Warper warp = new Warper();
    float[][] ag = applyH(nb,kb,b,g);
    float[][] q0 = warp.applyS(u,ag);
    for (int ic=0,lagic=kc; ic<nc; ++ic,++lagic) {
      float[][] qi = delay(lagic,q0);
      double vi = dot(qi,r);
      double xij = dot(qi,qi);
      v.set(nb+ic,0,-vi);
      x.set(nb+ic,nb+ic,xij);
    }
    return new DMatrix[]{x,v};
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

  private void printEigenvalues(DMatrix x) {
    DMatrixEvd evd = new DMatrixEvd(x);
    double[][] d = evd.getD().get();
    int nd = d.length;
    double[] diag = new double[nd];
    for (int id=0; id<nd; ++id) {
      diag[id] = d[id][id];
    }
    trace("Eigenvalues: ");
    dump(diag);
  }

  private float getFirstEigenvalue(DMatrix x) {
    DMatrixEvd evd = new DMatrixEvd(x);
    double[][] d = evd.getD().get();
    return (float) d[0][0];
  }
  private float getLastEigenvalue(DMatrix x) {
    DMatrixEvd evd = new DMatrixEvd(x);
    double[][] d = evd.getD().get();
    int n = d.length;
    return (float) d[n-1][n-1];
  }

  private float twoNorm(float[] x) {
    return sqrt(sum(pow(x,2.0f)));
  }
  private float twoNorm(float[][] x) {
    int nx = x.length;
    float[][] xsq = pow(x,2.0f);
    float sum = sum(xsq);
    return sqrt(sum);
  }

  private float percentChange(float xf, float xi) {
    return (xf-xi)/xi*100.0f;
  }

  private void printDeltaB(float[] deltab) {
    trace("deltab");
    dump(deltab);
  }
  private void printDeltaC(float[] deltac) {
    trace("deltac");
    dump(deltac);
  }
  private void printC(float[] c) {
    trace("c");
    dump(c);
  }
  private void printB(float[] b) {
    trace("b");
    dump(b);
  }
  /**
   * Prints diagonal of a square matrix.
   * @param x A square matrix.
   */
  private void printDiagonal(DMatrix x, String xname) {
    int ni = x.getM();
    float[] diag = new float[ni];
    for (int i=0; i<ni; ++i) {
      diag[i] = (float) x.get(i,i);
    }
    trace(xname);
    dump(diag);
  }
  private void printXandV(DMatrix x,DMatrix v,String type) {
      trace(type+" X");
      trace(x.toString());
      trace(type+" V");
      trace(v.toString());
  }
  private void plotResiduals(float[] rf,int iter,float pc) {
    SimplePlot sp = new SimplePlot();
    sp.addPoints(rf);
    sp.addTitle("Iteration = "+iter+" PercentChange = "+pc);
  }
  private void plotResiduals(float[][] rf,int iter,float pc) {
    SimplePlot sp = new SimplePlot();
    sp.addPixels(rf);
    sp.addTitle("Iteration = "+iter+" PercentChange = "+pc);
  }

  private static void trace(String s) {
    System.out.println(s);
  }
 }



