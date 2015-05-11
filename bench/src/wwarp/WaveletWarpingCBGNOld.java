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
 public class WaveletWarpingCBGNOld {
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
   * @param pcmax The maximum percent change needed to stop this iterative method.
   */
  public void setMaxPercentChange(float pcmax) {
    _pcmax = pcmax;
  }



  /**
   * Sets the stability factor by which to steplength the diagonal of the matrix X.
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
    int ncb = nc+nb;
    float[] b = copy(bguess);
    float[] c = copy(cguess);
    trace("bguess");
    dump(b);
    trace("cguess");
    dump(c);
    float[] bnew = new float[nb];
    float[] cnew = new float[nc];
    _rapprapp = new float[niter];
    _steplength = new float[niter];
    _riri = new float[niter];
    _rfrf = new float[niter];
    _vfvf = new float[niter];
    _deltamagc = new float[niter];
    _deltamagb = new float[niter];
    _condnum = new float[niter];
    _twonormb = new float[niter];
    _twonormc = new float[niter];

    //Gauss-Newton iterations
    for (int iter=0; iter<niter; ++iter) {
      trace("iter = "+iter);
      //compute initial residual
      float[] ri = computeResidual(nc,kc,c,nb,kb,b,u,f,g);
      float riri = rms(ri);
      _riri[iter] = riri;

      //Build X and V
      DMatrix[] xv = buildApproxHessianAndGradient(nc,kc,c,nb,kb,b,ri,u,g);
      //DMatrix[] xv = buildScaledIandV(nc,kc,c,nb,kb,b,ri,u,g);
      DMatrix x = xv[0];
      trace("Original x");
      trace(x.toString());
      DMatrix v = xv[1];
      float[] vf = convertDToF(v.getArray());
      float vfvf = twoNorm(vf);
      _vfvf[iter] = vfvf;

      //Add in constraint db0=0
      double xkbkb = x.get(-kb,-kb);
      v.set(-kb,0,0.0);
      for (int icb=0; icb<ncb; ++icb) {
        x.set(icb,-kb,0.0);
        x.set(-kb,icb,0.0);
      }
      x.set(-kb,-kb,xkbkb);
     
      //Stabilize x (Reduce rounding errors by multiplying 
      //diagonal of matrix X by some number above 1.0.
      double sc = 1.0+_sfac;
      trace("sc = "+sc);
      for (int icb=0; icb<ncb; ++icb) {
        double xii = x.get(icb,icb)*sc;
        x.set(icb,icb,xii);
      }
      trace("Stabilize x");
      trace(x.toString());

      /*//Penalize large deltabs away from zero lag
      double[] z = new double[ncb];
      double scale = abs(1.0/min(c));
      for (int ib=0; ib<nb; ++ib) {
        z[ib] = (double) (scale*abs(-kb-ib));
        trace("z = "+z[ib]);
      }*/
      //trace(DMatrix.diagonal(z).toString());
      //trace("before penalty");
      //trace(x.toString());
      //x = x.plus(DMatrix.diagonal(z));
      //trace("after penalty");
      //trace(x.toString());
      //trace("v");
      //trace(v.toString());

      //condition number
      float condnum = (float) x.cond();
      _condnum[iter] = condnum;

      //solve SPD linear system of equations for deltab and deltac.
      float[] deltab = new float[nb];
      float[] deltac = new float[nc];
      DMatrixChd chd = new DMatrixChd(x);
      DMatrix deltabc = chd.solve(v);
      trace("deltabc:");
      trace(deltabc.toString());
      double[] tempbc = deltabc.getArray();
      for (int ib=0; ib<nb; ++ib) {
        deltab[ib] = (float) tempbc[ib];
      }
      for (int ic=0; ic<nc; ++ic) {
        deltac[ic] = (float) tempbc[nb+ic];
      }

      ErrorFunction ef = new ErrorFunction();
      ef.setParameters(nc,kc,c,deltac,nb,kb,b,deltab,u,f,g);
      BrentMinFinder bf = new BrentMinFinder(ef);
      float steplength = (float) bf.findMin(_minsteplength,1.0,_minsteplength);
      _steplength[iter] = steplength;
      bnew = add(b,mul(steplength,deltab));
      cnew = add(c,mul(steplength,deltac));
      float[] rf = computeResidual(nc,kc,cnew,nb,kb,bnew,u,f,g);
      float rfrf = rms(rf);
      _rfrf[iter] = rfrf;
      /*if (iter==100) { 
        _r200 = rf;
        _rfrf200 = rms(rf);
        _c200 = cnew;
        _b200 = bnew;
      }
      if (iter==150) { 
        _r900 = rf;
        _rfrf900 = rms(rf);
        _c900 = cnew;
        _b900 = bnew;

        int ng = g.length;
        float[] bg200 = applyH(nb,kb,_b200,g);
        float[] bg900 = applyH(nb,kb,_b900,g);
        SimplePlot sp = new SimplePlot();
        SimplePlot sp1 = new SimplePlot();
        SimplePlot sp2 = new SimplePlot();
        SimplePlot sp3 = new SimplePlot();
        SimplePlot sp4 = new SimplePlot();
        PointsView pv7 = sp3.addPoints(div(bg200,bg900));
        pv7.setLineColor(Color.BLACK);
        PointsView pv8 = sp4.addPoints(makerms1(bg200));
        pv8.setLineColor(Color.BLACK);
        PointsView pv9 = sp4.addPoints(makerms1(bg900));
        pv9.setLineColor(Color.BLUE);
        PointsView pv1 = sp.addPoints(sub(_r200,_r900));
        pv1.setLineColor(Color.BLACK);
        PointsView pv2 = sp1.addPoints(makerms1(0,nb-1,_b200));
        pv2.setLineColor(Color.BLACK);
        PointsView pv3 = sp2.addPoints(makerms1(0,nc-1,_c200));
        pv3.setLineColor(Color.BLACK);
        PointsView pv5 = sp1.addPoints(makerms1(0,nb-1,_b900));
        pv5.setLineColor(Color.BLUE);
        PointsView pv6 = sp2.addPoints(makerms1(0,nc-1,_c900));
        pv6.setLineColor(Color.BLUE);
        sp.addTitle("r 200 (BLACK) - 900 (BLUE)");
        sp1.addTitle("rms is 1 b 200 (BLACK) 900 (BLUE)");
        sp2.addTitle("rms is 1 c 200 (BLACK) 900 (BLUE)");
        sp3.addTitle("rms is 1, Bg 200 (BLACK) / Bg 900 (BLUE)");
        sp4.addTitle("rms is 1, Bg 200 (BLACK), Bg 900 (BLUE)");
      }
      */

      trace("initial r'r , final r'r , 2norm neg grad , condnum , 1st eigval , steplength");
      trace(riri+" , "+rfrf+" , "+" , "+vfvf+" , "+condnum+" , "+getFirstEigenvalue(x)+" , "+steplength);
      trace("x00 = "+x.get(0,0));
      trace("xnbcnbc = "+x.get(nb+nc-1,nb+nc-1));

      _deltamagb[iter] = sqrt(sum(pow(deltab,2)));
      _deltamagc[iter] = sqrt(sum(pow(deltac,2)));
      float pc = percentChange(rfrf,riri);
      if (pc>0.0f) {
        trace("Stopped Iterating: Percent change is positive");
        trace("Maximum percent change before stopping = "+_pcmax);
        trace("Percent change = "+pc);
        _lastiter = iter-1;
        if (_lastiter<0)
          _lastiter=0;
        trace("last iter = "+_lastiter);
        return new float[][]{c,b};
      }
      
      //Update c and b to newest c and b for the next iteration.
      b = copy(bnew);
      c = copy(cnew);
      if (abs(pc)<_pcmax) {
        trace("Stopped Iterating: Percent change is below"+
        " maximum percent change required for stopping");
        trace("Maximum percent change before stopping = "+_pcmax);
        trace("Percent change = "+pc);
        _lastiter = iter;
        trace("last iter = "+_lastiter);
        return new float[][]{c,b};
      }
      _twonormb[iter] = twoNorm(b);
      _twonormc[iter] = twoNorm(c);
      _lastiter = iter;
      trace("last iter = "+_lastiter);
      trace("clast");
      dump(c);
      trace("deltac");
      dump(deltac);
      trace("blast");
      dump(b);
      trace("deltab");
      dump(deltab);
    }
    /*
    trace("rfrf200 = "+_rfrf200);
    trace("rfrf900 = "+_rfrf900);
    trace("c200");
    dump(_c200);
    trace("c900");
    dump(_c900);
    trace("b200");
    dump(_b200);
    trace("b900");
    dump(_b900);
    */

    return new float[][]{c,b};
  }

   public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bguess, int nc, int kc, float[] cguess,
    float[][] u, float[][] f, float[][] g, int niter)
  {
    int ncb = nc+nb;
    float[] b = copy(bguess);
    float[] c = copy(cguess);
    float[] bnew = new float[nb];
    float[] cnew = new float[nc];
    _rapprapp = new float[niter];
    _steplength = new float[niter];
    _riri = new float[niter];
    _rfrf = new float[niter];
    _vfvf = new float[niter];
    _deltamagc = new float[niter];
    _deltamagb = new float[niter];
    _condnum = new float[niter];
    _twonormb = new float[niter];
    _twonormc = new float[niter];

    //Gauss-Newton iterations
    for (int iter=0; iter<niter; ++iter) {
      trace("iter = "+iter);
      //compute initial residual
      trace("Compute initial residual");
      float[][] ri = computeResidual(nc,kc,c,nb,kb,b,u,f,g);
      float riri = rms(ri);
      _riri[iter] = riri;
      trace("End Compute initial residual");

      //Build X (approximated Hessian) and V (gradient)
      trace("Build X and v");
      DMatrix[] xv = buildApproxHessianAndGradient(nc,kc,c,nb,kb,b,ri,u,g);
      //DMatrix[] xv = buildScaledIandV(nc,kc,c,nb,kb,b,ri,u,g);
      trace("ri");
      trace("maxri = "+max(ri));
      trace("meanri = "+(sum(ri)/(ri.length*ri[0].length)));
      dump(ri[0]);
      trace(xv[0].toString());
      DMatrix x = xv[0];
      DMatrix v = xv[1];
      trace("End Build X and v");
      float[] vf = convertDToF(v.getArray());
      float vfvf = twoNorm(vf);
      _vfvf[iter] = vfvf;

      //Add in constraint db0=0
      trace("Add in constraint");
      double xkbkb = x.get(-kb,-kb);
      v.set(-kb,0,0.0);
      for (int icb=0; icb<ncb; ++icb) {
        x.set(icb,-kb,0.0);
        x.set(-kb,icb,0.0);
      }
      x.set(-kb,-kb,xkbkb);
      trace("End Add in constraint");
    
      trace("BeforeStabfirstEigenvalue = "+getFirstEigenvalue(x));
      //Stabilize x (Reduce rounding errors by multiplying 
      //diagonal of matrix X by some number above 1.0.
      trace("stabilization");
      double sc = 1.0+_sfac;
      trace("sc = "+sc);
      for (int icb=0; icb<ncb; ++icb) {
        double xii = x.get(icb,icb)*sc;//+sc;
        x.set(icb,icb,xii);
      }
      trace("end stabilization");

      //condition number
      float condnum = (float) x.cond();
      _condnum[iter] = condnum;
      //trace(x.toString());
      //trace(v.toString());
      trace("AfterStabfirstEigenvalue = "+getFirstEigenvalue(x));

      //solve SPD linear system of equations for deltab and deltac.
      float[] deltab = new float[nb];
      float[] deltac = new float[nc];
      trace("Cholesky Decomp");
      DMatrixChd chd = new DMatrixChd(x);
      DMatrix deltabc = chd.solve(v);
      double[] tempbc = deltabc.getArray();
      for (int ib=0; ib<nb; ++ib) {
        deltab[ib] = (float) tempbc[ib];
      }
      for (int ic=0; ic<nc; ++ic) {
        deltac[ic] = (float) tempbc[nb+ic];
      }
      trace("End Cholesky Decomp");

      trace("Line search");
      ErrorFunction2D ef = new ErrorFunction2D();
      ef.setParameters(nc,kc,c,deltac,nb,kb,b,deltab,u,f,g);
      BrentMinFinder bf = new BrentMinFinder(ef);
      float steplength = (float) bf.findMin(_minsteplength,1.0,0.00001);
      _steplength[iter] = steplength;
      bnew = add(b,mul(steplength,deltab));
      cnew = add(c,mul(steplength,deltac));
      trace("End Line search");

      float[][] rf = computeResidual(nc,kc,cnew,nb,kb,bnew,u,f,g);
      float rfrf = rms(rf);
      _rfrf[iter] = rfrf;
      trace("initial r'r , final r'r , 2norm neg grad , condnum , 1st eigval , steplength");
      trace(riri+" , "+rfrf+" , "+" , "+vfvf+" , "+condnum+" , "+getFirstEigenvalue(x)+" , "+steplength);

      trace("x00 = "+x.get(0,0));
      trace("xnbcnbc = "+x.get(nb+nc-1,nb+nc-1));
      _deltamagb[iter] = sqrt(sum(pow(deltab,2)));
      _deltamagc[iter] = sqrt(sum(pow(deltac,2)));
      float pc = percentChange(rfrf,riri);
      if (pc>0.0f) {
        trace("Stopped Iterating: Percent change is positive");
        trace("Maximum percent change before stopping = "+_pcmax);
        trace("Percent change = "+pc);
        _lastiter = iter-1;
        trace("last iter = "+_lastiter);
        return new float[][]{c,b};
      }
      //Update c and b to newest c and b for the next iteration.
      b = copy(bnew);
      c = copy(cnew);
      //if (abs(pc)<_pcmax && _steplengthgreater0pt9) {
      if (abs(pc)<_pcmax) {
        //trace("Steplength is above 0.9");
        trace("Stopped Iterating: Percent change is below"+
        " maximum percent change required for stopping");
        trace("Maximum percent change before stopping = "+_pcmax);
        trace("Percent change = "+pc);
        _lastiter = iter;
        trace("last iter = "+_lastiter);
        return new float[][]{c,b};
      }
      _twonormb[iter] = twoNorm(b);
      _twonormc[iter] = twoNorm(c);
      _lastiter = iter;
      trace("last iter = "+_lastiter);
      trace("clast");
      dump(c);
      trace("blast");
      dump(b);
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
    //caa[0] *= _sfac;
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
    int nc, int kc, int nb, int kb, float[] b,
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
    int nc, int kc, int nb, int kb, float[] b,
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


  public float[] getRfRf() {
    return _rfrf;
  }
  public float[] getVV() {
    return _vfvf;
  }
  public float[] getStepLength() {
    return _steplength;
  }
  public float[] getRappRapp() {
    return _rapprapp;
  }
  public float[] getRiRi() {
    return _riri;
  }
  public float[] getDeltaMagB() {
    return _deltamagb;
  }
  public float[] getDeltaMagC() {
    return _deltamagc;
  }
  public float[] getCondNum() {
    return _condnum;
  }
  public float[] getTwoNormb() {
    return _twonormb;
  }
  public float[] getTwoNormc() {
    return _twonormc;
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
   public float[] computeResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[] u, float[] f, float[] g) 
  {
    Warper warp = new Warper();
    float[] bg = applyH(nb,kb,b,g);
    float[] sbg = warp.applyS(u,bg);
    float[] csbg = applyH(nc,kc,c,sbg);
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
  public float[][] computeResidual(int nc, int kc, float[] c, int nb, int kb, float[] b,
    float[][] u, float[][] f, float[][] g) 
  {
    Warper warp = new Warper();
    float[][] csbg = applyH(nc,kc,c,warp.applyS(u,applyH(nb,kb,b,g)));
    return sub(csbg,f);
  }

  



///////////////////////////////////////////////////////////////////////////
  // private
  private float _pcmax = 0.000f;
  private float _minsteplength = 0.0001f;
  private double _sfac = 1.0;
  private int _itmin = -1;
  private int _itmax = -1;
  private int _ng0 = 0;//staring size of array/image.
  private int _lastiter = 0;
  private float[] _rfrf;
  private float[] _vfvf;
  private float[] _deltamagc;
  private float[] _deltamagb;
  private float[] _rapprapp;
  private float[] _condnum;
  private float[] _riri;
  private float[] _steplength;
  private float[] _a;
  private float[] _twonormb;
  private float[] _twonormc;
  private boolean _parallel = false;
  private boolean _steplengthgreater0pt9 = false;
  private float[] _r200;
  private float _rfrf200;
  private float[] _c200;
  private float[] _b200;
  private float[] _r900;
  private float _rfrf900;
  private float[] _c900;
  private float[] _b900;

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
          float[] r = computeResidual(_nc,_kc,c,_nb,_kb,b,_u,_f,_g);
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
          float[][] r = computeResidual(_nc,_kc,c,_nb,_kb,b,_u,_f,_g);
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
    //P'P and Q'P and P'Q
    for (int ib=0,lagi=kb; ib<nb; ++ib,++lagi) {
      float[] pi = applyH(nc,kc,c,warp.applyS(u,delay(lagi,g)));
      double vi = dot(pi,r);
      v.set(ib,0,-vi);
      for (int jb=0,lagj=kb; jb<=ib; ++jb,++lagj) {
        float[] pj = applyH(nc,kc,c,warp.applyS(u,delay(lagj,g)));
        double xij = dot(pi,pj);
        x.set(ib,jb,xij);
        x.set(jb,ib,xij);
      }
      for (int jc=0,lagj=kc; jc<nc; ++jc,++lagj) {
        float[] qi = delay(lagj,q0);
        vi = dot(qi,r);
        v.set(nb+jc,0,-vi);
        double xij = dot(qi,pi);
        x.set(nb+jc,ib,xij);
        x.set(ib,nb+jc,xij);
      }
    }

    //Q'Q
    for (int ic=nb,lagi=kc; ic<nbc; ++ic,++lagi) {
      float[] qi = delay(lagi,q0);
      for (int jc=nb,lagj=kc; jc<=ic; ++jc,++lagj) {
        float[] qj = delay(lagj,q0);
        double xij = dot(qi,qj);
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
    int nbc = nb+nc;
    DMatrix v = new DMatrix(nbc,1);
    DMatrix x = new DMatrix(nbc,nbc);
    Warper warp1 = new Warper();
    Warper warp3 = new Warper();
    Warper warp4 = new Warper();
    float[][] q0 = warp1.applyS(u,applyH(nb,kb,b,g));
    float sump = 0.0f;
    float sumq = 0.0f;
    float maxp = 0.0f;
    float maxq = 0.0f;

    //P'r and Q'r and P'Q
    for (int ib=0; ib<nb; ++ib) {
        int lagi = kb+ib;
        float[][] pi = applyH(nc,kc,c,warp1.applyS(u,delay(lagi,g)));
        sump = sump + sum(pi);
        if (max(pi)>maxp)
          maxp = max(pi);
        double vi = dot(pi,r);
        v.set(ib,0,-vi);
        for (int jb=0,lagj=kb; jb<=ib; ++jb,++lagj) {
          //float[][] pj = zerofloat(nt,u.length);
          float[][] pj = applyH(nc,kc,c,warp1.applyS(u,delay(lagj,g)));
          double xij = dot(pi,pj);
          x.set(ib,jb,xij);
          x.set(jb,ib,xij);
        }
        for (int jc=0,lagj=kc; jc<nc; ++jc,++lagj) {
          float[][] qi = delay(lagj,q0);
          sumq = sumq + sum(qi);
          if (max(qi)>maxq)
            maxq = max(qi);
          vi = dot(qi,r);
          v.set(nb+jc,0,-vi);
          double xij = dot(qi,pi);
          x.set(nb+jc,ib,xij);
          x.set(ib,nb+jc,xij);
        }
    }
    trace("meanp = "+sump/(q0.length*q0[0].length*nb));
    trace("meanq = "+sumq/(q0.length*q0[0].length*nc));
    trace("maxp = "+maxp);
    trace("maxq = "+maxq);

    //Q'Q
    for (int ic=nb,lagi=kc; ic<nbc; ++ic,++lagi) {
      float[][] qi = delay(lagi,q0);
      for (int jc=nb,lagj=kc; jc<=ic; ++jc,++lagj) {
        float[][] qj = delay(lagj,q0);
        double xij = dot(qi,qj);
        x.set(ic,jc,xij);
        x.set(jc,ic,xij);
      }
    }
    //trace(x.toString());

    //scale x and v by the maximum value in v
    /*
    double maxv = max(v.getArray());
    double omaxv = 1.0/maxv;
    x = x.times(omaxv);
    v = v.times(omaxv);
    */

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

  /*
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
  */

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
    int nx = x.length;
    float[] xsq = pow(x,2.0f);
    float sum = sum(xsq);
    return sqrt(sum);
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


  private static void trace(String s) {
    System.out.println(s);
  }
 }


