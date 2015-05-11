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
 * distortion enables this algorithm to estimate the wavelets. To estimate
 * these wavelets, a variation in the amount of warping should exist.
 *
 * For images, convolution with the wavelet is assumed to be in only 
 * the 1st dimension. For definiteness, this 1st dimension is assumed
 * to be time.
 *
 * @author Chris Graziano, Colorado School of Mines
 * @version 2013.05.10
 */
public class WaveletWarpingCBGN {

  /**
   * Sets the min-max range of times used 
   * to estimate wavelets and all residual 
   * calculations.
   * @param itmin minimum time, in samples.
   * @param itmax maximum time, in samples.
   */
  public void setTimeRange(int itmin, int itmax) 
  {
    _itmin = itmin;
    _itmax = itmax;
  }

  /**
   * The Gauss-Newton method is based on estimating a search direction and 
   * a step length. The newest update of the wavelets c and b
   * will be the search direction scaled by the step length. 
   * A line search along the search direction determines the 
   * step length.
   * This method defines the lower bound of values to be tested in the line
   * search. 
   * This lower bound also serves as the tolerance value for the line search 
   * algorithm. A value will be returned by the line search algorithm when
   * the distance between the current value of the line search and the true minimum is 
   * smaller than the tolerance value.
   * The upper bound of the values to be searched by the line search is set to be 1.0.
   * The default lower bound and tolerance value is 0.0.
   * @param min minimum value in line search.
   */
  public void setLineSearchMin(double min) 
  {
    _minStepLength = min;
  }

  /**
   * Sets the smallest percent change of the rms of the residuals before the
   * Gauss-Newton method will stop iterating. The default percent change required before 
   * the Gauss-Newton method will stop iterating is 0.01 percent.
   * Note that you should NOT include the negative in the percent change, we are assuming 
   * that you wish to have a decrease in the residuals.
   * @param minRmsPercentChange The minimum percent change needed to stop this iterative method.
   */
  public void setMinPercentChange(float minRmsPercentChange) {
    _minRmsPercentChange = minRmsPercentChange;
  }

  /**
   * Estimates the wavelet c and the inverse wavelet b using the Gauss-Newton algorithm.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index of b[0].
   * @param bguess array of coefficients for the first guess of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index of c[0].
   * @param cguess array of coefficients for the first guess of the wavelet c.
   * @param niter number of iterations the Gauss-Newton method will have.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   */
  public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bGuess, int nc, int kc, float[] cGuess,
    float[] u, float[] f, float[] g, int niter)
  {
    Check.argument(-nb<kb,"-nb<kb");
    Check.argument(kb<=0,"kb<=0");
    Check.argument(-nc<kc,"-nc<kc");
    Check.argument(kc<=0,"kc<=0");

    //trace("bGuess");//db
    //dump(bGuess);//db
    //trace("cGuess");//db
    //dump(cGuess);//db

    _dataResRmsInitS = new float[niter];//db
    _dataResRmsFinaS = new float[niter];//db
    _dataRes2NormSqInitS = new float[niter];//db
    _dataRes2NormSqFinaS = new float[niter];//db
    _bPenaltyRes2NormSqInitS = new float[niter];//db
    _bPenaltyRes2NormSqFinaS = new float[niter];//db
    _cPenaltyRes2NormSqInitS = new float[niter];//db
    _cPenaltyRes2NormSqFinaS = new float[niter];//db
    _allResRmsInitS = new float[niter];//db
    _allResRmsFinaS = new float[niter];//db
    _allResRmsAllS = new float[niter];//db
    _allRes2NormSqInitS = new float[niter];//db
    _allRes2NormSqFinaS = new float[niter];//db
    _alphaBS = new float[niter];
    _alphaCS = new float[niter];

    _stepLengthS = new float[niter];//db
    _conditionNumberS = new float[niter];//db
    _b2NormS = new float[niter];//db
    _c2NormS = new float[niter];//db
    _deltaB2NormS = new float[niter];//db
    _deltaC2NormS = new float[niter];//db
    _gradient2NormS = new float[niter];//db
    _rmsPercentChangeS = new float[niter];//db

    _zDiag= new double[nc];
    _wDiag = new double[nb];
    _alphaB = 0.0;
    _alphaC = 0.0;
    _scaleW = sqrt(dot(f,f));
    _scaleZ = sqrt(dot(f,f));

    int nbc = nb+nc;
    float[] b = copy(bGuess);
    float[] c = copy(cGuess);
    float[] bNew = new float[nb];
    float[] cNew = new float[nc];
    float[] deltaBDeltaC = new float[nbc];
    float[] deltaB = new float[nb];
    float[] deltaC = new float[nc];
    float[] dataResInit = new float[f.length];
    float[] bPenaltyResInit = new float[nb];
    float[] cPenaltyResInit = new float[nc];
    float[] dataResFina = new float[f.length];
    float[] bPenaltyResFina = new float[nb];
    float[] cPenaltyResFina = new float[nc];
    float[] minC = new float[nc];
    float[] minB = new float[nb];
    float minAllResRmsFina = Float.MAX_VALUE;
    int minIter = 0;
    float allResRmsInit = 0.0f;
    float allResRmsFina = 0.0f;
    float allRes2NormSqInit = 0.0f;
    float allRes2NormSqFina = 0.0f;
    float dataResRmsInit = 0.0f;
    float dataResRmsFina = 0.0f;
    float dataRes2NormSqInit = 0.0f;
    float dataRes2NormSqFina = 0.0f;
    float bPenaltyRes2NormSqInit = 0.0f;
    float bPenaltyRes2NormSqFina = 0.0f;
    float cPenaltyRes2NormSqInit = 0.0f;
    float cPenaltyRes2NormSqFina = 0.0f;
    float stepLength = 0.0f;
    float b2Norm = 0.0f;
    float c2Norm = 0.0f;
    float gradient2Norm = 0.0f;
    float rmsPercentChange = 0.0f;
    DMatrix[] xv = new DMatrix[2];//Approximated Hessian (x) gradient (v)
    DMatrix v = new DMatrix(nbc,1);
    DMatrix x = new DMatrix(nbc,nbc);

    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
    //Gauss-Newton iterations
    for (int iter=0; iter<niter; ++iter) {
      trace("###############iter = "+iter+"#####################");//db
      trace("alphaB = "+_alphaB);
      trace("current b:");
      dump(b);
      trace("current c:");
      dump(c);
      trace("W diag:");
      dump(_wDiag);
      trace("Z diag:");
      dump(_zDiag);
      //Compute Initial Residuals
      dataResInit = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
      bPenaltyResInit = computeBPenaltyResidual(convertFToD(b),_wDiag);
      trace("bPenaltyResInit:");
      dump(bPenaltyResInit);
      cPenaltyResInit = computeCPenaltyResidual(convertFToD(c),_zDiag);
      //Compute Initial Measures of Residuals
      dataResRmsInit = rms(dataResInit);
      dataRes2NormSqInit = (float) dot(dataResInit,dataResInit);
      bPenaltyRes2NormSqInit = (float) dotAll(bPenaltyResInit,bPenaltyResInit);
      cPenaltyRes2NormSqInit = (float) dotAll(cPenaltyResInit,cPenaltyResInit);
      allResRmsInit = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
      trace("dataRes2NormSqInit = "+dataRes2NormSqInit);
      trace("bPenaltyRes2NormSqInit = "+bPenaltyRes2NormSqInit);
      trace("cPenaltyRes2NormSqInit = "+cPenaltyRes2NormSqInit);
      allRes2NormSqInit = dataRes2NormSqInit+bPenaltyRes2NormSqInit+cPenaltyRes2NormSqInit;

      if (iter==0)
        _firstAll2NormSq = allRes2NormSqInit;
      _firstAll2NormSqTemp = allRes2NormSqInit; 
      
      ////////////Build Hessian and gradient/////////////////
      xv = buildApproxHessianAndGradient(nc,kc,c,nb,kb,b,dataResInit,u,g);
      x = xv[0];
      v = xv[1];

      ///////////Add Penalization/////////////////////
      xv = penalize(_wDiag,_zDiag,convertFToD(b),convertFToD(c),x,v);
      x = xv[0];
      v = xv[1];
    
      //////////Constrain deltab0 to be 0////////////
      xv = constrainDeltaB0ToBe0(nc,nb,kb,x,v);
      x = xv[0];
      v = xv[1];

      //Solve for deltab and deltac
      deltaBDeltaC = solveForSearchDirection(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleW);
      //deltaBDeltaC = solveForSearchDirection(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleZ);
      deltaB = copy(nb,0,deltaBDeltaC);
      deltaC = copy(nc,nb,deltaBDeltaC);

      //Solve for scaled deltab and deltac
      stepLength = solveForStepLength(nb,kb,b,deltaB,nc,kc,c,deltaC,u,f,g,_wDiag,_zDiag);
      deltaB = mul(stepLength,deltaB);
      deltaC = mul(stepLength,deltaC);
      bNew = add(b,deltaB);
      cNew = add(c,deltaC);

      //Compute Final Residuals
      dataResFina = computeDataResidual(nc,kc,cNew,nb,kb,bNew,u,f,g);
      bPenaltyResFina = computeBPenaltyResidual(convertFToD(bNew),_wDiag);
      cPenaltyResFina = computeCPenaltyResidual(convertFToD(cNew),_zDiag);
      //Compute Final Measures of Residuals
      dataResRmsFina = rms(dataResFina);
      dataRes2NormSqFina = (float) dot(dataResFina,dataResFina);
      bPenaltyRes2NormSqFina = (float) dotAll(bPenaltyResFina,bPenaltyResFina);
      cPenaltyRes2NormSqFina = (float) dotAll(cPenaltyResFina,cPenaltyResFina);
      allResRmsFina = rmsOfObjectiveFunction(dataResFina,bPenaltyResFina,cPenaltyResFina);
      allRes2NormSqFina = dataRes2NormSqFina+bPenaltyRes2NormSqFina+cPenaltyRes2NormSqFina;
      rmsPercentChange = percentChange(allResRmsFina,allResRmsInit);

      //Debugging Information
      if (iter == 0)
        _allResRmsAllS[0] = allResRmsInit;
      else
        _allResRmsAllS[iter] = allResRmsInit;
      _allResRmsAllS[iter] = allResRmsInit;
      //_allResRmsAllS[iter+1] = allResRmsFina;
      _dataResRmsInitS[iter] = dataResRmsInit;
      _dataResRmsFinaS[iter+1] = dataResRmsFina;
      _dataRes2NormSqInitS[iter] = dataRes2NormSqInit;
      _dataRes2NormSqFinaS[iter+1] = dataRes2NormSqFina;
      _bPenaltyRes2NormSqInitS[iter] = bPenaltyRes2NormSqInit;
      _bPenaltyRes2NormSqFinaS[iter+1] = bPenaltyRes2NormSqFina;
      _cPenaltyRes2NormSqInitS[iter] = cPenaltyRes2NormSqInit;
      _cPenaltyRes2NormSqFinaS[iter+1] = cPenaltyRes2NormSqFina;
      _allResRmsInitS[iter] = allResRmsInit;
      _allResRmsFinaS[iter+1] = allResRmsFina;
      _allRes2NormSqInitS[iter] = allRes2NormSqInit;
      _allRes2NormSqFinaS[iter+1] = allRes2NormSqFina;
      _stepLengthS[iter] = stepLength;
      _conditionNumberS[iter] = (float)x.cond();
      _b2NormS[iter] = twoNorm(b);
      _c2NormS[iter] = twoNorm(c);
      _deltaB2NormS[iter] = twoNorm(deltaB);
      _deltaC2NormS[iter] = twoNorm(deltaC);
      _gradient2NormS[iter] = twoNorm(v);
      _alphaBS[iter] = (float)_alphaB;
      _alphaCS[iter] = (float)_alphaC;
      _rmsPercentChangeS[iter] = rmsPercentChange;

      trace("Initial RMS of All Residuals: "+allResRmsInit);
      trace("Final RMS of All Residuals: "+allResRmsFina);
      trace("Initial RMS of Data Residuals: "+dataResRmsInit);
      trace("Final RMS of Data Residuals: "+dataResRmsFina);
      trace("Initial 2NormSq of All Residuals: "+allRes2NormSqInit);
      trace("Final 2NormSq of All Residuals: "+allRes2NormSqFina);
      trace("Initial 2NormSq of Data Residuals: "+dataRes2NormSqInit);
      trace("Final 2NormSq of Data Residuals: "+dataRes2NormSqFina);
      trace("Initial 2NormSq of b Penalty Residuals: "+bPenaltyRes2NormSqInit);
      trace("Final 2NormSq of b Penalty Residuals: "+bPenaltyRes2NormSqFina);
      trace("Initial 2NormSq of c Penalty Residuals: "+cPenaltyRes2NormSqInit);
      trace("Final 2NormSq of c Penalty Residuals: "+cPenaltyRes2NormSqFina);
      trace("Step Length: "+stepLength);
      trace("Condition Number: "+_conditionNumberS[iter]);
      trace("RMS Percent Change: "+rmsPercentChange);

      //printB(bNew);
      //printC(cNew);
      //printDeltaB(deltab);
      //printDeltaC(deltac);

      //Update c and b to newest c and b for the next iteration.
      b = copy(bNew);
      c = copy(cNew);
      _lastIter = iter;

      //Store c and b corresponding to smallest RMS of all residuals so far.
      if (allResRmsFina<minAllResRmsFina) {
        minIter = iter;
        minAllResRmsFina = allResRmsFina;
        minC = c;
        minB = b;
      }
      if (abs(rmsPercentChange)<_minRmsPercentChange) {
        trace("iter = "+iter);
        trace("iter+1 = "+(iter+1));
        _allResRmsAllS[iter] = allResRmsInit;
        _allResRmsAllS[iter+1] = allResRmsFina;
        _dataResRmsInitS[iter] = dataResRmsInit;
        _dataResRmsFinaS[iter+1] = dataResRmsFina;
        _dataRes2NormSqInitS[iter] = dataRes2NormSqInit;
        trace("link _dataRes2NormSqInitS = "+_dataRes2NormSqInitS[iter]);
        _dataRes2NormSqFinaS[iter+1] = dataRes2NormSqFina;
        _bPenaltyRes2NormSqInitS[iter] = bPenaltyRes2NormSqInit;
        _bPenaltyRes2NormSqFinaS[iter+1] = bPenaltyRes2NormSqFina;
        _cPenaltyRes2NormSqInitS[iter] = cPenaltyRes2NormSqInit;
        _cPenaltyRes2NormSqFinaS[iter+1] = cPenaltyRes2NormSqFina;
        _allResRmsInitS[iter] = allResRmsInit;
        _allResRmsFinaS[iter+1] = allResRmsFina;
        _allRes2NormSqInitS[iter] = allRes2NormSqInit;
        _allRes2NormSqFinaS[iter+1] = allRes2NormSqFina;
        _stepLengthS[iter] = stepLength;
        _conditionNumberS[iter] = (float)x.cond();
        _b2NormS[iter] = twoNorm(b);
        _c2NormS[iter] = twoNorm(c);
        _deltaB2NormS[iter] = twoNorm(deltaB);
        _deltaC2NormS[iter] = twoNorm(deltaC);
        _gradient2NormS[iter] = twoNorm(v);
        _alphaBS[iter] = (float)_alphaB;
        _alphaCS[iter] = (float)_alphaC;
        _rmsPercentChangeS[iter] = rmsPercentChange;
        ++iter;


        if (_end == true) {
          if (minAllResRmsFina<allResRmsFina) {
            trace("Reverting back to iteration "+minIter+"."); 
            trace("Produced the smallest Rms of all residuals");
            _lastIter = minIter;
            return new float[][]{minC,minB};
          }
          return new float[][]{c,b};
        }
        if (_penalize) {
          trace("###############iter = "+iter+"#####################");//db
          trace("Stopped Iterating: Percent change is below"+
          " maximum percent change required for stopping");
          trace("Maximum percent change before stopping = "+_minRmsPercentChange);
          trace("Percent change = "+rmsPercentChange);
          trace("Penalize on b, solve for non-zero alpha");
          trace("link alphaB = "+_alphaB);
          solveForAlphaB(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleW);
          //solveForAlphaC(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleZ);
          //solveForAlphaCAlphaB(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleZ);
          //return new float[][]{c,b};
          _penalize = false;
          trace("penalize = "+_penalize);
          trace("zelda alphaB = "+_alphaB);
        }
        else {
          trace("###############iter = "+iter+"#####################");//db
          trace("Eliminate penalization on b, alpha=0");
          _alphaB = 0.0;
          setWDiagAndZDiag(nc,kc,c,nb,kb,b);
          _end = true;

        }
      }
    }
    if (minAllResRmsFina<allResRmsFina) {
      trace("Reverting back to iteration "+minIter+"."); 
      trace("Produced the smallest Rms of all residuals");
      return new float[][]{minC,minB};
    }
    return new float[][]{c,b};
  }
  public float[][] getWaveletCInverseB(
    int nb, int kb, float[] bGuess, int nc, int kc, float[] cGuess,
    float[][] u, float[][] f, float[][] g, int niter)
  {
    Check.argument(-nb<kb,"-nb<kb");
    Check.argument(kb<=0,"kb<=0");
    Check.argument(-nc<kc,"-nc<kc");
    Check.argument(kc<=0,"kc<=0");

    //trace("bguess");//db
    //dump(bguess);//db
    //trace("cguess");//db
    //dump(cguess);//db

    _dataResRmsInitS = new float[niter];//db
    _dataResRmsFinaS = new float[niter];//db
    _dataRes2NormSqInitS = new float[niter];//db
    _dataRes2NormSqFinaS = new float[niter];//db
    _bPenaltyRes2NormSqInitS = new float[niter];//db
    _bPenaltyRes2NormSqFinaS = new float[niter];//db
    _cPenaltyRes2NormSqInitS = new float[niter];//db
    _cPenaltyRes2NormSqFinaS = new float[niter];//db
    _allResRmsInitS = new float[niter];//db
    _allResRmsFinaS = new float[niter];//db
    _allResRmsAllS = new float[niter];//db
    _allRes2NormSqInitS = new float[niter];//db
    _allRes2NormSqFinaS = new float[niter];//db
    _alphaBS = new float[niter];
    _alphaCS = new float[niter];

    _stepLengthS = new float[niter];//db
    _conditionNumberS = new float[niter];//db
    _b2NormS = new float[niter];//db
    _c2NormS = new float[niter];//db
    _deltaB2NormS = new float[niter];//db
    _deltaC2NormS = new float[niter];//db
    _gradient2NormS = new float[niter];//db
    _rmsPercentChangeS = new float[niter];//db

    _zDiag= new double[nc];
    _wDiag = new double[nb];
    _alphaB = 0.0;//015482059452917565f;
    _alphaC = 0.0;
    _scaleW = sqrt(dot(f,f));
    _scaleZ = sqrt(dot(f,f));

    int nt = f[0].length;
    int nx = f.length;
    int nbc = nb+nc;
    float[] b = copy(bGuess);
    float[] c = copy(cGuess);
    float[] bNew = new float[nb];
    float[] cNew = new float[nc];
    float[] deltaBDeltaC = new float[nbc];
    float[] deltaB = new float[nb];
    float[] deltaC = new float[nc];
    float[][] dataResInit = new float[nx][nt];
    float[] bPenaltyResInit = new float[nb];
    float[] cPenaltyResInit = new float[nc];
    float[][] dataResFina = new float[nx][nt];
    float[] bPenaltyResFina = new float[nb];
    float[] cPenaltyResFina = new float[nc];
    float[] minC = new float[nc];
    float[] minB = new float[nb];
    float minAllResRmsFina = Float.MAX_VALUE;
    int minIter = 0;
    float allResRmsInit = 0.0f;
    float allResRmsFina = 0.0f;
    float allRes2NormSqInit = 0.0f;
    float allRes2NormSqFina = 0.0f;
    float dataResRmsInit = 0.0f;
    float dataResRmsFina = 0.0f;
    float dataRes2NormSqInit = 0.0f;
    float dataRes2NormSqFina = 0.0f;
    float bPenaltyRes2NormSqInit = 0.0f;
    float bPenaltyRes2NormSqFina = 0.0f;
    float cPenaltyRes2NormSqInit = 0.0f;
    float cPenaltyRes2NormSqFina = 0.0f;
    float stepLength = 0.0f;
    float b2Norm = 0.0f;
    float c2Norm = 0.0f;
    float gradient2Norm = 0.0f;
    float rmsPercentChange = 0.0f;
    DMatrix[] xv = new DMatrix[2];//Approximated Hessian (x) gradient (v)
    DMatrix v = new DMatrix(nbc,1);
    DMatrix x = new DMatrix(nbc,nbc);

    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
    //Gauss-Newton iterations
    for (int iter=0; iter<niter; ++iter) {
      trace("###############iter = "+iter+"#####################");//db
      trace("alphaB = "+_alphaB);
      trace("current b:");
      dump(b);
      //SimplePlot sp = new SimplePlot();
      //sp.addTitle("b = "+iter+" alphaB = "+_alphaB);
      //sp.asPoints(b);
      trace("current c:");
      dump(c);
      trace("W diag:");
      dump(_wDiag);
      trace("Z diag:");
      dump(_zDiag);
      //Compute Initial Residuals
      dataResInit = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
      bPenaltyResInit = computeBPenaltyResidual(convertFToD(b),_wDiag);
      cPenaltyResInit = computeCPenaltyResidual(convertFToD(c),_zDiag);
      //Compute Initial Measures of Residuals
      dataResRmsInit = rms(dataResInit);
      dataRes2NormSqInit = (float) dot(dataResInit,dataResInit);
      bPenaltyRes2NormSqInit = (float) dotAll(bPenaltyResInit,bPenaltyResInit);
      cPenaltyRes2NormSqInit = (float) dotAll(cPenaltyResInit,cPenaltyResInit);
      allResRmsInit = rmsOfObjectiveFunction(dataResInit,bPenaltyResInit,cPenaltyResInit);
      allRes2NormSqInit = dataRes2NormSqInit+bPenaltyRes2NormSqInit+cPenaltyRes2NormSqInit;
      trace("all2normsq = "+allRes2NormSqInit);

      if (iter==0)
        _firstAll2NormSq = allRes2NormSqInit;


      ////////////Build Hessian and gradient/////////////////
      xv = buildApproxHessianAndGradient(nc,kc,c,nb,kb,b,dataResInit,u,g);
      x = xv[0];
      v = xv[1];

      ///////////Add Penalization/////////////////////
      xv = penalize(_wDiag,_zDiag,convertFToD(b),convertFToD(c),x,v);
      x = xv[0];
      v = xv[1];
    
      //////////Constrain deltab0 to be 0////////////
      xv = constrainDeltaB0ToBe0(nc,nb,kb,x,v);
      x = xv[0];
      v = xv[1];

      //Solve for deltab and deltac
      deltaBDeltaC = solveForSearchDirection(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleW,_firstAll2NormSq);
      //deltaBDeltaC = solveForSearchDirection(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleZ);
      deltaB = copy(nb,0,deltaBDeltaC);
      deltaC = copy(nc,nb,deltaBDeltaC);

      //Solve for scaled deltab and deltac
      stepLength = solveForStepLength(nb,kb,b,deltaB,nc,kc,c,deltaC,u,f,g,_wDiag,_zDiag);
      deltaB = mul(stepLength,deltaB);
      deltaC = mul(stepLength,deltaC);
      bNew = add(b,deltaB);
      cNew = add(c,deltaC);

      //Compute Final Residuals
      dataResFina = computeDataResidual(nc,kc,cNew,nb,kb,bNew,u,f,g);
      bPenaltyResFina = computeBPenaltyResidual(convertFToD(bNew),_wDiag);
      cPenaltyResFina = computeCPenaltyResidual(convertFToD(cNew),_zDiag);
      //Compute Final Measures of Residuals
      dataResRmsFina = rms(dataResFina);
      dataRes2NormSqFina = (float) dot(dataResFina,dataResFina);
      bPenaltyRes2NormSqFina = (float) dotAll(bPenaltyResFina,bPenaltyResFina);
      cPenaltyRes2NormSqFina = (float) dotAll(cPenaltyResFina,cPenaltyResFina);
      allResRmsFina = rmsOfObjectiveFunction(dataResFina,bPenaltyResFina,cPenaltyResFina);
      allRes2NormSqFina = dataRes2NormSqFina+bPenaltyRes2NormSqFina+cPenaltyRes2NormSqFina;
      rmsPercentChange = percentChange(allResRmsFina,allResRmsInit);

      //Debugging Information
      if (iter == 0)
        _allResRmsAllS[0] = allResRmsInit;
      else
        _allResRmsAllS[iter] = allResRmsInit;
      _allResRmsAllS[iter] = allResRmsInit;
      _allResRmsAllS[iter+1] = allResRmsFina;
      _dataResRmsInitS[iter] = dataResRmsInit;
      _dataResRmsFinaS[iter+1] = dataResRmsFina;
      _dataRes2NormSqInitS[iter] = dataRes2NormSqInit;
      _dataRes2NormSqFinaS[iter+1] = dataRes2NormSqFina;
      _bPenaltyRes2NormSqInitS[iter] = bPenaltyRes2NormSqInit;
      _bPenaltyRes2NormSqFinaS[iter+1] = bPenaltyRes2NormSqFina;
      _cPenaltyRes2NormSqInitS[iter] = cPenaltyRes2NormSqInit;
      _cPenaltyRes2NormSqFinaS[iter+1] = cPenaltyRes2NormSqFina;
      _allResRmsInitS[iter] = allResRmsInit;
      _allResRmsFinaS[iter+1] = allResRmsFina;
      _allRes2NormSqInitS[iter] = allRes2NormSqInit;
      _allRes2NormSqFinaS[iter+1] = allRes2NormSqFina;
      _stepLengthS[iter] = stepLength;
      _conditionNumberS[iter] = (float)x.cond();
      _b2NormS[iter] = twoNorm(b);
      _c2NormS[iter] = twoNorm(c);
      _deltaB2NormS[iter] = twoNorm(deltaB);
      _deltaC2NormS[iter] = twoNorm(deltaC);
      _gradient2NormS[iter] = twoNorm(v);
      _alphaBS[iter] = (float)_alphaB;
      _alphaCS[iter] = (float)_alphaC;
      _rmsPercentChangeS[iter] = rmsPercentChange;

      trace("Initial RMS of All Residuals: "+allResRmsInit);
      trace("Final RMS of All Residuals: "+allResRmsFina);
      trace("Initial RMS of Data Residuals: "+dataResRmsInit);
      trace("Final RMS of Data Residuals: "+dataResRmsFina);
      trace("Initial 2NormSq of All Residuals: "+allRes2NormSqInit);
      trace("Final 2NormSq of All Residuals: "+allRes2NormSqFina);
      trace("Initial 2NormSq of Data Residuals: "+dataRes2NormSqInit);
      trace("Final 2NormSq of Data Residuals: "+dataRes2NormSqFina);
      trace("Initial 2NormSq of b Penalty Residuals: "+bPenaltyRes2NormSqInit);
      trace("Final 2NormSq of b Penalty Residuals: "+bPenaltyRes2NormSqFina);
      trace("Initial 2NormSq of c Penalty Residuals: "+cPenaltyRes2NormSqInit);
      trace("Final 2NormSq of c Penalty Residuals: "+cPenaltyRes2NormSqFina);
      trace("Step Length: "+stepLength);
      trace("Condition Number: "+_conditionNumberS[iter]);
      trace("RMS Percent Change: "+rmsPercentChange);

      //printB(bNew);
      //printC(cNew);
      //printDeltaB(deltab);
      //printDeltaC(deltac);

      //Update c and b to newest c and b for the next iteration.
      b = copy(bNew);
      c = copy(cNew);
      //SimplePlot sp = new SimplePlot();
      //sp.addTitle("b iter = "+iter);
      //sp.addPoints(b);
      //SimplePlot sp1 = new SimplePlot();
      //sp1.addTitle("c iter = "+iter);
      //sp1.addPoints(c);
      _lastIter = iter;

      //Store c and b corresponding to smallest RMS of all residuals so far.
      if (allResRmsFina<minAllResRmsFina) {
        minIter = iter;
        minAllResRmsFina = allResRmsFina;
        minC = c;
        minB = b;
      }
      if (abs(rmsPercentChange)<_minRmsPercentChange) {
        _allResRmsAllS[iter] = allResRmsInit;
        _allResRmsAllS[iter+1] = allResRmsFina;
        _dataResRmsInitS[iter] = dataResRmsInit;
        _dataResRmsFinaS[iter+1] = dataResRmsFina;
        _dataRes2NormSqInitS[iter] = dataRes2NormSqInit;
        _dataRes2NormSqFinaS[iter+1] = dataRes2NormSqFina;
        _bPenaltyRes2NormSqInitS[iter] = bPenaltyRes2NormSqInit;
        _bPenaltyRes2NormSqFinaS[iter+1] = bPenaltyRes2NormSqFina;
        _cPenaltyRes2NormSqInitS[iter] = cPenaltyRes2NormSqInit;
        _cPenaltyRes2NormSqFinaS[iter+1] = cPenaltyRes2NormSqFina;
        _allResRmsInitS[iter] = allResRmsInit;
        _allResRmsFinaS[iter+1] = allResRmsFina;
        _allRes2NormSqInitS[iter] = allRes2NormSqInit;
        _allRes2NormSqFinaS[iter+1] = allRes2NormSqFina;
        _stepLengthS[iter] = stepLength;
        _conditionNumberS[iter] = (float)x.cond();
        _b2NormS[iter] = twoNorm(b);
        _c2NormS[iter] = twoNorm(c);
        _deltaB2NormS[iter] = twoNorm(deltaB);
        _deltaC2NormS[iter] = twoNorm(deltaC);
        _gradient2NormS[iter] = twoNorm(v);
        _alphaBS[iter] = (float)_alphaB;
        _alphaCS[iter] = (float)_alphaC;
        _rmsPercentChangeS[iter] = rmsPercentChange;
        ++iter;
        if (_end == true) {
          if (minAllResRmsFina<allResRmsFina) {
            trace("Reverting back to iteration "+minIter+"."); 
            trace("Produced the smallest Rms of all residuals");
            trace("minC");
            _lastIter = minIter;
            dump(minC);
            trace("minB");
            dump(minB);
            return new float[][]{minC,minB};
          }
          trace("lastC");
          dump(c);
          trace("lastB");
          dump(b);
          return new float[][]{c,b};
        }
        if (_penalize) {
          trace("###############iter = "+iter+"#####################");//db
          trace("Stopped Iterating: Percent change is below"+
          " maximum percent change required for stopping");
          trace("Maximum percent change before stopping = "+_minRmsPercentChange);
          trace("Percent change = "+rmsPercentChange);
          trace("Penalize on b, solve for non-zero alpha");
          solveForAlphaB(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleW,_firstAll2NormSq);
          //solveForAlphaC(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleZ);
          //solveForAlphaCAlphaB(nb,kb,b,nc,kc,c,u,f,g,x,v,_scaleZ);
          //return new float[][]{c,b};
          _penalize = false;
          trace("penalize = "+_penalize);
        }
        else {
          trace("###############iter = "+iter+"#####################");//
          trace("Eliminate penalization on b, alpha=0");
          _alphaB = 0.0;
          setWDiagAndZDiag(nc,kc,c,nb,kb,b);
          _end = true;

        }
      }
    }
    if (minAllResRmsFina<allResRmsFina) {
      trace("Reverting back to iteration "+minIter+"."); 
      trace("Produced the smallest Rms of all residuals");
      trace("minC");
      dump(minC);
      trace("minB");
      dump(minB);
      return new float[][]{minC,minB};
    }
    trace("lastC");
    dump(c);
    trace("lastB");
    dump(b);

    return new float[][]{c,b};
  }

  public int getLastIter() {
    trace("last iter = "+_lastIter);
    return _lastIter;

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
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index for b[0].
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

  public float[] getDataResRmsInitialS() {
    return _dataResRmsInitS;
  }
  public float[] getDataResRmsFinalS() {
    return _dataResRmsFinaS;
  }

  public float[] getDataRes2NormSqInitialS() {
    return _dataRes2NormSqInitS;
  }
  public float[] getDataRes2NormSqFinalS() {
    return _dataRes2NormSqFinaS;
  }

  public float[] getBPenaltyRes2NormSqInitialS() {
    return _bPenaltyRes2NormSqInitS;
  }
  public float[] getBPenaltyRes2NormSqFinalS() {
    return _bPenaltyRes2NormSqFinaS;
  }

  public float[] getCPenaltyRes2NormSqInitialS() {
    return _cPenaltyRes2NormSqInitS;
  }
  public float[] getCPenaltyRes2NormSqFinalS() {
    return _cPenaltyRes2NormSqFinaS;
  }

  public float[] getAllResRmsInitialS() {
    return _allResRmsInitS;
  }
  public float[] getAllResRmsFinalS() {
    return _allResRmsFinaS;
  }
  public float[] getAllResRmsAllS() {
    dump(_allResRmsAllS);
    return _allResRmsAllS;
  }

  public float[] getAllRes2NormSqInitialS() {
    return _allRes2NormSqInitS;
  }
  public float[] getAllRes2NormSqFinalS() {
    return _allRes2NormSqFinaS;
  }

  public float[] getStepLengthS() {
    return _stepLengthS;
  }

  public float[] getConditionNumberS() {
    return _conditionNumberS;
  }

  public float[] getB2NormS() {
    return _b2NormS;
  }

  public float[] getC2NormS() {
    return _c2NormS;
  }

  public float[] getDeltaB2NormS() {
    return _deltaB2NormS;
  }

  public float[] getDeltaC2NormS() {
    return _deltaC2NormS;
  }

  public float[] getGradient2NormS() {
    return _gradient2NormS;
  }

  public float[] getRmsPercentChangeS() {
    return _rmsPercentChangeS;
  }

  public float[] getAlphaB() {
    return _alphaBS;
  }
  public float[] getAlphaC() {
    return _alphaCS;
  }
  
   /**
   * Applies the specified wavelet C.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index for c[0].
   * @param c array of coefficients for the wavelet c.
   * @param f array with input sequence f(t).
   * @return array with filtered output sequence.
   */
  public float[] applyC(int nc, int kc, float[] c, float[] f) {
    return convolve(nc,kc,c,f);
  }
  public float[][] applyC(int nc, int kc, float[] c, float[][] f) {
    return convolve(nc,kc,c,f);
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


    
  /**
   * Returns the rms value of the image/trace.
   * The range of time samples included in 
   * the rms calculation is set by the setTimeRange(itmin,itmax)
   * method or by simply inputting your range into 
   * the appropriate method.
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

  /**
   * The overall residual contains a data residual and a 
   * penalization residual. To compute the rms of the total
   * residual, the individual residuals are needed.
   * @param dataResiduals data residuals 
   * @param bPenaltyResiduals penalization residual corresponding to b
   * @param cPenaltyResiduals penalization residual corresponding to c
   */
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
  
///////////////////////////////////////////////////////////////////////////
  // private
  private double _alphaB = 0.00;
  private double _alphaC = 0.00;
  private double _scaleW;
  private double _scaleZ;
  private double _minStepLength;
  private double _firstAll2NormSq;
  private double _firstAll2NormSqTemp;
  private double _alphaBManualValue;
  private float _minRmsPercentChange = 0.000f;
  private int _itmin = -1;
  private int _itmax = -1;
  private int _ng0 = 0;//staring size of array/image.
  private int _lastIter = 0;
  private float[] _dataResRmsInitS; 
  private float[] _dataResRmsFinaS; 
  private float[] _dataRes2NormSqInitS;
  private float[] _dataRes2NormSqFinaS;
  private float[] _bPenaltyRes2NormSqInitS;
  private float[] _bPenaltyRes2NormSqFinaS;
  private float[] _cPenaltyRes2NormSqInitS;
  private float[] _cPenaltyRes2NormSqFinaS;
  private float[] _allResRmsInitS;
  private float[] _allResRmsFinaS;
  private float[] _allResRmsAllS;
  private float[] _allRes2NormSqInitS;
  private float[] _allRes2NormSqFinaS;
  private float[] _stepLengthS; 
  private float[] _conditionNumberS;
  private float[] _b2NormS;
  private float[] _c2NormS;
  private float[] _deltaB2NormS;
  private float[] _deltaC2NormS;
  private float[] _gradient2NormS;
  private float[] _alphaBS;
  private float[] _alphaCS;
  private float[] _rmsPercentChangeS;
  private double[] _zDiag;
  private double[] _wDiag;
  private boolean _parallel = false;
  private boolean _penalize = true;
  private boolean _end = false;
  private boolean _setAlphaB = false;

  /**
   * Given an approximated Hessian and a gradient corresponding 
   * to previous estimates of c and b, we can compute the 
   * search direction for the Gauss-Newton method.
   * Contains no built in capabilities to catch the case where
   * Cholesky decomposition fails.
   * @param approximatedHessian 
   * @param gradient
   */
  private float[] solveForSearchDirection(
    DMatrix approximatedHessian, DMatrix gradient)
  {
    //solve SPD linear system of equations for deltab and deltac.
    DMatrixChd chd = new DMatrixChd(approximatedHessian);
    DMatrix deltaBDeltaC = chd.solve(gradient);
    return convertDToF(deltaBDeltaC.getArray());
  }

  /**
   * Given an approximated Hessian and a gradient corresponding 
   * to previous estimates of c and b, we can compute the 
   * search direction for the Gauss-Newton method. If
   * Cholesky decomposition fails, a new search direction
   * is computed based on a stabilized system of equations. 
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index of b[0].
   * @param b array of coefficients for the previous estimate of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index of c[0].
   * @param c array of coefficients for the previous estimate of the wavelet c.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   * @param approximatedHessian 
   * @param gradient
   * @param scaleW scale factor used to make stabilization effective.
   */
  private float[] solveForSearchDirection(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[] u, float[] f, float[] g,
                                          DMatrix approximatedHessian, DMatrix gradient, double scaleW) 
  {
    //solve SPD linear system of equations for deltab and deltac.
    DMatrixChd chd = new DMatrixChd(approximatedHessian);
    DMatrix deltaBDeltaC = new DMatrix(nb+nc,1);
    try {
      deltaBDeltaC = chd.solve(gradient);
    }
    catch (IllegalStateException e) {
      trace("caught exception");
      float preConditionNumber = (float) approximatedHessian.cond();
      trace("Condition Number Before Stabilization"+preConditionNumber);
      return solveForSearchDirectionAndAlphaB(nb,kb,b,nc,kc,c,u,f,g,approximatedHessian,gradient,scaleW);
      //return solveForSearchDirectionAndAlphaC(nb,kb,b,nc,kc,c,u,f,g,approximatedHessian,gradient,scaleW);
      //return solveForSearchDirectionAndAlphaCAlphaB(nb,kb,b,nc,kc,c,u,f,g,approximatedHessian,gradient,scaleW);
    }
    return convertDToF(deltaBDeltaC.getArray());
  }
  private float[] solveForSearchDirection(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[][] u, float[][] f, float[][] g,
                                          DMatrix approximatedHessian, DMatrix gradient, 
                                          double scaleW, double firstAll2NormSq) 
  {
    //solve SPD linear system of equations for deltab and deltac.
    DMatrixChd chd = new DMatrixChd(approximatedHessian);
    DMatrix deltaBDeltaC = new DMatrix(nb+nc,1);
    try {
      deltaBDeltaC = chd.solve(gradient);
    }
    catch (IllegalStateException e) {
      trace("caught exception");
      float preConditionNumber = (float) approximatedHessian.cond();
      trace("Condition Number Before Stabilization"+preConditionNumber);
      return solveForSearchDirectionAndAlphaB(nb,kb,b,nc,kc,c,u,f,g,approximatedHessian,gradient,scaleW,firstAll2NormSq);
      //return solveForSearchDirectionAndAlphaC(nb,kb,b,nc,kc,c,u,f,g,approximatedHessian,gradient,scaleW);
      //return solveForSearchDirectionAndAlphaCAlphaB(nb,kb,b,nc,kc,c,u,f,g,approximatedHessian,gradient,scaleW);
    }
    return convertDToF(deltaBDeltaC.getArray());
  }

  /**
   * Given an approximated Hessian and a gradient corresponding 
   * to previous estimates of c and b, we can compute the 
   * search direction for the Gauss-Newton method. If
   * Cholesky decomposition fails, a new search direction
   * is computed based on a stabilized system of equations. 
   * This stabilization is based on penalizing outer lags in b.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index of b[0].
   * @param b array of coefficients for the previous estimate of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index of c[0].
   * @param c array of coefficients for the previous estimate of the wavelet c.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   * @param approximatedHessian 
   * @param gradient
   * @param scaleW scale factor used to make stabilization effective.
   */
  private float[] solveForSearchDirectionAndAlphaB(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[] u, float[] f, float[] g,
                                          DMatrix approximatedHessian, DMatrix gradient, double scaleW) 
  {
    ObjectiveFunctionWithAlphaB objFuncWithAlphaB = new ObjectiveFunctionWithAlphaB();
    objFuncWithAlphaB.setParameters(nc,kc,c,nb,kb,b,u,f,g,approximatedHessian,gradient,scaleW);
    double maxAlphaB = objFuncWithAlphaB.getMaxAlphaB();
    double minAlphaB = objFuncWithAlphaB.getMinAlphaB();
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaB);
    _alphaB = bf.findMin(minAlphaB,maxAlphaB,minAlphaB);//justify these parameters.
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
    return objFuncWithAlphaB.getDeltaBDeltaC();
  }
  private float[] solveForSearchDirectionAndAlphaB(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[][] u, float[][] f, float[][] g,
                                          DMatrix approximatedHessian, DMatrix gradient, double scaleW, double firstAll2NormSq) 
  {
    ObjectiveFunctionWithAlphaB2D objFuncWithAlphaB = new ObjectiveFunctionWithAlphaB2D();
    objFuncWithAlphaB.setParameters(nc,kc,c,nb,kb,b,u,f,g,approximatedHessian,gradient,scaleW,firstAll2NormSq);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaB);
    _alphaB = bf.findMin(_alphaB,1.0,0.0001);//justify these parameters.
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
    return objFuncWithAlphaB.getDeltaBDeltaC();
  }

  /**
   * Experimental *
   * Given an approximated Hessian and a gradient corresponding 
   * to previous estimates of c and b, we can compute the 
   * search direction for the Gauss-Newton method. If
   * Cholesky decomposition fails, a new search direction
   * is computed based on a stabilized system of equations. 
   * This stabilization is based on penalizing outer lags in c.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index of b[0].
   * @param b array of coefficients for the previous estimate of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index of c[0].
   * @param c array of coefficients for the previous estimate of the wavelet c.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   * @param approximatedHessian 
   * @param gradient
   * @param scaleZ scale factor used to make stabilization effective.
   */
   /*
  private float[] solveForSearchDirectionAndAlphaC(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[] u, float[] f, float[] g,
                                          DMatrix approximatedHessian, DMatrix gradient, double scaleZ) 
  {
    ObjectiveFunctionWithAlphaC objFuncWithAlphaC = new ObjectiveFunctionWithAlphaC();
    objFuncWithAlphaC.setParameters(nc,kc,c,nb,kb,b,u,f,g,approximatedHessian,gradient,scaleZ);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaC);
    _alphaC = bf.findMin(_alphaC,1.0,0.0001);//justify these parameters.
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
    return objFuncWithAlphaC.getDeltaBDeltaC();
  }
  */
  /*private float[] solveForSearchDirectionAndAlphaC(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[][] u, float[][] f, float[][] g,
                                          DMatrix x, DMatrix v, double scaleZ) 
  {
    ObjectiveFunctionWithAlphaC objFuncWithAlphaC = new ObjectiveFunctionWithAlphaC();
    objFuncWithAlphaC.setParameters(nc,kc,c,nb,kb,b,u,f,g,x,v,scaleZ);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaC);
    _alphaC = bf.findMin(_alphaC,0.1,0.0001);//justify these parameters.
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
    return objFuncWithAlphaC.getDeltaBDeltaC();
  }*/

  /**
   * Experimental *
   * Given an approximated Hessian and a gradient corresponding 
   * to previous estimates of c and b, we can compute the 
   * search direction for the Gauss-Newton method. If
   * Cholesky decomposition fails, a new search direction
   * is computed based on a stabilized system of equations. 
   * This stabilization is based on penalizing outer lags in b and c.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index of b[0].
   * @param b array of coefficients for the previous estimate of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index of c[0].
   * @param c array of coefficients for the previous estimate of the wavelet c.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   * @param approximatedHessian 
   * @param gradient
   * @param scaleWandZ scale factor used to make stabilization effective.
   */
   /*
  private float[] solveForSearchDirectionAndAlphaCAlphaB(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[] u, float[] f, float[] g,
                                          DMatrix approximatedHessian, DMatrix gradient, double scaleWandZ) 
  {
    ObjectiveFunctionWithAlphaCAlphaB objFuncWithAlphaCAlphaB = new ObjectiveFunctionWithAlphaCAlphaB();
    objFuncWithAlphaCAlphaB.setParameters(nc,kc,c,nb,kb,b,u,f,g,approximatedHessian,gradient,scaleWandZ);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaCAlphaB);
    _alphaC = bf.findMin(_alphaC,1.0,0.0001);//justify these parameters.
    _alphaB = _alphaC;
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
    return objFuncWithAlphaCAlphaB.getDeltaBDeltaC();
  }
  */
 /*private float[] solveForSearchDirectionAndAlphaCAlphaB(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[][] u, float[][] f, float[][] g,
                                          DMatrix x, DMatrix v, double scaleWandZ) 
  {
    ObjectiveFunctionWithAlphaCAlphaB objFuncWithAlphaCAlphaB = new ObjectiveFunctionWithAlphaCAlphaB();
    objFuncWithAlphaCAlphaB.setParameters(nc,kc,c,nb,kb,b,u,f,g,x,v,scaleWandZ);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaCAlphaB);
    _alphaC = bf.findMin(_alphaC,0.1,0.0001);//justify these parameters.
    _alphaB = _alphaC;
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
    return objFuncWithAlphaCAlphaB.getDeltaBDeltaC();
  }*/



  /**
   * Solves for the degree at which to penalize outer lags of b.
   * @param nb number of samples in the inverse wavelet b.
   * @param kb the sample index of b[0].
   * @param b array of coefficients for the previous estimate of the inverse wavelet b.
   * @param nc number of samples in the wavelet c.
   * @param kc the sample index of c[0].
   * @param c array of coefficients for the previous estimate of the wavelet c.
   * @param u relates PP time to PS time (in samples).
   * @param f the PP trace.
   * @param g the PS trace.
   * @param approximatedHessian 
   * @param gradient
   * @param scaleWandZ scale factor used to make stabilization effective.
   */
  private void solveForAlphaB(int nb, int kb, float[] b,
                                 int nc, int kc, float[] c,
                                 float[] u, float[] f, float[] g,
                                 DMatrix approximatedHessian, DMatrix gradient, double scaleW) 
  {
    float[] dataResInit = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
    ////////////Build Hessian and gradient/////////////////
    DMatrix[] approxHessAndGrad = buildApproxHessianAndGradient(nc,kc,c,nb,kb,b,dataResInit,u,g);
    approximatedHessian = approxHessAndGrad[0];
    gradient = approxHessAndGrad[1];

    ///////////Add Penalization/////////////////////
    approxHessAndGrad = penalize(_wDiag,_zDiag,convertFToD(b),convertFToD(c),approximatedHessian,gradient);
    approximatedHessian = approxHessAndGrad[0];
    gradient = approxHessAndGrad[1];
  
    //////////Constrain deltab0 to be 0////////////
    approxHessAndGrad = constrainDeltaB0ToBe0(nc,nb,kb,approximatedHessian,gradient);
    approximatedHessian = approxHessAndGrad[0];
    gradient = approxHessAndGrad[1];
    trace("b");
    dump(b);
    trace("alphaB = "+_alphaB);
    ObjectiveFunctionWithAlphaB objFuncWithAlphaB = new ObjectiveFunctionWithAlphaB();
    objFuncWithAlphaB.setParameters(nc,kc,c,nb,kb,b,u,f,g,approximatedHessian,gradient,scaleW);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaB);
    double maxAlphaB = objFuncWithAlphaB.getMaxAlphaB();
    double minAlphaB = objFuncWithAlphaB.getMinAlphaB();
    trace("maxAlphaBTest = "+maxAlphaB);
    trace("minAlphaBTest = "+minAlphaB);
    if (_setAlphaB==true)
      _alphaB = _alphaBManualValue;
    else
      _alphaB = bf.findMin(minAlphaB,maxAlphaB,minAlphaB);//justify these parameters.
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);//updates parameters that are affected by alphaB
  }
  private void solveForAlphaB(int nb, int kb, float[] b,
                                 int nc, int kc, float[] c,
                                 float[][] u, float[][] f, float[][] g,
                                 DMatrix approximatedHessian, DMatrix gradient, double scaleW, double firstAll2NormSq) 
  {
    float[][] dataResInit = computeDataResidual(nc,kc,c,nb,kb,b,u,f,g);
    ////////////Build Hessian and gradient/////////////////
    DMatrix[] approxHessAndGrad = buildApproxHessianAndGradient(nc,kc,c,nb,kb,b,dataResInit,u,g);
    approximatedHessian = approxHessAndGrad[0];
    gradient = approxHessAndGrad[1];

    ///////////Add Penalization/////////////////////
    approxHessAndGrad = penalize(_wDiag,_zDiag,convertFToD(b),convertFToD(c),approximatedHessian,gradient);
    approximatedHessian = approxHessAndGrad[0];
    gradient = approxHessAndGrad[1];
  
    //////////Constrain deltab0 to be 0////////////
    approxHessAndGrad = constrainDeltaB0ToBe0(nc,nb,kb,approximatedHessian,gradient);
    approximatedHessian = approxHessAndGrad[0];
    gradient = approxHessAndGrad[1];
    ObjectiveFunctionWithAlphaB2D objFuncWithAlphaB = new ObjectiveFunctionWithAlphaB2D();
    objFuncWithAlphaB.setParameters(nc,kc,c,nb,kb,b,u,f,g,approximatedHessian,gradient,scaleW,firstAll2NormSq);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaB);
    double maxAlphaB = objFuncWithAlphaB.getMaxAlphaB();
    double minAlphaB = objFuncWithAlphaB.getMinAlphaB();
    trace("maxAlphaBTest = "+maxAlphaB);
    trace("minAlphaBTest = "+minAlphaB);
    if (_setAlphaB==true)
      _alphaB = _alphaBManualValue;
    else
      _alphaB = bf.findMin(minAlphaB,maxAlphaB,minAlphaB);//justify these parameters.
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);//updates parameters that are affected by alphaB
  }
  public void setAlphaB(double alphaB) {
    _setAlphaB=true;
    _alphaBManualValue=alphaB;
  }

  /*private void solveForAlphaC(int nb, int kb, float[] b,
                                 int nc, int kc, float[] c,
                                 float[] u, float[] f, float[] g,
                                 DMatrix x, DMatrix v, double scaleZ) 
  {
    ObjectiveFunctionWithAlphaC objFuncWithAlphaC = new ObjectiveFunctionWithAlphaC();
    objFuncWithAlphaC.setParameters(nc,kc,c,nb,kb,b,u,f,g,x,v,scaleZ);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaC);
    _alphaC = bf.findMin(_alphaC,0.1,_alphaC);//justify these parameters.
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);//updates parameters that are affected by alphaC
  }*/
  /*private void solveForAlphaC(int nb, int kb, float[] b,
                                 int nc, int kc, float[] c,
                                 float[][] u, float[][] f, float[][] g,
                                 DMatrix x, DMatrix v, double scaleZ) 
  {
    ObjectiveFunctionWithAlphaC objFuncWithAlphaC = new ObjectiveFunctionWithAlphaC();
    objFuncWithAlphaC.setParameters(nc,kc,c,nb,kb,b,u,f,g,x,v,scaleZ);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaC);
    _alphaC = bf.findMin(_alphaC,0.1,_alphaC);//justify these parameters.
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);//updates parameters that are affected by alphaC
  }

  private void solveForAlphaCAlphaB(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[] u, float[] f, float[] g,
                                          DMatrix x, DMatrix v, double scaleWandZ) 
  {
    ObjectiveFunctionWithAlphaCAlphaB objFuncWithAlphaCAlphaB = new ObjectiveFunctionWithAlphaCAlphaB();
    objFuncWithAlphaCAlphaB.setParameters(nc,kc,c,nb,kb,b,u,f,g,x,v,scaleWandZ);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaCAlphaB);
    //_alphaC = bf.findMin(_alphaC,0.001,0.01);//justify these parameters.
    _alphaB = bf.findMin(0.00001,1,0.00);//justify these parameters.
    _alphaB = _alphaC;
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
  }
  */
  /*private void solveForAlphaCAlphaB(int nb, int kb, float[] b,
                                          int nc, int kc, float[] c,
                                          float[][] u, float[][] f, float[][] g,
                                          DMatrix x, DMatrix v, double scaleWandZ) 
  {
    ObjectiveFunctionWithAlphaCAlphaB objFuncWithAlphaCAlphaB = new ObjectiveFunctionWithAlphaCAlphaB();
    objFuncWithAlphaCAlphaB.setParameters(nc,kc,c,nb,kb,b,u,f,g,x,v,scaleWandZ);
    BrentMinFinder bf = new BrentMinFinder(objFuncWithAlphaCAlphaB);
    _alphaC = bf.findMin(_alphaC,0.1,0.0001);//justify these parameters.
    _alphaB = _alphaC;
    setWDiagAndZDiag(nc,kc,c,nb,kb,b);
  }*/


  private float solveForStepLength(int nb, int kb, float[] b, float[] deltaB,
    int nc, int kc, float[] c, float[] deltaC,
    float[] u, float[] f, float[] g,
    double[] wDiag, double[] zDiag)
  {
    ObjectiveFunction of = new ObjectiveFunction();
    of.setParameters(nc,kc,c,deltaC,nb,kb,b,deltaB,u,f,g,wDiag,zDiag);
    BrentMinFinder bf = new BrentMinFinder(of);
    return (float) bf.findMin(_minStepLength,1.0,_minStepLength);
  }
  private float solveForStepLength(int nb, int kb, float[] b, float[] deltaB,
    int nc, int kc, float[] c, float[] deltaC,
    float[][] u, float[][] f, float[][] g,
    double[] wDiag, double[] zDiag)
  {
    ObjectiveFunction2D of = new ObjectiveFunction2D();
    of.setParameters(nc,kc,c,deltaC,nb,kb,b,deltaB,u,f,g,wDiag,zDiag);
    BrentMinFinder bf = new BrentMinFinder(of);
    return (float) bf.findMin(_minStepLength,1.0,_minStepLength);
  }

  private float[][] separateDeltaBDeltaC(int nb, int nc, DMatrix deltaBDeltaC) {
    double[] tempBC = deltaBDeltaC.getArray();
    float[] deltaB = new float[nb];
    float[] deltaC = new float[nc];
    for (int ib=0; ib<nb; ++ib) {
      deltaB[ib] = (float) tempBC[ib];
    }
    for (int ic=0; ic<nc; ++ic) {
      deltaC[ic] = (float) tempBC[nb+ic];
    }
    return new float[][]{deltaB,deltaC};
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
  
  private void setWDiagAndZDiag(int nc, int kc, float[] c, int nb, int kb, 
    float[] b) 
  {
    trace("_alphaB = "+_alphaB);
    _wDiag = computeWDiagonal(nb,kb,_alphaB,_scaleW);
    _zDiag = computeZDiagonal(nc,kc,_alphaC,_scaleZ);
  }

  private double[] computeWDiagonal(int nb, int kb, double alphaB, double scaleW) {
    double[] wDiag = new double[nb];
    int lag = 0;
    for (int ib=0; ib<nb; ++ib) {
      lag = abs(-kb-ib);
      wDiag[ib] = sqrt(alphaB*scaleW)*lag;
    }
    return wDiag;
  }
  private double[] computeWDiagonalNoAlphaB(int nb, int kb, double scaleW) {
    double[] wDiag = new double[nb];
    int lag = 0;
    for (int ib=0; ib<nb; ++ib) {
      lag = abs(-kb-ib);
      wDiag[ib] = sqrt(scaleW)*lag;
    }
    return wDiag;
  }
  private double[] computeWTransposeWDiagonal(double[] wDiag) {
    int nb = wDiag.length;
    double[] wTransWDiag = new double[nb];
    for (int ib=0; ib<nb; ++ib) 
      wTransWDiag[ib] = wDiag[ib]*wDiag[ib];
    return wTransWDiag;
  }
  private double[] computeWTransposeWB(double[] wTransWDiag, double[] b) {
    int nb = b.length;
    double[] wTransWB = new double[nb];
    for (int ib=0; ib<nb; ++ib) {
      wTransWB[ib] = b[ib]*wTransWDiag[ib];
    }
    return wTransWB;
  }
  private double[] computeZTransposeZC(double[] zTransZDiag, double[] c) {
    int nc = c.length;
    double[] zTransZC = new double[nc];
    for (int ic=0; ic<nc; ++ic) {
      zTransZC[ic] = c[ic]*zTransZDiag[ic];
    }
    return zTransZC;
  }

  private double[] computeZDiagonal(int nc, int kc, double alphaC, double scaleZ) {
    double[] zDiag = new double[nc];
    int lag = 0;
    for (int ic=0; ic<nc; ++ic) {
      lag = abs(-kc-ic);
      zDiag[ic] = sqrt(alphaC*scaleZ)*lag;
    }
    return zDiag;
  }
  private double[] computeZTransposeZDiagonal(double[] zDiag) {
    int nc = zDiag.length;
    double[] zTransZDiag = new double[nc];
    for (int ic=0; ic<nc; ++ic) 
      zTransZDiag[ic] = zDiag[ic]*zDiag[ic];
    return zTransZDiag;
  }
  private DMatrix addWWDiagAndZZDiagToX(DMatrix x, double[] wTransWDiag, double[] zTransZDiag) {
    int nb = wTransWDiag.length;
    int nc = zTransZDiag.length;
    int nbc = nb+nc;
    DMatrix zero = new DMatrix(nbc,nbc);
    DMatrix xNew = (x.plus(zero));
    for (int ib=0; ib<nb; ++ib) {
      xNew.set(ib,ib,(x.get(ib,ib)+wTransWDiag[ib]));
    }
    for (int ic=0; ic<nc; ++ic) {
      xNew.set(nb+ic,nb+ic,(x.get(nb+ic,nb+ic)+zTransZDiag[ic]));
    }
    return xNew;
  }
  private DMatrix subtractWWbAndZZcFromV(DMatrix v, double[] wTransWB, double[] zTransZC) {
    int nb = wTransWB.length;
    int nc = zTransZC.length;
    int nbc = nb+nc;
    DMatrix zero = new DMatrix(nbc,1);
    DMatrix vNew = v.plus(zero);
    for (int ib=0; ib<nb; ++ib) {
      vNew.set(ib,0,(v.get(ib,0)-wTransWB[ib]));
    }
    for (int ic=0; ic<nc; ++ic) {
      vNew.set(nb+ic,0,(v.get(nb+ic,0)-zTransZC[ic]));
    }
    return vNew;
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
  private float[] computeBPenaltyResidual(double[] b, double[] wDiag) {
    int nb = b.length;
    float[] r = new float[nb];
    for (int ib=0; ib<nb; ++ib) {
      r[ib] = (float) (b[ib]*wDiag[ib]);
    }
    return r;
  }
  private float[] computeCPenaltyResidual(double[] c, double[] zDiag) {
    int nc = c.length;
    float[] r = new float[nc];
    for (int ic=0; ic<nc; ++ic) {
      r[ic] = (float) (c[ic]*zDiag[ic]);
    }
    return r;
  }

  /**
   * The objective function changes when outer lags of b are being penalized.
   * This class computes this objective function and also computes the 
   * optimal degree at which to penalize the outer lags of b.
   */
  private class ObjectiveFunctionWithAlphaB implements BrentMinFinder.Function {
    public void setParameters(int nc, int kc, float[] c,
                              int nb, int kb, float[] b,
                              float[] u, float[] f, float[] g, 
                              DMatrix x, DMatrix v,
                              double scaleW) 
    {
      _nb = nb;
      _nc = nc;
      _kc = kc;
      _kb = kb;
      _b = copy(b);
      _c = copy(c);
      _u = u;
      _f = f;
      _g = g;
      _x = x;
      _v = v;
      _scaleW = scaleW;
    }
    public double evaluate(double alphaB) {
      //compute penalizations
      trace("alphaBTest = "+alphaB);
      trace("scaleW = "+_scaleW);
      _wDiag = computeWDiagonal(_nb,_kb,alphaB,_scaleW);
      trace("b");
      dump(_b);
      _zDiag = new double[_nc];
      _wTransWDiag = computeWTransposeWDiagonal(_wDiag);//Diagonal of W'W
      _zTransZDiag = new double[_nc];
      _wTransWB = computeWTransposeWB(_wTransWDiag,convertFToD(_b));//W'Wb
      _zTransZC = new double[_nc];

      //add penalizations
      _xNew = addWWDiagAndZZDiagToX(_x,_wTransWDiag,_zTransZDiag);
      _vNew = subtractWWbAndZZcFromV(_v,_wTransWB,_zTransZC);

      //solve system
      DMatrixChd chd = new DMatrixChd(_xNew);
      _deltaBDeltaC = new float[_nb+_nc];
      try {
        _deltaBDeltaC = solveForSearchDirection(_xNew,_vNew);
      }      
      catch (IllegalStateException e) {
        return Float.MAX_VALUE;
      }
      _deltaB = copy(_nb,0,_deltaBDeltaC);
      _deltaC = copy(_nc,_nb,_deltaBDeltaC);
      trace("_deltaB = ");
      dump(_deltaB);
      trace("_deltaC = ");
      dump(_deltaC);

      //update deltab and deltac and c and b
      _stepLength = solveForStepLength(_nb,_kb,_b,_deltaB,_nc,_kc,_c,_deltaC,_u,_f,_g,_wDiag,_zDiag);
      _deltaB = mul(_stepLength,_deltaB);
      _deltaC = mul(_stepLength,_deltaC);
      _bNew = add(_b,_deltaB);
      _cNew = add(_c,_deltaC);
      
      //compute residuals
      _dataRes = computeDataResidual(_nc,_kc,_cNew,_nb,_kb,_bNew,_u,_f,_g);
      trace("dataRes");
      dump(_dataRes);
      _bPenaltyRes = computeBPenaltyResidual(convertFToD(_bNew),_wDiag);
      _cPenaltyRes = computeCPenaltyResidual(convertFToD(_cNew),_zDiag);
      _data2NormSq = dot(_dataRes,_dataRes);
      _bPenalty2NormSq = dotAll(_cPenaltyRes,_cPenaltyRes);
      _cPenalty2NormSq = dotAll(_bPenaltyRes,_bPenaltyRes);
      trace("data2NormSq = "+_data2NormSq);
      trace("bpenalty2normsq= "+_bPenalty2NormSq);
      trace("alphaB = "+alphaB);
      trace("dottotal = "+(_data2NormSq+_bPenalty2NormSq+_cPenalty2NormSq));
      trace("rmsrftotal = "+rmsOfObjectiveFunction(_dataRes,_bPenaltyRes,_cPenaltyRes));
      return _data2NormSq+_bPenalty2NormSq+_cPenalty2NormSq;
    }
    public double getMaxAlphaB() {
      trace("Max **************************");
      double eps = 1e-10;
      double[] wDiag = computeWDiagonalNoAlphaB(_nb,_kb,_scaleW);
      double[] zDiag = new double[_nc];
      double[] wTransWDiag = computeWTransposeWDiagonal(wDiag);//Diagonal of W'W
      double[] zTransZDiag = new double[_nc];
      double[] wTransWB = computeWTransposeWB(wTransWDiag,convertFToD(_b));//W'Wb
      double[] zTransZC = new double[_nc];

      float[] dataResInit = computeDataResidual(_nc,_kc,_c,_nb,_kb,_b,_u,_f,_g);
      float[] bPenaltyResInit = computeBPenaltyResidual(convertFToD(_b),wDiag);
      float[] cPenaltyResInit = computeCPenaltyResidual(convertFToD(_c),zDiag);
      double data2NormSqInit = dot(dataResInit,dataResInit);
      double bPenalty2NormSqInit = dotAll(bPenaltyResInit,bPenaltyResInit)+eps;
      double cPenalty2NormSqInit = dotAll(cPenaltyResInit,cPenaltyResInit);
      double sumInit = data2NormSqInit+bPenalty2NormSqInit+cPenalty2NormSqInit;
      trace("firstAll2NormSq = "+_firstAll2NormSq);
      trace("data2NormSqInit = "+data2NormSqInit);
      trace("bPenalty2NormSqInit = "+bPenalty2NormSqInit);
      trace("cPenalty2NormSqInit = "+cPenalty2NormSqInit);
      double partOfDifference = 0.25;
      double difference = _firstAll2NormSq-data2NormSqInit;
      trace("difference = "+difference);
      if (difference<0)
        difference = _firstAll2NormSqTemp-data2NormSqInit;
      trace("_firstAll2NormSqTemp = "+_firstAll2NormSqTemp);
      trace("difference = "+difference);

      double alpha = ((partOfDifference*difference+data2NormSqInit)-data2NormSqInit)/bPenalty2NormSqInit;
      trace("**************************");
      return alpha;
    }
    public double getMinAlphaB() {
      double eps = 1e-10;
      double[] wDiag = computeWDiagonalNoAlphaB(_nb,_kb,_scaleW);
      double[] zDiag = new double[_nc];
      double[] wTransWDiag = computeWTransposeWDiagonal(wDiag);//Diagonal of W'W
      double[] zTransZDiag = new double[_nc];
      double[] wTransWB = computeWTransposeWB(wTransWDiag,convertFToD(_b));//W'Wb
      double[] zTransZC = new double[_nc];
      trace("wDiag:");
      dump(wDiag);
      trace("b:");
      dump(_b);

      float[] dataResInit = computeDataResidual(_nc,_kc,_c,_nb,_kb,_b,_u,_f,_g);
      float[] bPenaltyResInit = computeBPenaltyResidual(convertFToD(_b),wDiag);
      trace("bPen");
      dump(bPenaltyResInit);
      float[] cPenaltyResInit = computeCPenaltyResidual(convertFToD(_c),zDiag);
      double data2NormSqInit = dot(dataResInit,dataResInit);
      double bPenalty2NormSqInit = dotAll(bPenaltyResInit,bPenaltyResInit)+eps;
      double cPenalty2NormSqInit = dotAll(cPenaltyResInit,cPenaltyResInit);
      double sumInit = data2NormSqInit+bPenalty2NormSqInit+cPenalty2NormSqInit;
      double partOfDifference = 0.001;
      double difference = _firstAll2NormSq-data2NormSqInit;
      trace("difference = "+difference);
      if (difference<0)
        difference = _firstAll2NormSqTemp-data2NormSqInit;
      trace("_firstAll2NormSqTemp = "+_firstAll2NormSqTemp);
      trace("difference = "+difference);
      trace("bepenalty = "+bPenalty2NormSqInit);
      double alpha = ((partOfDifference*difference+data2NormSqInit)-data2NormSqInit)/bPenalty2NormSqInit;
      trace("Hey listen!:"+alpha);

      return alpha;
    }

    public float[] getDeltaBDeltaC() {
      return _deltaBDeltaC;
    }

    private int _nb,_nc,_kb,_kc;
    private double _scaleW;
    private float _stepLength;
    private double _data2NormSq,_bPenalty2NormSq,_cPenalty2NormSq;
    private float[] _b,_bNew,_c,_cNew,_deltaB,_deltaC,_deltaBDeltaC;
    private float[] _dataRes,_bPenaltyRes,_cPenaltyRes;
    private float[] _u,_f,_g;
    private double[] _wDiag,_zDiag,_wTransWDiag,_zTransZDiag,_wTransWB,_zTransZC;
    private DMatrix _x,_v,_xNew,_vNew;
  }
  private class ObjectiveFunctionWithAlphaB2D implements BrentMinFinder.Function {
    public void setParameters(int nc, int kc, float[] c,
                              int nb, int kb, float[] b,
                              float[][] u, float[][] f, float[][] g, 
                              DMatrix x, DMatrix v,
                              double scaleW, double firstAll2NormSq) 
    {
      _nb = nb;
      _nc = nc;
      _kc = kc;
      _kb = kb;
      _b = copy(b);
      _c = copy(c);
      _u = u;
      _f = f;
      _g = g;
      _x = x;
      _v = v;
      _firstAll2NormSq = firstAll2NormSq;
      _scaleW = scaleW;
    }
    public double evaluate(double alphaB) {
      //compute penalizations
      _wDiag = computeWDiagonal(_nb,_kb,alphaB,_scaleW);
      _zDiag = new double[_nc];
      _wTransWDiag = computeWTransposeWDiagonal(_wDiag);//Diagonal of W'W
      _zTransZDiag = new double[_nc];
      _wTransWB = computeWTransposeWB(_wTransWDiag,convertFToD(_b));//W'Wb
      _zTransZC = new double[_nc];

      _dataResInit = computeDataResidual(_nc,_kc,_c,_nb,_kb,_b,_u,_f,_g);
      _bPenaltyResInit = computeBPenaltyResidual(convertFToD(_b),_wDiag);
      _cPenaltyResInit = computeCPenaltyResidual(convertFToD(_c),_zDiag);
      double _data2NormSqInit = dot(_dataResInit,_dataResInit);
      double _bPenalty2NormSqInit = dotAll(_bPenaltyResInit,_bPenaltyResInit);
      double _cPenalty2NormSqInit = dotAll(_cPenaltyResInit,_cPenaltyResInit);
      double sumInit = _data2NormSqInit+_bPenalty2NormSqInit+_cPenalty2NormSqInit;

      //add penalizations
      _xNew = addWWDiagAndZZDiagToX(_x,_wTransWDiag,_zTransZDiag);
      _vNew = subtractWWbAndZZcFromV(_v,_wTransWB,_zTransZC);
      trace("old condition number = "+_x.cond());
      trace("new condition number = "+_xNew.cond());
      printDiagonal(_x, "old X");
      printDiagonal(_xNew, "New X");
      trace("old v");
      trace(_v.toString());
      trace("new v");
      trace(_vNew.toString());

      //solve system
      DMatrixChd chd = new DMatrixChd(_xNew);
      float[] deltaBDeltaC = new float[_nb+_nc];
      try {
        _deltaBDeltaC = solveForSearchDirection(_xNew,_vNew);
      }      
      catch (IllegalStateException e) {
        return Float.MAX_VALUE;
      }
      _deltaB = copy(_nb,0,_deltaBDeltaC);
      _deltaC = copy(_nc,_nb,_deltaBDeltaC);
      trace("_deltaB");
      dump(_deltaB);
      trace("_deltaC");
      dump(_deltaC);

      //update deltab and deltac and c and b
      _stepLength = solveForStepLength(_nb,_kb,_b,_deltaB,_nc,_kc,_c,_deltaC,_u,_f,_g,_wDiag,_zDiag);
      trace("step length = "+_stepLength);
      _deltaB = mul(_stepLength,_deltaB);
      _deltaC = mul(_stepLength,_deltaC);
      trace("_deltaB");
      dump(_deltaB);
      trace("_deltaC");
      dump(_deltaC);
      _bNew = add(_b,_deltaB);
      _cNew = add(_c,_deltaC);
      trace("_b");
      dump(_b);
      trace("_c");
      dump(_c);
      trace("_bNew");
      dump(_bNew);
      trace("_cNew");
      dump(_cNew);
      
      //compute residuals
      _dataRes = computeDataResidual(_nc,_kc,_cNew,_nb,_kb,_bNew,_u,_f,_g);
      _bPenaltyRes = computeBPenaltyResidual(convertFToD(_bNew),_wDiag);
      _cPenaltyRes = computeCPenaltyResidual(convertFToD(_cNew),_zDiag);
      _data2NormSq = dot(_dataRes,_dataRes);
      _bPenalty2NormSq = dotAll(_bPenaltyRes,_bPenaltyRes);
      _cPenalty2NormSq = dotAll(_cPenaltyRes,_cPenaltyRes);
      double sum = _data2NormSq+_bPenalty2NormSq+_cPenalty2NormSq;
      trace("sumInit = "+sumInit);
      trace("firstAll2NormSq = "+_firstAll2NormSq);
      if (sumInit>_firstAll2NormSq) {
        trace("Too much penalization");
        return Float.MAX_VALUE;
      }

      trace("alphaB = "+alphaB);
      trace("dot pen = "+_bPenalty2NormSq);
      trace("dot data = "+_data2NormSq);
      trace("dot total = "+(_data2NormSq+_bPenalty2NormSq+_cPenalty2NormSq));
      trace("rmsrftotal = "+rmsOfObjectiveFunction(_dataRes,_bPenaltyRes,_cPenaltyRes));
      return sum;
    }
    public double getMaxAlphaB() {
      double eps = 1e-10;
      double[] wDiag = computeWDiagonalNoAlphaB(_nb,_kb,_scaleW);
      double[] zDiag = new double[_nc];
      double[] wTransWDiag = computeWTransposeWDiagonal(wDiag);//Diagonal of W'W
      double[] zTransZDiag = new double[_nc];
      double[] wTransWB = computeWTransposeWB(wTransWDiag,convertFToD(_b));//W'Wb
      double[] zTransZC = new double[_nc];

      float[][] dataResInit = computeDataResidual(_nc,_kc,_c,_nb,_kb,_b,_u,_f,_g);
      float[] bPenaltyResInit = computeBPenaltyResidual(convertFToD(_b),wDiag);
      float[] cPenaltyResInit = computeCPenaltyResidual(convertFToD(_c),zDiag);
      double data2NormSqInit = dot(dataResInit,dataResInit);
      double bPenalty2NormSqInit = dotAll(bPenaltyResInit,bPenaltyResInit)+eps;
      double cPenalty2NormSqInit = dotAll(cPenaltyResInit,cPenaltyResInit);
      double sumInit = data2NormSqInit+bPenalty2NormSqInit+cPenalty2NormSqInit;
      trace("firstAll2NormSq = "+_firstAll2NormSq);
      trace("data2NormSqInit = "+data2NormSqInit);
      trace("bPenalty2NormSqInit = "+bPenalty2NormSqInit);
      trace("cPenalty2NormSqInit = "+cPenalty2NormSqInit);
      double partOfDifference = 0.25;
      double difference = _firstAll2NormSq-data2NormSqInit+1e-10;
      double alpha = ((partOfDifference*difference+data2NormSqInit)-data2NormSqInit)/bPenalty2NormSqInit;
      return (_firstAll2NormSq-data2NormSqInit)/(2.0*(bPenalty2NormSqInit+cPenalty2NormSqInit));
    }
    public double getMinAlphaB() {
      double eps = 1e-10;
      double[] wDiag = computeWDiagonalNoAlphaB(_nb,_kb,_scaleW);
      double[] zDiag = new double[_nc];
      double[] wTransWDiag = computeWTransposeWDiagonal(wDiag);//Diagonal of W'W
      double[] zTransZDiag = new double[_nc];
      double[] wTransWB = computeWTransposeWB(wTransWDiag,convertFToD(_b));//W'Wb
      double[] zTransZC = new double[_nc];
      trace("wDiag:");
      dump(wDiag);
      trace("b:");
      dump(_b);

      float[][] dataResInit = computeDataResidual(_nc,_kc,_c,_nb,_kb,_b,_u,_f,_g);
      float[] bPenaltyResInit = computeBPenaltyResidual(convertFToD(_b),wDiag);
      trace("bPen");
      dump(bPenaltyResInit);
      float[] cPenaltyResInit = computeCPenaltyResidual(convertFToD(_c),zDiag);
      double data2NormSqInit = dot(dataResInit,dataResInit);
      double bPenalty2NormSqInit = dotAll(bPenaltyResInit,bPenaltyResInit)+eps;
      double cPenalty2NormSqInit = dotAll(cPenaltyResInit,cPenaltyResInit);
      double sumInit = data2NormSqInit+bPenalty2NormSqInit+cPenalty2NormSqInit;
      double partOfDifference = 0.001;
      double difference = _firstAll2NormSq-data2NormSqInit+1e-10;
      double alpha = ((partOfDifference*difference+data2NormSqInit)-data2NormSqInit)/bPenalty2NormSqInit;
      //double alpha = (_firstAll2NormSq-data2NormSqInit)/(2.0*bPenalty2NormSqInit);
      return alpha;
    }
    public float[] getDeltaBDeltaC() {
      return _deltaBDeltaC;
    }

    private int _nb,_nc,_kb,_kc;
    private double _scaleW;
    private float _stepLength;
    private double _data2NormSq,_bPenalty2NormSq,_cPenalty2NormSq,_firstAll2NormSq;
    private float[] _b,_bNew,_c,_cNew,_deltaB,_deltaC,_deltaBDeltaC;
    private float[] _bPenaltyRes,_bPenaltyResInit,_cPenaltyRes,_cPenaltyResInit;
    private float[][] _u,_f,_g,_dataRes,_dataResInit;
    private double[] _wDiag,_zDiag,_wTransWDiag,_zTransZDiag,_wTransWB,_zTransZC;
    private DMatrix _x,_v,_xNew,_vNew;
  }

  /*
  private class ObjectiveFunctionWithAlphaC implements BrentMinFinder.Function {
    public void setParameters(int nc, int kc, float[] c,
                              int nb, int kb, float[] b,
                              float[] u, float[] f, float[] g, 
                              DMatrix x, DMatrix v,
                              double scaleZ) 
    {
      _nb = nb;
      _nc = nc;
      _kc = kc;
      _kb = kb;
      _b = copy(b);
      _c = copy(c);
      _u = u;
      _f = f;
      _g = g;
      _x = x;
      _v = v;
      _scaleZ = scaleZ;
    }
    public double evaluate(double alphaC) {
      //compute penalizations
      _zDiag = computeZDiagonal(_nc,_kc,alphaC,_scaleZ);
      _wDiag = new double[_nb];
      _zTransZDiag = computeZTransposeZDiagonal(_zDiag);//Diagonal of Z'Z
      _wTransWDiag = new double[_nb];
      _zTransZC = computeZTransposeZC(_zTransZDiag,convertFToD(_c));//Z'Zc
      _wTransWB = new double[_nb];

      //add penalizations
      _xNew = addWWDiagAndZZDiagToX(_x,_wTransWDiag,_zTransZDiag);
      _vNew = subtractWWbAndZZcFromV(_v,_wTransWB,_zTransZC);

      //solve system
      DMatrixChd chd = new DMatrixChd(_xNew);
      _deltaBDeltaC = new float[_nb+_nc];
      try {
        _deltaBDeltaC = solveForSearchDirection(_xNew,_vNew);
      }      
      catch (IllegalStateException e) {
        return Float.MAX_VALUE;
      }
      _deltaB = copy(_nb,0,_deltaBDeltaC);
      _deltaC = copy(_nc,_nb,_deltaBDeltaC);

      //update deltab and deltac and c and b
      _stepLength = solveForStepLength(_nb,_kb,_b,_deltaB,_nc,_kc,_c,_deltaC,_u,_f,_g,_wDiag,_zDiag);
      _deltaB = mul(_stepLength,_deltaB);
      _deltaC = mul(_stepLength,_deltaC);
      _bNew = add(_b,_deltaB);
      _cNew = add(_c,_deltaC);
      
      //compute residuals
      _dataRes = computeDataResidual(_nc,_kc,_cNew,_nb,_kb,_bNew,_u,_f,_g);
      _bPenaltyRes = computeBPenaltyResidual(convertFToD(_bNew),_wDiag);
      _cPenaltyRes = computeCPenaltyResidual(convertFToD(_cNew),_zDiag);
      _data2NormSq = dot(_dataRes,_dataRes);
      _bPenalty2NormSq = dotAll(_cPenaltyRes,_cPenaltyRes);
      _cPenalty2NormSq = dotAll(_bPenaltyRes,_bPenaltyRes);
      trace("alphaC = "+alphaC);
      trace("dottotal = "+(_data2NormSq+_bPenalty2NormSq+_cPenalty2NormSq));
      trace("rmsrftotal = "+rmsOfObjectiveFunction(_dataRes,_bPenaltyRes,_cPenaltyRes));
      return _data2NormSq+_bPenalty2NormSq+_cPenalty2NormSq;
    }
    public float[] getDeltaBDeltaC() {
      return _deltaBDeltaC;
    }

    private int _nb,_nc,_kb,_kc;
    private double _scaleZ;
    private float _stepLength;
    private double _data2NormSq,_bPenalty2NormSq,_cPenalty2NormSq;
    private float[] _b,_bNew,_c,_cNew,_deltaB,_deltaC,_deltaBDeltaC;
    private float[] _dataRes,_bPenaltyRes,_cPenaltyRes;
    private float[] _u,_f,_g;
    private double[] _wDiag,_zDiag,_wTransWDiag,_zTransZDiag,_wTransWB,_zTransZC;
    private DMatrix _x,_v,_xNew,_vNew;
  }

  private class ObjectiveFunctionWithAlphaCAlphaB implements BrentMinFinder.Function {
    public void setParameters(int nc, int kc, float[] c,
                              int nb, int kb, float[] b,
                              float[] u, float[] f, float[] g, 
                              DMatrix x, DMatrix v,
                              double scaleWandZ) 
    {
      _nb = nb;
      _nc = nc;
      _kc = kc;
      _kb = kb;
      _b = copy(b);
      _c = copy(c);
      _u = u;
      _f = f;
      _g = g;
      _x = x;
      _v = v;
      _scaleWandZ = scaleWandZ;
    }
    public double evaluate(double alphaBAndC) {
      //compute penalizations
      _zDiag = computeZDiagonal(_nc,_kc,alphaBAndC,_scaleWandZ);
      _wDiag = computeWDiagonal(_nb,_kb,alphaBAndC,_scaleWandZ);
      _zTransZDiag = computeZTransposeZDiagonal(_zDiag);//Diagonal of Z'Z
      _wTransWDiag = computeWTransposeWDiagonal(_wDiag);//Diagonal of W'W
      _zTransZC = computeZTransposeZC(_zTransZDiag,convertFToD(_c));//Z'Zc
      _wTransWB = computeWTransposeWB(_wTransWDiag,convertFToD(_b));//W'Wb

      //add penalizations
      _xNew = addWWDiagAndZZDiagToX(_x,_wTransWDiag,_zTransZDiag);
      _vNew = subtractWWbAndZZcFromV(_v,_wTransWB,_zTransZC);

      //solve system
      DMatrixChd chd = new DMatrixChd(_xNew);
      _deltaBDeltaC = new float[_nb+_nc];
      try {
        _deltaBDeltaC = solveForSearchDirection(_xNew,_vNew);
      }      
      catch (IllegalStateException e) {
        return Float.MAX_VALUE;
      }
      _deltaB = copy(_nb,0,_deltaBDeltaC);
      _deltaC = copy(_nc,_nb,_deltaBDeltaC);

      //update deltab and deltac and c and b
      _stepLength = solveForStepLength(_nb,_kb,_b,_deltaB,_nc,_kc,_c,_deltaC,_u,_f,_g,_wDiag,_zDiag);
      _deltaB = mul(_stepLength,_deltaB);
      _deltaC = mul(_stepLength,_deltaC);
      _bNew = add(_b,_deltaB);
      _cNew = add(_c,_deltaC);
      
      //compute residuals
      _dataRes = computeDataResidual(_nc,_kc,_cNew,_nb,_kb,_bNew,_u,_f,_g);
      _bPenaltyRes = computeBPenaltyResidual(convertFToD(_bNew),_wDiag);
      _cPenaltyRes = computeCPenaltyResidual(convertFToD(_cNew),_zDiag);
      _data2NormSq = dot(_dataRes,_dataRes);
      _bPenalty2NormSq = dotAll(_cPenaltyRes,_cPenaltyRes);
      _cPenalty2NormSq = dotAll(_bPenaltyRes,_bPenaltyRes);
      trace("alphaBAndC = "+alphaBAndC);
      trace("dottotal = "+(_data2NormSq+_bPenalty2NormSq+_cPenalty2NormSq));
      trace("rmsrftotal = "+rmsOfObjectiveFunction(_dataRes,_bPenaltyRes,_cPenaltyRes));
      return _data2NormSq+_bPenalty2NormSq+_cPenalty2NormSq;
    }
    public float[] getDeltaBDeltaC() {
      return _deltaBDeltaC;
    }

    private int _nb,_nc,_kb,_kc;
    private double _scaleWandZ;
    private float _stepLength;
    private double _data2NormSq,_bPenalty2NormSq,_cPenalty2NormSq;
    private float[] _b,_bNew,_c,_cNew,_deltaB,_deltaC,_deltaBDeltaC;
    private float[] _dataRes,_bPenaltyRes,_cPenaltyRes;
    private float[] _u,_f,_g;
    private double[] _wDiag,_zDiag,_wTransWDiag,_zTransZDiag,_wTransWB,_zTransZC;
    private DMatrix _x,_v,_xNew,_vNew;
  }
  */
  /**
   * This class is used to calculate the objective function when no
   * penalization term exists. It's main purpose is to be used with 
   * the BrentMinFinder to compute the optimal step length 
   * along the search direction in the Gauss-Newton method.
   */
  private class ObjectiveFunction implements BrentMinFinder.Function {
    public void setParameters(int nc, int kc, float[] c, float[] deltac, 
                              int nb, int kb, float[] b, float[] deltab,
                              float[] u, float[] f, float[] g,
                              double[] wDiag, double[] zDiag)
    {
      _nc = nc;
      _nb = nb;
      _kc = kc;
      _kb = kb;
      _c = copy(c);
      _b = copy(b);
      _deltaC = copy(deltac);
      _deltaB = copy(deltab);
      _u = u;
      _f = f;
      _g = g;
      _wDiag = wDiag;
      _zDiag = zDiag;
    }
    public double evaluate(double stepLength) {
      float[] deltab = mul((float) stepLength,_deltaB);
      float[] deltac = mul((float) stepLength,_deltaC);
      float[] b = add(_b,deltab);
      float[] c = add(_c,deltac);

      float[] dataRes = computeDataResidual(_nc,_kc,c,_nb,_kb,b,_u,_f,_g);
      float[] bPenaltyRes = computeBPenaltyResidual(convertFToD(b),_wDiag);
      float[] cPenaltyRes = computeCPenaltyResidual(convertFToD(c),_zDiag);
      double data2NormSq = dot(dataRes,dataRes);
      double bPenalty2NormSq = dotAll(bPenaltyRes,bPenaltyRes);
      double cPenalty2NormSq = dotAll(cPenaltyRes,cPenaltyRes);
      return data2NormSq+cPenalty2NormSq+bPenalty2NormSq;
    }
    private float[] _c,_b,_deltaC,_deltaB,_u,_f,_g;
    private int _nc,_nb,_kc,_kb;
  }

  private class ObjectiveFunction2D implements BrentMinFinder.Function {
    public void setParameters(int nc, int kc, float[] c, float[] deltac, 
                              int nb, int kb, float[] b, float[] deltab,
                              float[][] u, float[][] f, float[][] g,
                              double[] wDiag, double[] zDiag)
    {
      _nc = nc;
      _nb = nb;
      _kc = kc;
      _kb = kb;
      _c = copy(c);
      _b = copy(b);
      _deltaC = copy(deltac);
      _deltaB = copy(deltab);
      _u = u;
      _f = f;
      _g = g;
      _wDiag = wDiag;
      _zDiag = zDiag;
    }
    public double evaluate(double stepLength) {
      float[] deltab = mul((float) stepLength,_deltaB);
      float[] deltac = mul((float) stepLength,_deltaC);
      float[] b = add(_b,deltab);
      float[] c = add(_c,deltac);

      float[][] dataRes = computeDataResidual(_nc,_kc,c,_nb,_kb,b,_u,_f,_g);
      float[] bPenaltyRes = computeBPenaltyResidual(convertFToD(b),_wDiag);
      float[] cPenaltyRes = computeCPenaltyResidual(convertFToD(c),_zDiag);
      double data2NormSq = dot(dataRes,dataRes);
      double bPenalty2NormSq = dotAll(bPenaltyRes,bPenaltyRes);
      double cPenalty2NormSq = dotAll(cPenaltyRes,cPenaltyRes);
      return data2NormSq+cPenalty2NormSq+bPenalty2NormSq;
    }
    private float[] _c,_b,_deltaC,_deltaB;
    private float[][] _u,_f,_g;
    private int _nc,_nb,_kc,_kb;
  }

  /**
   * Builds the approximated Hessian and gradient for the Gauss-Newton method.
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
    float[] q0 = warp.applyS(u,applyC(nb,kb,b,g));
    float[] qi = new float[nt];
    float[] qj = new float[nt];
    float[] pi = new float[nt];
    float[] pj = new float[nt];
    //P'P and Q'P and P'Q
    for (int ib=0,lagi=kb; ib<nb; ++ib,++lagi) {
      pi = applyC(nc,kc,c,warp.applyS(u,delay(lagi,g)));
      double vi = dot(pi,r);
      v.set(ib,0,-vi);
      for (int jb=0,lagj=kb; jb<=ib; ++jb,++lagj) {
        pj = applyC(nc,kc,c,warp.applyS(u,delay(lagj,g)));
        double xij = dot(pi,pj);
        x.set(ib,jb,xij);
        x.set(jb,ib,xij);
      }
      for (int jc=0,lagj=kc; jc<nc; ++jc,++lagj) {
        qi = delay(lagj,q0);
        vi = dot(qi,r);
        v.set(nb+jc,0,-vi);
        double xij = dot(qi,pi);
        x.set(nb+jc,ib,xij);
        x.set(ib,nb+jc,xij);
      }
    }

    //Q'Q
    for (int ic=nb,lagi=kc; ic<nbc; ++ic,++lagi) {
      qi = delay(lagi,q0);
      for (int jc=nb,lagj=kc; jc<=ic; ++jc,++lagj) {
        qj = delay(lagj,q0);
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
    int nx = u.length;
    int nbc = nb+nc;
    Warper warp = new Warper();
    DMatrix v = new DMatrix(nbc,1);
    DMatrix x = new DMatrix(nbc,nbc);
    float[][] q0 = warp.applyS(u,applyC(nb,kb,b,g));
    float[][] qi = new float[nx][nt];
    float[][] qj = new float[nx][nt];
    float[][] pi = new float[nx][nt];
    float[][] pj = new float[nx][nt];
    double vi = 0.0;
    double xij = 0.0;
    //P'P and Q'P and P'Q
    for (int ib=0,lagi=kb; ib<nb; ++ib,++lagi) {
      trace("ib = "+ib);
      pi = applyC(nc,kc,c,warp.applyS(u,delay(lagi,g)));
      vi = dot(pi,r);
      v.set(ib,0,-vi);
      for (int jb=0,lagj=kb; jb<=ib; ++jb,++lagj) {
        pj = applyC(nc,kc,c,warp.applyS(u,delay(lagj,g)));
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
  
  /**
   * Modifies the approximated Hessian and gradient when the outer lags in b 
   * are penalized.
   */
  private DMatrix[] penalize(double[] wDiag, double[] zDiag, double[] b, double[] c,
    DMatrix x, DMatrix v) 
  {
    double[] wTransWDiag = computeWTransposeWDiagonal(wDiag);//Diagonal of W'W
    double[] zTransZDiag = computeZTransposeZDiagonal(zDiag);//Diagonal of Z'Z
    double[] wTransWB = computeWTransposeWB(wTransWDiag,b);//W'Wb
    double[] zTransZC = computeZTransposeZC(zTransZDiag,c);//Z'Zc
    x = addWWDiagAndZZDiagToX(x,wTransWDiag,zTransZDiag);
    v = subtractWWbAndZZcFromV(v,wTransWB,zTransZC);
    return new DMatrix[]{x,v};
  }

  /*
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
      float[] pi = applyC(nc,kc,c,sdig);
      double vi = dot(pi,r);
      double xij = dot(pi,pi);
      v.set(ib,0,-vi);
      x.set(ib,ib,xij);
    }
    Warper warp = new Warper();
    float[] ag = applyC(nb,kb,b,g);
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
      float[][] pi = applyC(nc,kc,c,sdig);
      double vi = dot(pi,r);
      double xij = dot(pi,pi);
      v.set(ib,0,-vi);
      x.set(ib,ib,xij);
    }
    Warper warp = new Warper();
    float[][] ag = applyC(nb,kb,b,g);
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
  */

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
  private float twoNorm(DMatrix x) {
    int nColumns = x.getN();
    Check.argument(nColumns>=1,"Number of columns in x >= 1");
    return (float)sqrt(sum(pow(x.getArray(),2.0f)));
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
  private void printX(DMatrix x,String type) {
      trace(type);
      trace(x.toString());
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






