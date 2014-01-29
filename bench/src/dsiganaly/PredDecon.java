/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dsiganaly;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.dsp.Conv.*;
import edu.mines.jtk.mosaic.SimplePlot;
import linalgebra.ToeplitzRecursion;
 


/**
 * Computes the prediction operator coefficients and therefore the 
 * prediction error operator coefficients for a lag distance of 1, which
 * means the prediction error coefficents are the inverse of the wavelet
 * in the seismic trace.
 * @author Chris Graziano, CWP
 * @version 2013.12.18
 */
public class PredDecon {
  
  /**
   * Calculates the prediction error operator
   * coefficients from the input signal and the number of coefficients
   * that wish to be calculated. 
   * Note: The first prediction error coefficient is 1.
   * This particular method only operates
   * on a single signal. 
   * @param x the input signal 
   * @param kx the sample index of x[0]
   * @param na number of prediction error coefficients to be calculated
   */
  public PredDecon(float[] x, int npec) {
    int nx = x.length;
    int kr = 0;
    //According to Taner and Robinson, only autocorrelation coefficients r0
    //through rn-1 are in the Toeplitz matrix, while the right hand side
    //vector includes coefficients r1 through rn. This yields n-1 prediction
    //operator coefficents, but these coefficients yield n prediction error
    //operator coefficients, which is why the number of autocorrelation
    //coefficients) is na.
    int nr = npec;
    float[] r = new float[nr];//0 to n lag autocorrelation coefficients
    nx = x.length;
    int kx = 0;
    xcor(nx,kx,x,nx,kx,x,nr,0,r);

    //the first row of the toeplitz matrix needed to calculate the
    //prediction operator coefficients, will contain r0 to rn-1
    int nrm1 = nr-1;
    float[] toepRow1 = new float[nrm1];

    //The right hand side of the linear equations Ra=b to compute
    //prediction operator coefficients, will contain r1 to rn 
    //autocorrelation coefficients
    float[] g = new float[nrm1];

    //Fill toepRow1 and b with the corresponding autocorrelation 
    //coefficients
    for (int i=0; i<nrm1; ++i) {
      toepRow1[i] = r[i  ];
      g[i]        = r[i+1];
    }

    //Solve for the prediction coefficients that minimize the 
    //squared error between the predicted x value and the actual x 
    //value for each sample of the x array.
    _pc = ToeplitzRecursion.solve(toepRow1,g);
    int npc = _pc.length;

    //create prediction error coefficients from the prediction
    //coefficients
    _pec = new float[npec];
    _pec[0] = 1.0f;
    for (int i=1; i<npec; ++i) 
      _pec[i] = -_pc[i-1];
  }

  /**
   * Calculates the prediction error operator
   * coefficients from the input gather and the number of coefficients
   * that wish to be calculated. 
   * Note: The first prediction error coefficient is 1.
   * This particular method only operates
   * on a single gather. 
   * @param x the input signal 
   * @param kx the sample index of x[0]
   * @param na number of prediction error coefficients to be calculated
   */
  public PredDecon(float[][] x, int na, int ka) {
    int npec = na;
    int nx = x.length;
    ka = 0;
    //According to Taner and Robinson, only autocorrelation coefficients r0
    //through rn-1 are in the Toeplitz matrix, while the right hand side
    //vector includes coefficients r1 through rn. This yields n-1 prediction
    //operator coefficents, but these coefficients yield n prediction error
    //operator coefficients, which is why the number of autocorrelation
    //coefficients) is na.
    int nr = npec;
    float[] r = getAutoCorrGather(x, nr);
//r = div(r,max(r));

    //the first row of the toeplitz matrix needed to calculate the
    //prediction operator coefficients, will contain r0 to rn-1
    int nrm1 = nr-1;
    float[] toepRow1 = new float[nrm1];

    //The right hand side of the linear equations Ra=b to compute
    //prediction operator coefficients, will contain r1 to rn 
    //autocorrelation coefficients
    float[] g = new float[nrm1];

    //Fill toepRow1 and b with the corresponding autocorrelation 
    //coefficients
    for (int i=0; i<nrm1; ++i) {
      toepRow1[i] = r[abs(i+ka)];
      g[i]        = r[abs(i+1+ka)];
    }

    //Solve for the prediction coefficients that minimize the 
    //squared error between the predicted x value and the actual x 
    //value for each sample of the x array.
    _pc = ToeplitzRecursion.solve(toepRow1,g);
    int npc = _pc.length;
    dump(_pc);

    //create prediction error coefficients from the prediction
    //coefficients
    _pec = new float[npec];
    for (int i=0; i<npec; ++i) {
      if (i==-ka) _pec[i] = 1.0f;
      else if (i < -ka) _pec[i] = -_pc[i];
      else if (i > -ka) _pec[i] = -_pc[i-1];
    }
  }


  public float[] getPredErrorCoef() {
    return _pec;
  }

  public float[] getPredCoef() {
    return _pc;
  }

  /**
   * This method defines the autocorrelation of the
   * gather as the sum of the autocorrelations of each
   * trace in the gather
   * @param f gather
   */
  private float[] getAutoCorrGather(float[][] f, int lag) {
    int nx = f.length;
    int nt = f[0].length;

    float[] r = new float[lag];
    float[] rsum = new float[lag];
    int kr = 0;
    int kt = 0;
    for (int x=0; x<nx; ++x) {
      xcor(nt,kt,f[x],nt,kt,f[x],lag,kr,r);
      add(rsum,r,rsum);
    }
    return rsum;
  }

  /**
   * This method defines the autocorrelation of the
   * gather as the autocorrelation of the sum of the gather 
   * int the 2nd dimension(offset)
   * @param f gather
   * @param lag number of autocorrelation coefficients
   */
  private float[] getAutoCorrGather2(float[][] f,int lag) {
    int nx = f.length;
    int nt = f[0].length;
    int kt = 0;

    float[] r = new float[lag];
    float[] fsum = new float[nt];
    int kr = 0;
    for (int x=0; x<nx; ++x) {
      add(fsum,f[x],fsum);
    }
    xcor(nt,kt,fsum,nt,kt,fsum,lag,kr,r);
    return r;
  }

  private float[] _pec, _pc;
}


