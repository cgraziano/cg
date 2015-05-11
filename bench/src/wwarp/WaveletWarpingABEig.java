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
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates two wavelets from alignment by warping sequences or images.
 * The two sequences or images are assumed to have been convolved with 
 * different wavelets. Warping of one sequence or image to align with 
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

 public class WaveletWarpingABEig {
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


  /**
   * Calculates two inverse wavelets (a and b) by warping one sequence to another.
   * The sequences are related by warping such that f[t] ~ g[u[t]] if both sequences had 
   * the same wavelet. 
   * This method solves for the eigenvector corresponding to the smallest eigenvalue
   * of matrix Y. 
   * Y = T'T
   * T = [F | -SG]  (partitioned matrix)
   * This eigenvector can be thought as a partitioned vector that contains 
   * a and b: [a]
   *          |-|
   *          [b]
   * The returned array is [a]
   *                       |-|
   *                       [b]
   */
  public float[] getInverseAB(
    int na, int ka, int nb, int kb, 
    float[] u, float[] f, float[] g)
  {
    Check.argument(-na<ka,"-na<ka");
    Check.argument(-nb<kb,"-nb<kb");
    Check.argument(ka<=0,"ka<=0");
    Check.argument(kb<=0,"kb<=0");

    //Set default weights
    if (_weights==false) {
      int nf = f.length;
      _w1D = fillfloat(1.0f,nf);
    }
    
    int nab = na+nb;
    DMatrix y = computeY(na,ka,nb,kb,u,f,g);
    _y = y;

    float[] ab = new float[nab];
    DMatrixEvd evd = new DMatrixEvd(y);
    //System.out.println("V");
    //System.out.println(evd.getV().toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println("D");
    //System.out.println(evd.getD().toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println(evd.getD().get(0,0));
    _eigvect = evd.getV().get();
    _eigval = evd.getD().get();

    for (int i=0; i<nab; ++i)
      ab[i] = (float) evd.getV().get(i,0);
    _ab = ab;

    return ab;
  }

  
  /**
   * Calculates two inverse wavelets (a and b) by warping one sequence to another.
   * The sequences are related by warping such that f[t] ~ g[u[t]] if both sequences had 
   * the same wavelet. 
   * This method solves for the eigenvector corresponding to the smallest eigenvalue
   * of matrix Y. 
   * Y = T'T
   * T = [F | -SG]  (partitioned matrix)
   * This eigenvector can be thought as a partitioned vector that contains 
   * a and b: [a]
   *          |-|
   *          [b]
   */
  public float[] getInverseAB(
    int na, int ka, int nb, int kb, 
    float[][] u, float[][] f, float[][] g)
  {
    Check.argument(-na<ka,"-na<ka");
    Check.argument(-nb<kb,"-nb<kb");
    Check.argument(ka<=0,"ka<=0");
    Check.argument(kb<=0,"kb<=0");
    
    int nab = na+nb;
    DMatrix y = computeY(na,ka,nb,kb,u,f,g);
    _y = y;
    
    

    float[] ab = new float[na];
    DMatrixEvd evd = new DMatrixEvd(y);
    //System.out.println(evd.getV().toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println("D");
    //System.out.println(evd.getD().toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println(evd.getD().get(0,0));
    _eigvect = evd.getV().get();
    _eigval = evd.getD().get();

    for (int i=0; i<nab; ++i)
      ab[i] = (float) evd.getV().get(i,0);
    _ab = ab;

    return ab;
  }

  /**
   * Gets the inverse wavelet a from ab.
   */
  public float[] getA(int na, float[] ab) {
    float[] a = new float[na];
    for (int iab=0; iab<na; ++iab)
      a[iab] = ab[iab];
    return a;
  }
  /**
   * Gets the inverse wavelet b from ab.
   */
  public float[] getB(int nb, float[] ab) {
    int nab = ab.length;
    int b0 = nab-nb;
    float[] b = new float[nb];
    for (int iab=b0; iab<nab; ++iab) 
      b[iab-b0] = ab[iab];
    return b;
  }


  public DMatrix computeY(int na, int ka, int nb, int kb, 
    float[] u, float[] f, float[] g) {
    int nab = na+nb;
    DMatrix y = new DMatrix(nab,nab);
    Warper warp = new Warper();
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      //F'F
      float[] dif = delay(lagi,f);
      for (int ja=0,lagj=ka; ja<na; ++ja,++lagj) {
        float[] djf = delay(lagj,f);
        double ffij = dot(dif,djf);
        y.set(ia,ja,ffij);
      }

      //-F'SG and -S'G'F
      for (int ja=na,lagj=kb; ja<nab; ++ja,++lagj) {
        float[] djg = delay(lagj,g);
        float[] sldjg = warp.applyS(u,djg);
        double fsgij = dot(dif,sldjg);
        y.set(ia,ja,-fsgij);
        y.set(ja,ia,-fsgij);
      }
    }

    //G'S'SG
    for (int ia=na,lagi=kb; ia<nab; ++ia,++lagi) {
      float[] dig = delay(lagi,g);
      for (int ja=na,lagj=kb; ja<nab; ++ja,++lagj) {
        float[] djg = delay(lagj,g);
        float[] sldig = warp.applyS(u,dig);
        float[] sldjg = warp.applyS(u,djg);
        double gsgsij = dot(sldig,sldjg);
        y.set(ia,ja,gsgsij);
      }
    }

    /*
    //F'F
    for (int ia=0,lagi=kaf; ia<naf; ++ia,++lagi) {
      float[] dfi = delay(lagi,f);
      for (int ja=0,lagj=kaf; ja<naf; ++ja,++lagj) {
        float[] dfj = delay(lagj,f);
        double ffij = dot(dfi,dfj);
        y.set(ia,ja,ffij);
      }
    }

    //-F'SG and -S'G'F
    for (int ia=0,lagi=kaf; ia<naf; ++ia,++lagi) {
      float[] dfi = delay(lagi,f);
      for (int ja=naf,lagj=kag; ja<na; ++ja,++lagj) {
        float[] dgj = delay(lagj,g);
        float[] sldgj = applyDLSU(u,dgj);
        double fsgij = dot(dfi,sldgj);
        y.set(ia,ja,-fsgij);
        y.set(ja,ia,-fsgij);
      }
    }

    //G'S'SG
    for (int ia=naf,lagi=kag; ia<na; ++ia,++lagi) {
      float[] dgi = delay(lagi,g);
      for (int ja=naf,lagj=kag; ja<na; ++ja,++lagj) {
        float[] dgj = delay(lagj,g);
        float[] sldgi = applyDLSU(u,dgi);
        float[] sldgj = applyDLSU(u,dgj);
        double gsgsij = dot(sldgi,sldgj);
        y.set(ia,ja,gsgsij);
      }
    }
    */

    //Stabilize Y to be symmetric positive definite
    for (int ia=0; ia<nab; ++ia)
      y.set(ia,ia,y.get(ia,ia)*_sfac);
    return y;
  }

  public DMatrix computeY(int na, int ka, int nb, int kb, 
    float[][] u, float[][] f, float[][] g) {
    int nab = na+nb;
    DMatrix y = new DMatrix(nab,nab);
    Warper warp = new Warper();
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      //F'F
      float[][] dif = delay(lagi,f);
      for (int ja=0,lagj=ka; ja<na; ++ja,++lagj) {
        float[][] djf = delay(lagj,f);
        double ffij = dot(dif,djf);
        y.set(ia,ja,ffij);
      }
      //-F'SG and -S'G'F
      for (int ja=na,lagj=kb; ja<nab; ++ja,++lagj) {
        float[][] djg = delay(lagj,g);
        float[][] sldjg = warp.applyS(u,djg);
        double fsgij = dot(dif,sldjg);
        y.set(ia,ja,-fsgij);
        y.set(ja,ia,-fsgij);
      }
    }

    //G'S'SG
    for (int ia=na,lagi=kb; ia<nab; ++ia,++lagi) {
      float[][] dig = delay(lagi,g);
      for (int ja=na,lagj=kb; ja<nab; ++ja,++lagj) {
        float[][] djg = delay(lagj,g);
        float[][] sldig = warp.applyS(u,dig);
        float[][] sldjg = warp.applyS(u,djg);
        double gsgsij = dot(sldig,sldjg);
        y.set(ia,ja,gsgsij);
      }
    }

    //Stabilize Y to be symmetric positive definite
    for (int ia=0; ia<nab; ++ia)
      y.set(ia,ia,y.get(ia,ia)*_sfac);
    return y;
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
   * Returns the rms value of the image/trace.
   */
  public float rms(float[] x) {
    return (float)sqrt(dot(x,x)/x.length);
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

  public DMatrix getY() {
    return _y;
  }

  public String printY() {
    return _y.toString();
  }

  public float[] getEigVector(int i) {
    int nx = _eigvect[i].length;
    float[] eigvect = new float[nx];

    for(int ix=0; ix<nx; ++ix) {
      eigvect[ix] = (float) _eigvect[ix][i];
    }
    return eigvect;
  }

  public float getEigVal(int i) {
    return (float) _eigval[i][i];
  }

  /**
   * An indicator of if the null space of Ya=0 consists of one vector is 
   * if the smallest eigenvalue of matrix Y is significantly smaller than the 
   * second smallest eigenvalue of matrix Y. If this ratio is less than one, 
   * the null space of Ya=0 consists of one vecotr. If this ratio is near, equal, or
   * above one, the null space of Ya=0 will consist of more than one vector, which that
   * there exists an infinitie number of possible answers in the null space of Ya=0.
   */
  public float getEig01Ratio() {
    return (float) (_eigval[0][0]/_eigval[1][1]);
  }

  public float[] getEigVals() {
    int na = _eigval.length;
    float[] eigvals = new float[na];
    for (int i=0; i<na; ++i)
      eigvals[i] = (float) _eigval[i][i];
    return eigvals;
  }

  /**
   * Returns the sum of the squared differences between Fa and SGa
   * for specific n, nag, kaf, and kag values (inverse wavelet sampling parameters).
   * 
   * 
   * @param naf the number of inverse wavelet coefficients in af.
   * @param nag the number of inverse wavelet coefficients in ag.
   * @param kaf the sample index of af corresponding to 0 seconds.
   * @param kag the sample index of ag corresponding to 0 seconds.
   */
  public float getSumSqDiff(int nab) {
    double[] abd = convertFToD(_ab);
    DMatrix a = new DMatrix(nab,1,abd);
    DMatrix at = a.transpose();
    return (float) (at.times(_y)).times(a).get(0,0);
  }

  /**
   * Returns the sum of the squared differences between Fa and SGa
   * if a is an impulse centered at -ka. 
   * @param na the number of inverse wavelet coefficients in a.
   * @param ka the sample index of a[0].
   */
  public float getSumSqDiffNoWaveletEst(int na, int ka, int kb) {
    return (float) (_y.get(-ka,-ka)-_y.get(-kb+na,-ka)-_y.get(-ka,-kb+na)+_y.get(-kb,-kb));
  }

  /**
   *  A value above 1.0 indicates a well-determined wavelet.
   */
  public double getWDMeasure() {
    int na = _eigval.length;
    double lambda0 = _eigval[0][0];
    double lambda1 = _eigval[1][1];
    double lambdan = _eigval[na-1][na-1];
    double l0n = lambda0/lambdan;
    double l1n = lambda1/lambdan;
    double eps1 = Math.ulp(1.0);
    System.out.println("lambda0/lambdan = "+l0n);
    System.out.println("lambda1/lambdan = "+l1n);
    return (l1n-l0n)/eps1;
  }

  /**
   * Returns y(t) = x(t-lag).
   * Should only be used to build delayed trace or image for plotting purposes.
   */
  public static float[] tempdelay(int lag, float[] x) {
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

  public void testSPD(int nrow, int ncol) {
    double[][] r = randdouble(ncol,nrow);
    r = mul(r,2.0);
    r = sub(r,1.0);
    DMatrix c = new DMatrix(r);
    DMatrix y = (c.transpose()).times(c);
    DMatrixEvd evd = new DMatrixEvd(y);
    double[][] eigval = evd.getD().get();
    System.out.println("////////testSPD/////////");
    System.out.println("If the following eigevalue is negative, the random matrix is not pd");
    System.out.println("1st eigenvalue = "+eigval[0][0]);
    System.out.println("////////////////////////");
  }



///////////////////////////////////////////////////////////////////////////
  // private

  private double _wha = 0.1;
  private double _sfac = 1.0;
  private int _itmin = -1;
  private int _itmax = -1;
  private int _ng0 = 0;//staring size of array/image.
  private float[] _w1D;
  private float[] _ab;
  private float[][] _w2D;
  private double[][] _eigvect;
  private double[][] _eigval;
  private boolean _weights;
  private DMatrix _y;


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
      sum += (double)(x[it])*(double)(y[it])*(double)(w[it]);
    return sum;
  }

  private double dot(float[] x, float[] y, float[] w) {
    int nt = x.length;
    int itlo = (_itmin>0)?_itmin:0;
    int ithi = (_itmax>0)?_itmax:nt-1;
    double sum = 0.0;
    for (int it=itlo; it<=ithi; ++it) 
      sum += (double)(x[it])*(double)(y[it])*(double)(w[it]);
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

  private static void trace(String s) {
    System.out.println(s);
  }





 }
