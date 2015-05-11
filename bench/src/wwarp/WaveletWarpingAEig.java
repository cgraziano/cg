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

 public class WaveletWarpingAEig {
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
   * Calculates one inverse wavelet (a) by warping one sequence to another.
   * The sequences are related by warping such that f[t] ~ g[u[t]] if both sequences had 
   * the same wavelet. 
   * This method solves for the eigenvector corresponding to the smallest eigenvalue
   * of matrix Z. 
   * Z = D'D
   * D = (F-SG) 
   * The eigenvector is the wavelet in f and g.
   */
  public float[] getInverseA(
    int na, int ka, float[] u, float[] f, float[] g)
  {
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    //Set default weights (equal weight to all samples)
    if (_weights==false) {
      int nf = f.length;
      _w1D = fillfloat(1.0f,nf);
    }
    
    DMatrix z = computeZ(na,ka,u,f,g);
    _z = z;
    
    float[] a = new float[na];
    DMatrixEvd evd = new DMatrixEvd(z);
    //System.out.println("z");
    //System.out.println(z.toString());
    //System.out.println(evd.getD().toString());
    //System.out.println(evd.getV().toString());
    double[][][] eigvalvect = sort(evd.getD().get(), evd.getV().get());
    _eigval = eigvalvect[0];
    _eigvect = eigvalvect[1];
    //dump(_eigval);
    //dump(_eigvect);

    for (int i=0; i<na; ++i) {
      a[i] = (float) _eigvect[i][0];
    }
    _a = a;
    //System.out.println("z");
    //System.out.println(z.toString());
    //System.out.println("V");
    //System.out.println(evd.getV().toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println(evd.getV().get(0,1));//Smallest eigenvalue position (0,0) in D
    //System.out.println("D");
    //System.out.println(evd.getD().toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println(evd.getD().get(0,0));

    return a;
  }


  /**
   * Calculates one inverse wavelets a by warping one sequence to another.
   * The sequences are related by warping such that f[t] ~ g[u[t]] if both sequences had 
   * the same wavelet. 
   * This method solves for the eigenvector corresponding to the smallest eigenvalue
   * of matrix Z. 
   * Z = D'D
   * D = (F-SG) 
   * The eigenvector is the wavelet in f and g.
   */
  public float[] getInverseA(
    int na, int ka, float[][] u, float[][] f, float[][] g)
  {
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");

    //Set default weights (equal weight to all samples)
    if (_weights==false) {
      int nf = f.length;
      _w1D = fillfloat(1.0f,nf);
    }

    Stopwatch sw = new Stopwatch(); 
    sw.start();
    DMatrix z = computeZ(na,ka,u,f,g);
    sw.stop();
    System.out.println("time = "+sw.time());
    _z = z;
    
    float[] a = new float[na];
    DMatrixEvd evd = new DMatrixEvd(z);
    double[][][] eigvalvect = sort(evd.getD().get(), evd.getV().get());
    _eigval = eigvalvect[0];
    _eigvect = eigvalvect[1];
    //System.out.println("z");
    //System.out.println(z.toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println("V");
    //System.out.println(evd.getV().toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println(evd.getV().get(0,1));//Smallest eigenvalue position (0,0) in D
    //System.out.println("D");
    //System.out.println(evd.getD().toString());//Smallest eigenvalue position (0,0) in D
    //System.out.println(evd.getD().get(0,0));
    //_eigval = evd.getD().get();//eigvalvect[0];
    //_eigvect = evd.getV().get();//eigvalvect[1];
    //dump(_eigval);

    for (int i=0; i<na; ++i) {
      a[i] = (float) _eigvect[i][0];
    }
    _a = a;

    return a;
  }

  public DMatrix computeZ(int na, int ka, float[] u, float[] f, float[] g) {
    DMatrix z = new DMatrix(na,na);
    //Z = (F-SG)'(F-SG) = F'F+G'S'SG+F'SG+G'S'F
    double ffij = 0.0;
    double gssgij = 0.0;
    double fsgij = 0.0;
    double zij = 0.0;
    double zji = 0.0;
    Warper warp = new Warper();
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      System.out.println("ia = "+ia);
      float[] dfi = delay(lagi,f);
      float[] dgi = delay(lagi,g);
      float[] sldgi = warp.applyS(u,dgi);
      for (int ja=0,lagj=ka; ja<na; ++ja,++lagj) {
        float[] dfj = delay(lagj,f);
        float[] dgj = delay(lagj,g);
        float[] sldgj = warp.applyS(u,dgj);
        ffij = dot(dfi,dfj);
        gssgij = dot(sldgi,sldgj);
        fsgij = -dot(dfi,sldgj);
        zij = ffij + gssgij + fsgij + z.get(ia,ja);
        z.set(ia,ja,zij);
        zji = fsgij + z.get(ja,ia);
        z.set(ja,ia,zji);
      }
    }

    //Stabilize Z to be symmetric positive definite
    for (int ia=0; ia<na; ++ia)
      z.set(ia,ia,z.get(ia,ia)*_sfac);
    return z;
  }

  public DMatrix computeZ(final int na, final int ka, 
    final float[][] u, final float[][] f, final float[][] g) {
    //Z = (F-SG)'(F-SG) = F'F+G'S'SG+F'SG+G'S'F
    Parallel.setParallel(true);
    final DMatrix z = new DMatrix(na,na);
    //final double[][] z = new double[na][na];
    //F'F+G'S'SG+F'SG+G'S'F
    final Warper warp = new Warper();
    final float r = (float)((int)(warp.squeezing(1.0f,u)+1.0f));//Round up to the nearest integer.
    BandPassFilter bpf = constructBPF(r);
    int nu2 = u.length;
    final int nu1 = u[0].length;
    int ng2 = g.length;
    int ng1 = g[0].length;
    final float[][] uu = warp.upSampleLinear(nu1,1.0f,0.0f,u,r);
    Stopwatch sw = new Stopwatch();
    sw.start( );
    /*Parallel.loop(na,new Parallel.LoopInt() {
      public void compute(int ia) {
        int lagi = ka+ia;
    //for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
        System.out.println("ia = "+ia);
        float[][] dfi = delay(lagi,f);
        float[][] dgi = delay(lagi,g);
        float[][] sdgi = warp.applyS(r,nu1,uu,dgi);
        for (int ja=0,lagj=ka; ja<(ia+1); ++ja,++lagj) {
          float[][] dfj = delay(lagj,f);
          float[][] dgj = delay(lagj,g);
          float[][] sdgj = warp.applyS(r,nu1,uu,dgj);
          double ffij = dot(dfi,dfj);
          double gssgij = dot(sdgi,sdgj);
          double fsgij = -dot(dfi,sdgj);
          double gsfij = -dot(sdgi,dfj);
          double zij = ffij + gssgij + gsfij + fsgij;//z.get(ia,ja);
          //z[ia][ja] = zij;//+z[ia][ja];
          z.set(ia,ja,zij);
        }
      }
    });
    //}
    //Z is symmetric
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      for (int ja=0,lagj=ka; ja<ia+1; ++ja,++lagj) {
        //z[ja][ia] = z[ia][ja];
        z.set(ja,ia,z.get(ia,ja));
      }
    }
    sw.stop();
    //System.out.println(z.toString());
    System.out.println("time to construct z is "+sw.time());
    System.out.println("Z is symmetric? "+z.isSymmetric());
    */
    //F'F
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      System.out.println("ia = "+ia);
      float[][] dfi = delay(lagi,f);
      for (int ja=0,lagj=ka; ja<na; ++ja,++lagj) {
        float[][] dfj = delay(lagj,f);
        double zij = dot(dfi,dfj) + z.get(ia,ja);
        z.set(ia,ja,zij);
      }
    }
    //G'S'SG
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      System.out.println("ia = "+ia);
      float[][] dgi = delay(lagi,g);
      float[][] sldgi = warp.applyS(u,dgi);
      for (int ja=0,lagj=ka; ja<na; ++ja,++lagj) {
        float[][] dgj = delay(lagj,g);
        float[][] sldgj = warp.applyS(u,dgj);
        double zij = dot(sldgi,sldgj) + z.get(ia,ja);
        z.set(ia,ja,zij);
      }
    }
    //-F'SG+-G'S'F
    for (int ia=0,lagi=ka; ia<na; ++ia,++lagi) {
      System.out.println("ia = "+ia);
      float[][] dfi = delay(lagi,f);
      float[][] dgi = delay(lagi,g);
      float[][] sldgi = warp.applyS(u,dgi);
      for (int ja=0,lagj=ka; ja<na; ++ja,++lagj) {
        float[][] dfj = delay(lagj,f);
        float[][] dgj = delay(lagj,g);
        float[][] sldgj = warp.applyS(u,dgj);
        double fsldgij = -dot(dfi,sldgj);
        double zij = fsldgij + z.get(ia,ja);
        z.set(ia,ja,zij);
        double zji = fsldgij + z.get(ja,ia);
        z.set(ja,ia,zji);
      }
    }

    //Stabilize Z to be symmetric positive definite
    for (int ia=0; ia<na; ++ia)
      //z[ia][ia] = _sfac*z[ia][ia];
      z.set(ia,ia,z.get(ia,ia)*_sfac);
    //DMatrix mz = new DMatrix(z);
    return z;
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

  public float[][] getZ() {
    return convertDToF(_z.get());
  }

  public void getZ(float[][] z) {
    double[][] a = _z.get();
    int n2a = a.length;
    for (int i=0; i<n2a; ++i) 
      for (int j=0; j<n2a; ++j) 
        z[i][j] = (float) a[i][j];
  }

  public String printZ() {
    return _z.toString();
  }

  public double[] getEigVector(int i) {
    int nx = _eigvect[i].length;

    double[] eigvect = new double[nx];
    for(int ix=0; ix<nx; ++ix) {
      eigvect[ix] =  _eigvect[ix][i];
    }
    return eigvect;
  }

  public double getEigVal(int i) {
    return _eigval[i][i];
  }

  public static double ulp1() {
    System.out.println(Math.ulp(1.0));
    return Math.ulp(1.0);
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

  public float[] getEigVals() {
    int na = _eigval.length;
    float[] eigvals = new float[na];
    for (int i=0; i<na; ++i)
      eigvals[i] = (float) _eigval[i][i];
    return eigvals;
  }

  /**
   * Returns the sum of the squared differences between Fa and SGa
   * for specific na and kag values (inverse wavelet sampling parameters).
   * @param na the number of inverse wavelet coefficients in a.
   */
  public float getSumSqDiff(int na) {
    double[] ad = convertFToD(_a);
    DMatrix a = new DMatrix(na,1,ad);
    DMatrix at = a.transpose();
    return (float) (at.times(_z)).times(a).get(0,0) ;
  }

  /**
   * Returns the sum of the squared differences between Fa and SGa
   * if a is an impulse centered at -ka. 
   * @param na the number of inverse wavelet coefficients in a.
   * @param ka the sample index of a[0].
   */
  public float getSumSqDiffNoWaveletEst(int ka) {
    return (float) (_z.get(-ka,-ka));
  }

  public static void goTestUniqMeas() {
    System.out.println("Math.ulp(1.0) = "+Math.ulp(1.0));
    System.out.println("Math.ulp(1.0)/2.0 = "+Math.ulp(1.0)/2.0);
    System.out.println("1.0+Math.ulp(1.0) = "+(1.0+Math.ulp(1.0)));
    System.out.println("1.0+Math.ulp(1.0)/2.0 = "+(1.0+Math.ulp(1.0)/2.0));
    System.out.println("Math.ulp(100.0) = "+Math.ulp(100.0));
    System.out.println("Math.ulp(100.0)*0.9 = "+Math.ulp(100.0)*0.9);
    System.out.println("100.0+Math.ulp(100.0) = "+(100.0+Math.ulp(100.0)));
    System.out.println("100.0+Math.ulp(100.0)*0.9 = "+(100.0+Math.ulp(100.0)*0.9));
    System.out.println("Math.ulp(100.0) = "+Math.ulp(100.0));
    System.out.println("100.0*Math.ulp(1.0) = "+(100.0*Math.ulp(1.0)));
    System.out.println("100.0+Math.ulp(100.0) = "+(100.0+Math.ulp(100.0)));
    System.out.println("100.0+100.0*Math.ulp(1.0) = "+(100.0+100.0*Math.ulp(1.0)));
    System.out.println("Math.ulp(100.0) = "+Math.ulp(100.0));
    System.out.println("100.0*Math.ulp(1.0) = "+(100.0*Math.ulp(1.0)));
    System.out.println("100.0+Math.ulp(100.0) = "+(100.0+Math.ulp(100.0)));
    System.out.println("100.0+100.0*Math.ulp(1.0) = "+(100.0+100.0*Math.ulp(1.0)));
  }

///////////////////////////////////////////////////////////////////////////
  // private
  private double _sfac = 1.0;
  private int _itmin = -1;
  private int _itmax = -1;
  private int _ng0 = 0;//staring size of array/image.
  private float[] _a;
  private float[] _w1D;
  private float[][] _w2D;
  private double[][] _eigvect;
  private double[][] _eigval;
  private boolean _weights;
  private DMatrix _z;


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

  /**
   * Rearranges the 2D eigenvalue array and the corresponding 2D eigenvector array,
   * so that the diagonal of the eigenvalue array is in ascending order.
   */
  private double[][][] sort(double[][] eigval, double[][] eigvect) {
    int n = eigval.length;
    double[][][] eigvalvect;
    double[] diag = new double[n];
    diag = getDiagonal(eigval);
    if (!isSorted(diag)) {
      //System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!");
      //System.out.println("Eigenvalues are not sorted");
      double[] diagsort = new double[n];
      double[][] eigvalsort = new double[n][n];
      double[][] eigvectsort = new double[n][n];
      int[] indices = rampint(0,1,n);
      quickIndexSort(diag,indices);
      for (int i=0; i<n; ++i) {
        diagsort[i] = diag[indices[i]];
      }
      eigvalsort = DMatrix.diagonal(diagsort).get();

      for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
          eigvectsort[i][j] = eigvect[i][indices[j]];
        }
      }
      eigvalvect = new double[][][]{eigvalsort,eigvectsort};
      //System.out.println("Eigenvalues are now sorted");
    }
    else {
      //System.out.println("Eigenvalues are sorted");
      eigvalvect = new double[][][]{eigval,eigvect};
    }
    return eigvalvect;
  }

  private double[] getDiagonal(double[][] eigval) {
    int n = eigval.length;
    double[] diag = new double[n];
    for (int i=0; i<n; ++i) {
      diag[i] = eigval[i][i];
    }
    return diag;
  }

  private boolean isSorted(double[] x) {
    int nx = x.length;
    for (int i=1; i<nx; ++i) {
      if (x[i]<x[i-1])
        return false;
    }
    return true;
  }
  





  private static void trace(String s) {
    System.out.println(s);
  }
 }


