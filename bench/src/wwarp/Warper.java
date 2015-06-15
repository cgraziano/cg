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
 * Warps a uniformly-sampled signal, g(t), from time t to time u. 
 * The relation between time u and time t is defined by the uniformly-sampled function u(t). 
 * We are assuming that u(t) is a one-to-one function.
 * Two warping methods exist to warp g(t) to time u.
 * The method simpleWarp involves applying an anti-aliasing filter to g(t) and 
 * warping the anti-aliased g(t) to time u. If u'(t) is constant for all t, then the 
 * anti-aliasing filter does not eliminate information in g(t).
 * If u'(t) varies with t, then the anti-aliasing
 * filter will use the maximum value of u'(t) to ensure that g(t) is not aliased. This 
 * results in a loss of information in g(t) because not all t in g(t) need the maximum amount 
 * of anti-aliasing applied. In the applySimpleS method, the anti-alias filter is created using
 * the maximum u'(t), this filter is then applied to g(t), and then sinc interpolation 
 * is used to resample g(t) at times u. In order for impulses to be warped without wavelet 
 * distortion, g(u(t)) needs to be scaled by u'(t). 
 * The anti-alias filter and sinc interpolation can be 
 * represented by the linear operators L and W, respectively. The scaling will be included in the
 * operator W. Applying the simpleWarp method to
 * g can simply be represented as the linear operator S, which is WLg.
 * The other method, applyS, does not eliminate information in g(t) when u'(t) is varying. 
 * The steps and corresponding linear operators in the method warp are 
 * upsampling (U), warping (W), anti-alias filtering (L), and downsampling (D). U, L, and D will
 * be determined by the maximum u'(t).
 * The method warp is represented by S, which is DLWU.
 * 
 * @author Chris Graziano, Colorado School of Mines
 * @version 2014.09.08
 */

public class Warper {
  /**
   * Sets the maximum amount of warping allowed. If a u'(t) at some time is above this maximum 
   * amount of warp (_rmax), this u'(t) value will not influence the anti-alias filter. The default
   * _rmax is 10.0.
   */
  public void setRmax(float rmax) {
     _rmax = rmax;
  }

  public void setScale(boolean scale) {
    _scale = scale;
  }
  
  /**
   * If plotting is set to true, then the amplitude spectrums of uu, ug, wug, lwug, and dlwug 
   * and the traces uu, ug, wug, lwug, and dlwug will be plotted. This will only enable 
   * plotting within the applyS method that has single arrays of floats as inputs.
   * This method is for testing only and SHOULD NOT BE INCLUDED IN THE FINAL SHIPMENT OF 
   * CODE. The default plotting setting is no plotting will occur.
   */
  public void plotAmplitudeSpectrums(boolean plot) {
    _plot = plot;
  }

  /**
   * Returns uu, wug, lwug, and dlwug.
   */
  public float[][] getWarpingStages() {
    return new float[][]{_ug,_wug,_lwug,_dlwug};
  }

  public float getSamplingIntervalScaleFactor() {
    return _r;
  }


   /**
   * A method to warp sequences without aliasing.
   * There are four steps to this method:
   * (1) Upsample sequences g and u (U).
   * (2) Warp the upsampled sequence using sinc interpolation (W).
   * (3) Apply an anti-aliasing filter to the previously warped sequence (L)
   * (4) Downsample the anti-aliased sequence (D).
   * @param r the largest u'(t) in u. This value needs to be an integer represented as a 
              float. For example, if r is 3.01 or 3.99, the input for r is 4.0.
   * @param u the relationship between the times that events occur at in g and the times 
   *          where these events will occur at after warping.
   * @param g the sequence to be warped
   * @return DLSUg
   */
  public float[] applyS(float[] u, float[] g) {
    float r = (float)((int)(squeezing(1.0f,u)+0.999f));//Round up to the nearest integer.
    _r = r;
//r = 3.0f;
//trace("r = "+r);
    return applyS(r,u,g);
  }
  public float[] applyS(float[] u, float[] uu, float[] g) {
    float r = (float)((int)(squeezing(1.0f,u)+0.999f));//Round up to the nearest integer.
    _r = r;
//r = 3.0f;
//trace("r = "+r);
    return applyS(r,u,g);
  }

  public float[] applyS(float r,float[] u, float[] g) {
    int ng = g.length;
    int nu = u.length;
    
    if (r>1.0f) {
      int nug = (int)(r*(ng-1)+1);//Number of samples in upsampled sequence.
      int nuu = (int)(r*(nu-1)+1);//Number of samples in upsampled sequence.
      float dtug = 1.0f/r;

      float[] uu = upSampleLinear(nu,1.0f,0.0f,u,nuu,dtug,0.0f);//upsampled u
      float[] ug = upSample(ng,1.0f,0.0f,g,nug,dtug,0.0f);//upsampled g
      float[] wug  = applyW(dtug,uu,ug);//warp upsampled g
      float[] lwug = applyL(r,uu,wug);//anti-alias warped, upsampled g.
      float[] dlwug = subSample(r,nu,lwug);
      _ug = ug;
      _wug = wug;
      _lwug = lwug;
      _dlwug = dlwug;
      
      return dlwug;
    }

    else {
      float[] wg  = applyW(u,g);
      return wg;
    }
  }

  public float[][] applyS(float[][] u, float[][] g) {
    int n2 = u.length;
    int n1 = u[0].length;
    float[][] sg = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      sg[i2] = applyS(u[i2],g[i2]);
    }
    return sg;
  }
  public float[][][] applyS(final float[][][] u, final float[][][] g) {
    final int n3 = u.length;
    final int n2 = u[0].length;
    final int n1 = u[0][0].length;
    final float[][][] sg = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() 
    {
      public void compute(int i3) 
      {

    //for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        sg[i3][i2] = applyS(u[i3][i2],g[i3][i2]);
      }
      }
    });
    return sg;
  }

  /**
   * A method to warp sequences without aliasing.
   * There are four steps to this method:
   * (1) Upsample sequences g and u (U).
   * (2) Warp the upsampled sequence using sinc interpolation (W).
   * (3) Apply an anti-aliasing filter to the previously warped sequence (L)
   * (4) Downsample the anti-aliased sequence (D).
   * @param r the largest u'(t) in u. This value needs to be an integer represented as a 
              float. For example, if r is 3.01 or 3.99, the input for r is 4.0.
   * @param nu1 the number of samples in the fast dimension of u before being upsampled.
   * @param uu the relationship between the times that events occur at in g and the times 
   *            where these events will occur at after warping. This sequence has 
               already been upsampled.
   * @param g the sequence to be warped
   * @return DLSUg
   */
  public float[][] applyS(float r, int nu1, float[][] uu, float[][] g) {
    int ng2 = g.length;
    int ng1 = g[0].length;
    int nuu1 = uu[0].length;
    //int nu1 = (int)((nuu1+1)/r+1);//length of u before it was upsampled to create uu. 
    float[][] dlsug = new float[ng2][nu1];
    Stopwatch sw1 = new Stopwatch();
    sw1.start();
    
    if (r>1.0f) {
      int nug = (int)(r*(ng1-1)+1);//Number of samples in upsampled sequence.
      int nuu = (int)(r*(nu1-1)+1);//Number of samples in upsampled sequence.
      float dtug = 1.0f/r;
      java.util.Random ra = new java.util.Random(5);
      //sw.start();
      //float[][] uu = upSampleLinear(nu1,1.0f,0.0f,u,nuu,dtug,0.0f);//upsampled u
      //float[][] uu = fillfloat(0.0f,nug,ng2);
      //rand(ra,uu);
      //sw.stop();
      //System.out.println("upSample u time = "+sw.time());
      //sw.restart();
      float[][] ug = upSample(ng1,1.0f,0.0f,g,nug,dtug,0.0f);
      //float[][] ug = fillfloat(0.0f,nug,ng2);
      //rand(ra,ug);
      //sw.stop();
      //System.out.println("upSample g time = "+sw.time());
      //sw.restart();
      float[][] sug  = applyW(dtug,uu,ug);//warp upsampled g
      //float[][] sug = fillfloat(0.0f,nug,ng2);
      //rand(ra,sug);
      //sw.stop();
      //System.out.println("warp time = "+sw.time());
      //sw.restart();
      float[][] lsug = applyL(r,uu,sug);//anti-alias warped, upsampled g.
      //float[][] lsug = fillfloat(0.0f,nug,ng2);
      //rand(ra,lsug);
      //sw.stop();
      //System.out.println("anti alias time = "+sw.time());
      //sw.restart();
      dlsug = subSample(r,nu1,lsug);
      //dlsug = fillfloat(0.0f,nu1,ng2);
      //rand(ra,dlsug);

      //sw.stop();
      //System.out.println("subSample time = "+sw.time());
      //sw1.stop();
      //System.out.println("time to complete one W= "+sw1.time());
      return dlsug;
    }

    else {
      float[][] u = subSample(r,nu1,uu);
      float[][] sg  = applyW(u,g);
      return sg;
    }
  }

  /**
   * Upsamples a sequence from one sampling (ntf, dtf, ftf) to another
   * sampling (ntg, dtg, ftg). Note that dtg should be a multiple of dtf. 
   * This method uses sinc interpolation to upsample f.
   * @param ntf the number of samples in f
   * @param dtf the sampling interval of f
   * @param ftf the time of the first sample in f
   * @param f the sequence to be upsampled
   * @param ntg the number of samples in the upsampled sequence g (ntg = (int)((ntf+1)/r+)),
                where r is a whole number and the amount you want f upsampled by r = dtg/dtf
   * @param dtg the new sampling interval (should be a multiple of dtf)
   * @param ftg the time of the first sample in f
   */
  public float[] upSample(
    int ntf, float dtf, float ftf, float[] f, 
    int ntg, float dtg, float ftg)  
  {
      float[] g = new float[ntg];
      int r = (int)(dtf/dtg);
      float odr = 1.0f/r;
      float[][] gtemp = new float[r][ntg];
      gtemp[0] = f;
      for (int i=1; i<r; ++i)
        _si.interpolate(ntf,1.0,0.0,f,ntf,1.0,odr*i,gtemp[i]);
      for (int it=0; it<ntf; ++it) {
        g[it*r] = gtemp[0][it];
      }
      for (int i=1; i<r; ++i) {
        for (int it=0; it<ntf-1; ++it) {
          g[i+it*r] = gtemp[i][it];
        }
      }
      return g;
  }
  public float[][] upSample(
      int ntf, float dtf, float ftf, float[][] f, 
      int ntg, float dtg, float ftg)  
  {
      int n2 = f.length;
      int n1 = f[0].length;
      float[][] uf = new float[n2][n1];
      for (int i=0; i<n2; ++i) 
        uf[i] = upSample(ntf,dtf,ftf,f[i],ntg,dtg,ftg);
     
      return uf;
  }

  /**
   * Upsamples a sequence from one sampling (ntf, dtf, ftf) to another
   * sampling (ntg, dtg, ftg). Note that dtg should be a multiple of dtf. 
   * This method uses linear interpolation to upsample f.
   * @param ntf the number of samples in f
   * @param dtf the sampling interval of f
   * @param ftf the time of the first sample in f
   * @param f the sequence to be upsampled
   * @param ntg the number of samples in the upsampled sequence g (ntg = (int)((ntf+1)/r+)),
                where r is a whole number and the amount you want f upsampled by r = dtg/dtf
   * @param dtg the new sampling interval (should be a multiple of dtf)
   * @param ftg the time of the first sample in f
   */
  public float[] upSampleLinear(
      int ntf, float dtf, float ftf, float[] f, 
      int ntg, float dtg, float ftg)  
  {
      float[] g = new float[ntg];
      _li.setUniform(ntf,dtf,ftf,f);
      _li.interpolate(ntg,dtg,ftg,g);
      return g;
  }
  public float[][] upSampleLinear(
      int ntf, float dtf, float ftf, float[][] f, 
      int ntg, float dtg, float ftg)  
  {
      int n2 = f.length;
      int n1 = f[0].length;
      float[][] g = new float[n2][n1];

      for (int i=0; i<n2; ++i) {
        g[i] = upSampleLinear(ntf,dtf,ftf,f[i],ntg,dtg,ftg);
      }
      return g;
  }
  
  /**
   * Upsamples a sequence by a factor of r.
   * This method uses linear interpolation to upsample f.
   * @param ntf the number of samples in f
   * @param dtf the sampling interval of f
   * @param ftf the time of the first sample in f
   * @param f the sequence to be upsampled
   * @param r the amount to upsample f by (should be a whole number)
   */

  public float[][] upSampleLinear(int ntf, float dtf, float ftf, float[][] f, float r) {
    int n2 = f.length;
    int n1 = f[0].length;
    int ntg = (int)(r*(ntf-1)+1);//Number of samples in upsampled sequence.
    float dtg = 1.0f/r;
    float ftg = 0.0f;
    float[][] g = new float[n2][n1];

    for (int i=0; i<n2; ++i) {
      g[i] = upSampleLinear(ntf,dtf,ftf,f[i],ntg,dtg,ftg);
    }
    return g;

  }



  /**
   * Applies the warping operator W.
   * @param u array of new times to warp the current times of f to.
   * @param f array with input sequence f(t).
   * @return array with warped output sequence.
   */
  public float[] applyW(float[] u, float[] f) {
    return warp(u,f);
  }
  public float[] applyW(float dt, float[] u, float[] f) {
    return warp(dt,u,f);
  }
  public float[][] applyW(float[][] u, float[][] f) {
    return warp(u,f);
  }
  public float[][] applyW(float dt, float[][] u, float[][] f) {
    return warp(dt,u,f);
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
    return aaf(1.0f,_rmax,u,f);
  }
  public float[] applyL(float r, float[] u, float[] f) {
    return aaf(r,_rmax,u,f);
  }
  public float[][] applyL(float r, float[][] u, float[][] f) {
    return aaf(r,_rmax,u,f);
  }
  public float[][] applyL(BandPassFilter bpf, float r, float[][] u, float[][] f) {
    return aaf(bpf,r,u,f);
  }

  /** 
   * Subsamples the sequence f by a factor of r.
   * This method assumes that r is an whole number because
   * every sample that is not the rth sample is thrown away.
   */
  public float[] subSample(float dtf, float[] f, int ntg, float dtg) {
    float r = dtf/dtg;
    return subSample(r,ntg,f);
  }
  public float[] subSample(float r, int ntg, float[] f) {
    float[] g = new float[ntg];
    int ntgm1 = ntg-1;
    for (int i=0; i<ntgm1; ++i) {
      int ir = (int)(i*r);
      g[i] = f[ir];
    }
    return g;
  }
  public float[][] subSample(float r, int ntg, float[][] f) {
    int nf2 = f.length;
    float[][] g = new float[nf2][ntg];
    for (int i2=0; i2<nf2; ++i2)
      g[i2] = subSample(r,ntg,f[i2]);
    return g;
  }

   /**
   * Returns the largest squeezing r(t) = u'(t) not greater than rmax.
   * If less than or equal to one, then no squeezing is implied by u(t).
   * Assumes the sampling interval is 1.
   */
  public float squeezing(float dt, float[] u) {
    int nt = u.length;
    int itlo = 1;
    int ithi = nt-1;
    float r = 0.0f;
    for (int it=itlo; it<=ithi; ++it) {
      float du = (u[it]-u[it-1])/dt;
      if (r<du)
        r = du;
    }
    return min(r,_rmax);
  }
  public float squeezing(float dt, float[][] u) {
    int n = u.length;
    float r = 0.0f;
    for (int i=0; i<n; ++i)
      r = max(r,squeezing(dt,u[i]));
    return r;
  }

  /**
   * Returns y(t) = x(u(t)).
   * For paper purposes. Not meant to be used.
   */
  public float[] applyOldS(float[] u, float[] x) {
    float r = (float)((int)(squeezing(1.0f,u)+0.999f));//Round up to the nearest integer.
    float[] lx = applyL(r,u,x);
    float[] wlx = applyW(u,lx);
    return wlx;
  }



  
  ///////////////////////////////////////////////////////////////////////////
  // private
  private float _rmax = 10.0f;
  private float _r = 0.0f;
  private boolean _plot = false;
  private boolean _scale = true;
  private static final SincInterp _si = SincInterp.fromErrorAndFrequency(0.01,0.40);
  private static final LinearInterpolator _li = new LinearInterpolator();
  private float[] _ug, _wug, _lwug, _dlwug;
  

  /**
   * Returns y(t) = x(u(t)).
   */
  private float[] warp(float[] u, float[] x) {
    int nt = u.length;
    float[] y = new float[nt];
    _si.interpolate(x.length,1.0,0.0,x,nt,u,y);
    if (_scale==true) {
      y[0] *= u[1]-u[0];
      for (int it=1; it<nt; ++it)
        y[it] *= u[it]-u[it-1];
    }
    return y;
  }
  private float[] warp(float dt,float[] u, float[] x) {
    int nu = u.length;
    int nx = x.length;
    float[] y = new float[nu];
    float odt = 1.0f/dt;
    _si.interpolate(nx,dt,0.0,x,nu,u,y);
    if (_scale==true) {
      y[0] *= (u[1]-u[0])*odt;
      for (int it=1; it<nu; ++it) {
        y[it] *= (u[it]-u[it-1])*odt;
      }
    }
    return y;
  }
  private float[][] warp(float[][] u, float[][] x) {
    int n = u.length;
    float[][] y = new float[n][];
    for (int i=0; i<n; ++i)
      y[i] = warp(u[i],x[i]);
    return y;
  }
  private float[][] warp(float dt, float[][] u, float[][] x) {
    final int nu2 = u.length;
    final int nu1 = u[0].length;
    final int nx1 = x[0].length;
    final float[][] y = new float[nu2][nu1];
    final float odt = 1.0f/dt;
    for (int i=0; i<nu2; ++i) {
      _si.interpolate(nx1,dt,0.0,x[i],nu1,u[i],y[i]);
      if (_scale==true) {
        y[i][0] *= (u[i][1]-u[i][0])*odt;
        for (int it=1; it<nu1; ++it) {
          y[i][it] *= (u[i][it]-u[i][it-1])*odt;
        }
      }
    }
    return y;

  }
  /**
   * If necessary, applies an anti-alias filter to the sequence x(t).
   * An anti-alias filter is necessary if the warping includes squeezing.
   * Assumes the sampling interval is one.
   */
  private float[] aaf(float r, float rmax, float[] u, float[] x) {
    int nt = x.length;
    if (r>1.0) {
      float[] y = new float[nt];
      float width = 0.20f/r;
      BandPassFilter bpf = new BandPassFilter(0.0,0.5/r,width,0.01);
      bpf.apply(x,y);
      return y;
    } else {
      return copy(x);
    }
  }

  public float[][] aaf(float r, float rmax, float[][] u, float[][] x) {
    if (r>1.0) {
      int nx = x.length;
      int nt = x[0].length;
      float[][] y = new float[nx][nt];
      float width = 0.20f/r;
      BandPassFilter bpf = new BandPassFilter(0.0,0.5/r,width,0.01);
      for (int ix=0; ix<nx; ++ix)
        bpf.apply(x[ix],y[ix]);
      return y;
    } else {
      return copy(x);
    }
  }
  private float[][] aaf(BandPassFilter bpf, float r, float[][] u, float[][] x) {
    int nx2 = x.length;
    int nx1 = x[0].length;
    float width = 0.02f/r;
    float[][] y = new float[nx2][nx1];
    for (int ix=0; ix<nx2; ++ix)
      bpf.apply(x[ix],y[ix]);
    return y;
  }


  private static void trace(String s) {
    System.out.println(s);
  }







}


