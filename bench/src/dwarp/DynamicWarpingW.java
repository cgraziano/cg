/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dwarp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.lapack.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;
import dwarp.Warp;



public class DynamicWarpingW {

  public DynamicWarpingW(Warp warp) {
    SincInterp si = new SincInterp();
    _si = si;
    _warp = warp;
  }

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
* Sets the min-max range of frequencies in wavelet.
* @param fmin minimum frequency, in cycles/sample.
* @param fmax maximum frequency, in cycles/sample.
*/
  public void setFrequencyRange(double fmin, double fmax) {
    _bpf = new BandPassFilter(fmin,fmax,0.05,0.01);
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
  

  public float[] getInverseAWarpNoOSamp(
      int na, int ka,
      Sampling stf, Sampling stg,
      float awarp, float[] f, float[] g)
  {
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");
    int nt = f.length;
    double dt = stf.getDelta();
    
    float[][] d = computeDifferences(na,ka,_bpf,stf,stg,awarp,f,g);

    // The matrix C and right-hand-side vector b, for Ca = b.
    // For zero lag, we have a0 = a[-ka] = 1, so that
    // only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags.
    int ma = na-1;
    DMatrix c = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jc=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double cij = dot(d[ia],d[ja]);
        c.set(ic,jc,cij);
        ++jc;
      }
      c.set(ic,ic,c.get(ic,ic)*_sfac);
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    System.out.println("c=\n"+c);
    System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of C.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // lag 0, so a0 = 1
      } else {
        aa[ia] = (float)a.get(ic,0);
        ++ic;
      }
    }
    return aa;
  }


  /**
* For each lag of the inverse wavelet, computes differences between
* NMO-corrected gathers and the stacked-and-replicated versions of those
* gathers.
*/
  public float[][] computeDifferences(
    int na, int ka, BandPassFilter bpf,
    Sampling stf, Sampling stg,
    float awarp, float[] f, float[] g)
  {
    int nt = f.length;
    float[][] d = new float[na][nt];
    int ntg = stg.getCount();
    double dtg = stg.getDelta();
    double ftg = stg.getFirst();
    double t = dtg*(ntg-1);//time in trace
    //for viewing purposes only
    _df = new float[na][nt];
    _bdf = new float[na][ntg];
    _dg = new float[na][ntg];
    _sdg = new float[na][ntg];
    _bsdg = new float[na][ntg];
    _d = new float[na][nt];
    _stf = stf;
    _stg = stg;
    ////////////////////////
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[] df = delay(lag,f);
      float[] dg = delay(lag,g);
      float[] sdg = _warp.apply(stg,awarp,dg);
      _df[ia] = copy(df);
      _dg[ia] = copy(dg);
      _sdg[ia] = copy(sdg);
      
      //PlottingPlottingPlotting//
      String title = "df "+lag;
      int amax = 1;
      int tmin = 0;
      int tmax = 3;
      plotTrace(stf,df,tmin,tmax,amax,title);
      //PlottingPlottingPlotting//
      title = "sdg "+lag;
      amax = 1;
      tmin = 0;
      tmax = 3;
      plotTrace(stf,sdg,tmin,tmax,amax,title);

      for (int it=0; it<nt; ++it) {
        d[ia][it] = sdg[it]-df[it];
      }
      //PlottingPlottingPlotting//
      title = "d "+lag;
      amax = 1;
      tmin = 0;
      tmax = 3;
      plotTrace(stf,d[ia],tmin,tmax,amax,title);
      
      if (bpf!=null) {
        bpf.apply(d[ia],d[ia]);
      }
    }

    return d;
  }

  /*public float[] getInverseAWarp(
      int na, int ka,
      Sampling stf, Sampling stg,
      float[] shifts, float[] f, float[] g)
  {
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");
    int nt = f.length;
    double dt = stf.getDelta();

    float[][] d = computeDifferences(na,ka,_bpf,stf,stg,shifts,f,g);

    // The matrix C and right-hand-side vector b, for Ca = b.
    // For zero lag, we have a0 = a[-ka] = 1, so that
    // only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags.
    int ma = na-1;
    DMatrix c = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jc=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double cij = dot(d[ia],d[ja]);
        c.set(ic,jc,cij);
        ++jc;
      }
      c.set(ic,ic,c.get(ic,ic)*_sfac);
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    //System.out.println("c=\n"+c);
    //System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of C.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // lag 0, so a0 = 1
      } else {
        aa[ia] = (float)a.get(ic,0);
        ++ic;
      }
    }
    return aa;
  }
  */

  /*
  public float[] getInverseAWarp(
      int na, int ka,
      Sampling stf, Sampling stg, Sampling sx,
      float[][] shifts, float[][] f, float[][] g)
  {
    Check.argument(-na<ka,"-na<ka");
    Check.argument(ka<=0,"ka<=0");
    int nt = f[0].length;
    int nx = f.length;
    double dt = stf.getDelta();

    float[][][] d = computeDifferences(na,ka,_bpf,stf,stg,sx,shifts,f,g);

    // The matrix C and right-hand-side vector b, for Ca = b.
    // For zero lag, we have a0 = a[-ka] = 1, so that
    // only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags.
    int ma = na-1;
    DMatrix c = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jc=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double cij = dot(d[ia],d[ja]);
        c.set(ic,jc,cij);
        ++jc;
      }
      c.set(ic,ic,c.get(ic,ic)*_sfac);
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    System.out.println("c=\n"+c);
    System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of C.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // lag 0, so a0 = 1
      } else {
        aa[ia] = (float)a.get(ic,0);
        ++ic;
      }
    }
    return aa;
  }
  */
/**
* For each lag of the inverse wavelet, computes differences between
* NMO-corrected gathers and the stacked-and-replicated versions of those
* gathers.
*/
/*
  public float[][] computeDifferencesOSamp(
    int na, int ka, BandPassFilter bpf,
    Sampling stf, Sampling stg,
    float awarp, float[] f, float[] g)
  {
    int nt = f.length;
    float[][] d = new float[na][nt];
    int ntg = stg.getCount();
    double dtg = stg.getDelta();
    double ftg = stg.getFirst();
    double t = dtg*(ntg-1);//time in trace
    double dtgO = dtg*.5f;
    int ntgO = (int)(t/dtgO + 1.0);
    Sampling stgO = new Sampling(ntgO,dtgO,ftg);
    //for viewing purposes only
    _df = new float[na][nt];
    _dg = new float[na][ntg];
    _dgO = new float[na][ntgO];
    _sdgO = new float[na][ntgO];
    _sdg = new float[na][ntg];
    _d = new float[na][nt];
    _stf = stf;
    _stg = stg;
    _stgO = stgO;
    ////////////////////////
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[] df = delay(lag,f);
      float[] dg = delay(lag,g);
      float[] dgO = new float[ntgO];
      _si.interpolate(ntg,dtg,0.0,dg,ntgO,dtgO,0.0,dgO);
      float[] sdgO = _warp.apply(stgO,awarp,dgO);
      float[] sdg = new float[ntg];
      _si.interpolate(ntgO,dtgO,0.0,sdgO,ntg,dtg,0.0,sdg);
      
      _df[ia] = copy(df);
      _dg[ia] = copy(dg);
      _dgO[ia] = copy(dgO);
      _sdgO[ia] = copy(sdgO);
      _sdg[ia] = copy(sdg);

      if (bpf!=null) {
        bpf.apply(df,df);
        bpf.apply(sdg,sdg);
      }
      for (int it=0; it<nt; ++it) {
        d[ia][it] = sdg[it]-df[it];
      }
      _d[ia] = copy(d[ia]);
    }
    return d;
  }
  */




  /**
* For each lag of the inverse wavelet, computes differences between
* NMO-corrected gathers and the stacked-and-replicated versions of those
* gathers.
*/
/*
  public float[][] computeDifferences(
    int na, int ka, BandPassFilter bpf,
    Sampling stf, Sampling stg,
    float[] shifts, float[] f, float[] g)
  {
    int nt = f.length;
    float[][] d = new float[na][nt];
    int ntg = stg.getCount();
    double dtg = stg.getDelta();
    double ftg = stg.getFirst();
    double t = dtg*(ntg-1);//time in trace
    double dtg2 = dtg*.5f;
    int ntg2 = (int)(t/dtg2 + 1.0);
    Sampling stg2 = new Sampling(ntg2,dtg2,ftg);
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[] df = delay(lag,f);
      float[] dg = delay(lag,g);
      float[] dg2 = new float[ntg2];
      float[] shifts2 = new float[ntg2];
      _si.interpolate(ntg,dtg,0.0,dg,ntg2,dtg2,0.0,dg2);
      _si.interpolate(ntg,dtg,0.0,shifts,ntg2,dtg2,0.0,shifts2);
      float[] sdg2 = _warp.apply(stg2,shifts2,dg2);
      float[] sdg = new float[ntg];
      _si.interpolate(ntg2,dtg2,0.0,sdg2,ntg,dtg,0.0,sdg);
      
      

      if (bpf!=null) {
        bpf.apply(df,df);
        bpf.apply(sdg,sdg);
      }
      for (int it=0; it<nt; ++it) {
        d[ia][it] = sdg[it]-df[it];
      }
    }
    return d;
  }
  */


  /**
* For each lag of the inverse wavelet, computes differences between
* NMO-corrected gathers and the stacked-and-replicated versions of those
* gathers.
*/
/*
  public float[][][] computeDifferences(
    int na, int ka, BandPassFilter bpf,
    Sampling stf, Sampling stg, Sampling sx,
    float[][] shifts, float[][] f, float[][] g)
  {
    int nt = f[0].length;
    int nx = f.length;
    float[][][] d = new float[na][nx][nt];
    for (int ia=0,lag=ka; ia<na; ++ia,++lag) {
      float[][] df = delay(lag,f);
      float[][] dg = delay(lag,g);
      float[][] sdg = applyShifts(stg,dg,shifts);
      if (bpf!=null) {
        for (int ix=0; ix<nx; ++ix) {
          bpf.apply(df[ix],df[ix]);
          bpf.apply(sdg[ix],sdg[ix]);
        }
      }
      for (int ix=0; ix<nx; ++ix) {
        for (int it=0; it<nt; ++it) {
          d[ia][ix][it] = sdg[ix][it]-df[ix][it];
        }
      }
    }
    return d;
  }
  */

  /**
* Returns differences between NMO-corrected gathers and stacks.
* @param na number of samples in the inverse wavelet a.
* @param ka the sample index for a[0].
* @param st time sampling.
* @param sx offset sampling.
* @param vnmo array[nt] of NMO velocities.
* @param f array[nx][nt] with input CMP gather.
* @return array[na][nx][nt] of difference gathers.
*/
  /*public float[][][] getDifferenceGathers(
    int na, int ka,
    Sampling stf, Sampling stg, Sampling sx,
    float[][] shifts, float[][] f, float[][] g)
  {
    float[][][] d = computeDifferences(na,ka,_bpf,stf,stg,sx,shifts,f,g);
    for (int ia=0,lag=ka; ia<na; ++ia,++lag)
      d[ia] = delay(-lag,d[ia]);
    //for (int ia=1; ia<na; ++ia)
    // d[ia] = sub(d[ia],d[0]);
    return d;
  }
  */

  /**
* Returns differences between NMO-corrected gathers and stacks.
* @param na number of samples in the inverse wavelet a.
* @param ka the sample index for a[0].
* @param st time sampling.
* @param sx offset sampling.
* @param vnmo array[nt] of NMO velocities.
* @param f array[nx][nt] with input CMP gather.
* @return array[na][nx][nt] of difference gathers.
*/
/*
  public float[][] getDifferenceGathers(
    int na, int ka,
    Sampling stf, Sampling stg,
    float[] shifts, float[] f, float[] g)
  {
    float[][] d = computeDifferences(na,ka,_bpf,stf,stg,shifts,f,g);
    for (int ia=0,lag=ka; ia<na; ++ia,++lag)
      d[ia] = delay(-lag,d[ia]);
    //for (int ia=1; ia<na; ++ia)
    // d[ia] = sub(d[ia],d[0]);
    return d;
  }
  */

  /**
* Returns uniformly sampled warped sequence h(x1) = g(x1+u(x1)).
* @param sg sampling of the sequence g to be warped.(PS Image)
* @param g array for the sequence g to be warped.
* @param u array of shifts.
* @return array for the warped sequence h.
*/
  public float[] applyShifts(Sampling sg, float[] g, float[] u) {
    Sampling s1 = sg;
    int ng = sg.getCount();
    int n1 = s1.getCount();
    float[] h = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      double x1 = s1.getValue(i1)+u[i1];
      System.out.println("x1 = "+x1);
      h[i1] = _si.interpolate(sg,g,x1);
    }
    return h;
  }

  /**
* Returns uniformly sampled warped image h(x1,x2) = g(x1+u(x1,x2),x2).
* @param sg sampling of the sequence g to be warped.
* @param g array for the sequence g to be warped.
* @param u array of shifts.
* @return array for the warped sequence h.
*/
  public float[][] applyShifts(Sampling sg, float[][] g, float[][] u) {
    int n2 = g.length;
    float[][] h = new float[n2][];
    for (int i2=0; i2<n2; ++i2)
      h[i2] = applyShifts(sg,g[i2],u[i2]);
    return h;
  }

  /**
* Returns uniformly sampled warped image h(x1,x2,x3) = g(x1+u(x1,x2,x3),x2,x3).
* @param sg sampling of the sequence g to be warped.
* @param g array for the sequence g to be warped.
* @param u array of shifts.
* @return array for the warped sequence h.
*/
  public float[][][] applyShifts(Sampling sg,
      float[][][] g, float[][][] u) {
    int n3 = g.length;
    int n2 = g[0].length;
    float[][][] h = new float[n3][n2][];

    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        h[i3][i2] = applyShifts(sg,g[i3][i2],u[i3][i2]);
    return h;
  }

  public static float[] getAmplitudes(
    Sampling st, float smax, float[] f, float[] t)
  {
    int nt = st.getCount();
    float dt = (float)st.getDelta();
    float ft = (float)st.getFirst();
    float odt = 1.0f/dt;
    float dtmin = dt/smax;
    System.out.println("dtmin = "+dtmin);
    float[] a = new float[nt];


    // Time of first non-zero input sample.
    int nz = countLeadingZeros(f);
    float tnz = ft+nz*dt;

    // Number of leading zeros in output. A leading output sample is zero
    // if either (1) the corresponding input samples and all prior input
    // samples are zero, or (2) NMO stretch would exceed the maximum.
    nz = 0;
    if (t[0]<tnz || t[1]-t[0]<dtmin)
      ++nz;
    for (int it=1; it<nt; ++it) {
      if (t[it]<tnz || t[it]-t[it-1]<dtmin)
        ++nz;
    }

    // Compute only the non-zero amplitudes. These amplitudes are simply the
    // inverse of NMO stretch.
    if (nz==0) {
      a[0] = (t[1]-t[0])*odt;
      ++nz;
    }
    for (int it=nz; it<nt; ++it)
      a[it] = (t[it]-t[it-1])*odt;
    return a;
  }


  public float getVariancePef(
    int na, int ka, float[] a, float[][] g)
  {
    float[][] ga = applyFilter(na,ka,a,g);
    return pow(rms(ga),2.0f);
  }

  public float getVarianceDww(
    int na, int ka, float[] a,
    Sampling stf, Sampling stg, Sampling sx,
    float[][] shifts, float[][] f, float[][] g)
  {
    float[][] bda = applyBSA(na,ka,a,stf,stg,sx,shifts,f,g);
    return pow(rms(bda),2.0f);
  }

  /*public float getNormalizedVarianceDww(
int na, int ka, float[] a,
Sampling st, Sampling sx, float[] vnmo, float[][] f)
{
float[][] g = applyBNmoA(na,ka,a,st,sx,vnmo,f);
float[][] r = _nmo.stackAndReplicate(g);
return pow(rms(sub(g,r))/rms(g),2.0f);
}
*/

 public float[][] applyHSA(
    int na, int ka, float[] a,
    int nh, int kh, float[] h,
    Sampling stf, Sampling stg, Sampling sx,
    float[][] shifts, float[][] f, float[][] g)
  {
    int nt = f[0].length;
    int nx = f.length;
    float[][] ga = applyFilter(na,ka,a,g);
    float[][] sga = applyShifts(stg,g,shifts);
    return applyFilter(nh,kh,h,sga);
  }

  public float[][] applyBSA(
    int na, int ka, float[] a,
    Sampling stf, Sampling stg, Sampling sx,
    float[][] shifts, float[][] f, float[][] g)
  {
    int nx = sx.getCount();
    float[][] warp = applyShifts(stg,g,shifts);
    float[][] d = sub(warp,f);
    float[][] da = applyFilter(na,ka,a,d);
    float[][] bda = da;
    if (_bpf!=null) {
      for (int ix=0; ix<nx; ++ix)
        _bpf.apply(da[ix],bda[ix]);
    }
    return bda;
  }

  /**
* Returns inverse wavelet a estimated via PEF of gather.
* @param na number of samples in the inverse wavelet a.
* @param ka the sample index for a[0].
* @param f array[nx][nt] for CMP gather.
* @return array of coefficients for the inverse wavelet a.
*/
  public float[] getInverseAPef(int na, int ka, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;

    // CMP gather for different time shifts (band-pass filtered?).
    float[][][] d = new float[na][nx][nt];
    for (int ia=0; ia<na; ++ia) {
      d[ia] = delay(ka+ia,f);
      /* band-pass filter causes unstable estimates of wavelet
if (_bpf!=null) {
for (int ix=0; ix<nx; ++ix)
_bpf.apply(d[ia][ix],d[ia][ix]);
}
*/
    }

    // The matrix C and right-hand-side vector b, for Ca = b.
    // For zero lag, we have a0 = a[-ka] = 1, so that
    // only na-1 coefficients of a are unknown;
    // the unknown coefficients are those corresponding to non-zero lags.
    int ma = na-1;
    DMatrix c = new DMatrix(ma,ma);
    DMatrix b = new DMatrix(ma,1);
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) continue; // skip lag zero, because a0 = 1
      for (int ja=0,jc=0; ja<na; ++ja) {
        if (ja==-ka) continue; // skip lag zero, because a0 = 1
        double cij = dot(d[ia],d[ja]);
        c.set(ic,jc,cij);
        ++jc;
      }
      c.set(ic,ic,c.get(ic,ic)*_sfac);
      double bi = -dot(d[ia],d[-ka]);
      b.set(ic,0,bi);
      ++ic;
    }
    System.out.println("c=\n"+c);
    System.out.println("b=\n"+b);

    // Solve for inverse filter a using Cholesky decomposition of C.
    DMatrixChd chd = new DMatrixChd(c);
    DMatrix a = chd.solve(b);
    float[] aa = new float[na];
    for (int ia=0,ic=0; ia<na; ++ia) {
      if (ia==-ka) {
        aa[ia] = 1.0f; // lag 0, so a0 = 1
      } else {
        aa[ia] = (float)a.get(ic,0);
        ++ic;
      }
    }
    return aa;
  }

  /**
* Computes the rms statistic for a specified vertical range in
* the gather.
* @param st time sampling of the gather
* @param sx space sampling of the gather
* @param t1 time to begin rms computation (inclusive)
* @param t2 time to end rms computation (exclusive)
*/
  public float rms(Sampling st, Sampling sx,
      float t1, float t2, float[][] f) {
    int nt = st.getCount();
    int nx = sx.getCount();
    int it1 = st.indexOfNearest((double)(t1));
    int it2 = st.indexOfNearest((double)(t2));
    float sum = 0.0f;
    float ftx = 0.0f;
    int n = nx*(it2-it1+1);
    for (int ix=0; ix<nx; ++ix)
      for (int it=it1; it<=it2; ++it) {
        ftx = f[ix][it];
        sum += ftx*ftx;
      }
    float ms = sum/n;//ms is mean square
    return sqrt(ms);
  }

  /**
* rms of the entire gather
*/
  private float rms(float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    return (float)(sqrt(dot(f,f)/nx/nt));
  }

  /**
* Computes the mean for a specified vertical range in
* the gather.
* @param st time sampling of the gather
* @param sx space sampling of the gather
* @param t1 time to begin rms computation (inclusive)
* @param t2 time to end rms computation (exclusive)
*/
  public float mean(Sampling st, Sampling sx,
      float t1, float t2, float[][] f) {
    int nt = st.getCount();
    int nx = sx.getCount();
    int it1 = st.indexOfNearest((double)(t1));
    int it2 = st.indexOfNearest((double)(t2));
    float sum = 0.0f;
    float ftx = 0.0f;
    int n = nx*(it2-it1+1);
    for (int ix=0; ix<nx; ++ix)
      for (int it=it1; it<=it2; ++it) {
        sum += abs(f[ix][it]);
      }
    return sum/n;
  }

  public static void plotAmplitudeSpectrum(Sampling st, float[] p,
    boolean db,String title) {
    // Time sampling.
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();

    // Frequency sampling.
    int nfft = FftReal.nfftSmall(2*nt);
    int nf = nfft/2+1;
    double df = 1.0/(nfft*dt);
    double ff = 0.0;
    Sampling fs = new Sampling(nf,df,ff);
    float[] amp = computeAmplitudeSpectrum(st, fs, nfft, p, db);
    plotSpectrum(fs,amp,title);
  }

  // Computes the amplitude spectra for the specified signal.
  public static float[] computeAmplitudeSpectrum(Sampling st,
    Sampling fs, int nfft, float[] p, boolean db)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nf = fs.getCount();
    double df = fs.getDelta();
    double ff = fs.getFirst();
    // Real-to-complex fast Fourier transform.
    FftReal fft = new FftReal(nfft);
    float[] cf = new float[2*nf];
    copy(nt,p,cf);
    fft.realToComplex(-1,cf,cf);

    // Adjust phase for possibly non-zero time of first sample.
    float[] wft = rampfloat(0.0f,-2.0f*FLT_PI*(float)(df*ft),nf);
    cf = cmul(cf,cmplx(cos(wft),sin(wft)));

    // Amplitude spectrum, normalized.
    float[] af = cabs(cf);
    float amax = max(max(af),FLT_EPSILON);
    af = mul(1.0f/amax,af);
    if (db) {
      af = log10(af);
      af = mul(20.0f,af);
    }
    return af;
  }

  /*public void plotDifferencePlots(boolean bdf, boolean bdg, boolean bdgO,
boolean bsdgO, boolean bsdg, boolean d, float tmin, float tmax,
float amax)
{
int na = _df.length;
for (int ia=0; ia<na; ++ia) {
String title = "df lag = "+ia;
if (bdf) plotTrace(_stf,_df[ia],tmin,tmax,amax,title);
title = "dg lag = "+ia;
if (bdg) plotTrace(_stg,_dg[ia],tmin,tmax,amax,title);
title = "dgO lag = "+ia;
if (bdgO) plotTrace(_stgO,_dgO[ia],tmin,tmax,amax,title);
title = "sdgO lag = "+ia;
if (bsdgO) plotTrace(_stgO,_sdgO[ia],tmin,tmax,amax,title);
title = "sdg lag = "+ia;
if (bsdg) plotTrace(_stg,_sdg[ia],tmin,tmax,amax,title);
title = "d lag = "+ia;
if (d) plotTrace(_stf,_d[ia],tmin,tmax,amax,title);
}
}
*/

  public void plotDifferencePlots(boolean bdf, boolean bbdf, boolean bdg,
    boolean bsdg, boolean bbsdg, boolean d, float tmin, float tmax,
    float amax)
  {
    int na = _df.length;
    for (int ia=0; ia<na; ++ia) {
      String title = "df lag = "+ia;
      if (bdf) plotTrace(_stf,_df[ia],tmin,tmax,amax,title);
      title = "bdf lag = "+ia;
      if (bbdf) plotTrace(_stf,_bdf[ia],tmin,tmax,amax,title);
      title = "dg lag = "+ia;
      if (bdg) plotTrace(_stg,_dg[ia],tmin,tmax,amax,title);
      title = "sdg lag = "+ia;
      if (bsdg) plotTrace(_stg,_sdg[ia],tmin,tmax,amax,title);
      title = "bsdg lag = "+ia;
      if (bbsdg) plotTrace(_stg,_bsdg[ia],tmin,tmax,amax,title);
      title = "d lag = "+ia;
      if (d) plotTrace(_stf,_d[ia],tmin,tmax,amax,title);
    }
  }

  public void plotDifferenceSpectrums(boolean bdg,
    boolean bdgO, boolean bsdgO, boolean bsdg, int ia)
  {
    String title = "Amplitude dg lag = "+ia;
    if (bdg) plotAmplitudeSpectrum(_stg,_dg[ia],false,title);
    title = "Amplitude dgO lag = "+ia;
    if (bdgO) plotAmplitudeSpectrum(_stgO,_dgO[ia],false,title);
    title = "Amplitude sdgO lag = "+ia;
    if (bsdgO) plotAmplitudeSpectrum(_stgO,_sdgO[ia],false,title);
    title = "Amplitude sdg lag = "+ia;
    if (bsdg) plotAmplitudeSpectrum(_stg,_sdg[ia],false,title);
  }

  public void plotDifferenceSpectrums(boolean bdg, boolean bsdg, int ia)
  {
    String title = "Amplitude dg lag = "+ia;
    if (bdg) plotAmplitudeSpectrum(_stg,_dg[ia],false,title);
    title = "Amplitude sdg lag = "+ia;
    if (bsdg) plotAmplitudeSpectrum(_stg,_sdg[ia],false,title);
  }


  ///////////////////////////////////////////////////////////////
  //private
  private SincInterp _si;
  private Warp _warp;
  private double _sfac = 1.0;
  private int _itmin, _itmax;
  private BandPassFilter _bpf;
  private float[][] _df,_bdf,_dg,_dgO,_sdgO,_bsdg,_sdg,_d;
  private Sampling _stf,_stg,_stgO;

  

  /**
* For testing only
*/
  private static void plotTrace(Sampling st, float[] p,
      float tmin, float tmax, float amax, String title) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setVLabel("Time (s)");
    sp.setSize(400,750);
    sp.setVLimits(tmin,tmax);
    sp.setHLimits(-amax,amax);
    sp.addTitle(title);
    PointsView pv = sp.addPoints(st,p);
  }
  /**
* For testing only
*/
  private static void plotSpectrum(Sampling sf, float[] spec,
      String title) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.LOWER_LEFT);
    sp.setVLabel("Frequency");
    sp.setSize(750,400);
    sp.addTitle(title);
    PointsView pv = sp.addPoints(sf,spec);
  }
  /**
* For testing only
*/
  private static void plotGather(Sampling st,Sampling sx, float[][] p,
      float tmin, float tmax,float perc, String title) {
    SimplePlot sp = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
    sp.setHLabel("Offset (km)");
    sp.setVLabel("Time (s)");
    sp.setSize(400,750);
    sp.setVLimits(tmin,tmax);
    sp.addColorBar();
    PixelsView pv = sp.addPixels(st,sx,p);
    sp.addTitle(title);
    pv.setPercentiles(100-perc,perc);
  }

  /**
* Delays the CMP gather f by specified lag (which may be negative).
*/
  private static float[] delay(int lag, float[] f) {
    int nt = f.length;
    int itlo = max(0,lag); // 0 <= it-lag
    int ithi = min(nt,nt+lag); // it-lag < nt
    float[] g = new float[nt];
    for (int it=0; it<itlo; ++it)
      g[it] = 0.0f;
    for (int it=itlo; it<ithi; ++it)
      g[it] = f[it-lag];
    for (int it=ithi; it<nt; ++it)
      g[it] = 0.0f;
    return g;
  }
  /**
* Delays the CMP gather f by specified lag (which may be negative).
*/
  private static float[][] delay(int lag, float[][] f) {
    int nt = f[0].length;
    int nx = f.length;
    int itlo = max(0,lag); // 0 <= it-lag
    int ithi = min(nt,nt+lag); // it-lag < nt
    float[][] g = new float[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<itlo; ++it)
        g[ix][it] = 0.0f;
      for (int it=itlo; it<ithi; ++it)
        g[ix][it] = f[ix][it-lag];
      for (int it=ithi; it<nt; ++it)
        g[ix][it] = 0.0f;
    }
    return g;
  }

  private double dot(float[] f, float[] g) {
    int nt = f.length;
    double sum = 0.0;
    for (int it=_itmin; it<=_itmax; ++it)
      sum += f[it]*g[it];
    return sum;
  }
  
  private double dot(float[][] f, float[][] g) {
    int nt = f[0].length;
    int nx = f.length;
    double sum = 0.0;
    for (int ix=0; ix<nx; ++ix)
      for (int it=_itmin; it<=_itmax; ++it)
        sum += f[ix][it]*g[ix][it];
    return sum;
  }

  public static float[][] applyFilter(
    int nh, int kh, float[] h, float[][] f)
  {
    int nt = f[0].length;
    int nx = f.length;
    float[][] g = new float[nx][nt];
    applyFilter(nh,kh,h,f,g);
    return g;
  }

  private static void applyFilter(
    int nh, int kh, float[] h, float[][] f, float[][] g)
  {
    int nt = f[0].length;
    int nx = f.length;
    for (int ix=0; ix<nx; ++ix)
      conv(nh,kh,h,nt,0,f[ix],nt,0,g[ix]);
    preserveLeadingZeros(f,g);
  }

  private static void preserveLeadingZeros(float[][] f, float[][] g) {
    int nx = f.length;
    for (int ix=0; ix<nx; ++ix)
      preserveLeadingZeros(f[ix],g[ix]);
  }

  private static void preserveLeadingZeros(float[] f, float[] g) {
    int nt = f.length;
    int nz = countLeadingZeros(f);
    for (int it=0; it<nz; ++it)
      g[it] = 0.0f;
  }

  private static int countLeadingZeros(float[] f) {
    int n = f.length;
    int nz = 0;
    for (int i=0; i<n && f[i]==0.0f; ++i)
      ++nz;
    return nz;
  }


  

}
