/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package wwarp;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.mosaic.*;

/**
 * Warp correction.
 * For a gather of input seismic traces f(t,x), warp correction is defined
 * by
 * the transformation g(u,x) = a(u,x)*f(t(u,x),x), where t(u,x) are times at
 * which to evaluate f(t,x), a(u,x) are amplitude scale factors, and g(u,x) is
 * an warp-corrected gather of output traces.
 * <p>
 * The times t(u,x) may be specified directly, or computed from a 
 * dynamic warping algorithm
 * Amplitudes a(u,x) may also be specified,
 * or computed to preserve mutes and/or limit the maximum warp stretch
 * or squeeze. 
 * Note
 * that times t(u,x), amplitudes a(u,x), and the output gather g(u,x) are
 * functions of output time u, whereas f(t,x) is a function of input time t.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.12.17
 */
public class Warp {

  /**
   * Sets the maximum stretch or squeeze factor.
   * Used to zero amplitudes for which NMO stretch is excessive.
   * @param smax the maximum stretch or squeeze factor.
   */
  public void setStretchMax(double smax) {
    _smax = (float)smax;
  }

  /**
   * Applies this correction for specified shifts.
   * @param sft uniform time sampling of g.
   * @param sht uniform time sampling of shifts(h)
   * @param g array[nt] to be shifted.
   * @param h array[nt] of shifts.
   * @return array[nt] Shift-corrected output trace.
   */
   //CHECKEDCHECKEDCHECKEDCHECKEDCHECKEDCHECKEDCHECKED
  public float[] applyByShifts(Sampling sft, Sampling sgt,
    float[] g, float[] h) 
  {
    float[][] ta = timesAndAmplitudes(sft,sgt,g,h);
    float[] t = ta[0];
    float[] a = ta[1];
    return apply(sgt,t,a,g);
  }

  /**
   * Applies this correction for specified shifts.
   * @param sft uniform time sampling.
   * @param sht uniform time sampling.
   * @param sx uniform space sampling.
   * @param h array[nx][nt] of shifts.
   * @param f array[nx][nt] for input trace.
   * @return array[nx][nt] Shift-corrected output trace.
   */
   /*
  public float[][] applyByShifts(Sampling sft, Sampling sht, 
    Sampling sx, float[][] h, float[][] f) 
  {
    float[][][] ta = timesAndAmplitudes(sft,sht,sx,f,h);
    float[][] t = ta[0];
    float[][] a = ta[1];
    return apply(sft,sx,t,a,f);
  }
  */

  /**
   * Returns y(t) = x(u(t)).
   * Applies this correction for specified warped times.
   * @param u array[nt] of warped times.
   * @param f array[nt] for input trace.
   * @return array[nt] Shift-corrected output trace.
   */
  public float[] applyByTimes(float[] u, float[] f) {
    int nt = u.length;
    float[] y = new float[nt];
    _si.interpolate(nt,1.0,0.0,f,nt,u,y);
    y[0] *= u[1]-u[0];
    for (int it=1; it<nt; ++it)
      y[it] *= u[it]-u[it-1];
    return y;
  }

  /**
   * Returns y(t) = x(u(t)).
   * Applies this correction for specified warped times.
   * @param sx uniform space sampling.
   * @param u array[nt] of warped times.
   * @param f array[nt] for input trace.
   * @return array[nt] Shift-corrected output trace.
   */
  public float[][] applyByTimes(Sampling sx,
    float[][] u, float[][] f) 
  {
    int nx = sx.getCount();
    int nt = u[0].length;
    float[][] y = new float[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      _si.interpolate(nt,1.0,0.0,f[ix],nt,u[ix],y[ix]);
      y[ix][0] *= u[ix][1]-u[ix][0];
      for (int it=1; it<nt; ++it)
        y[ix][it] *= u[ix][it]-u[ix][it-1];
    }
    return y;
  }

  /**
   * Returns arrays of times and amplitudes for warp correction.
   * Sets to zero any amplitudes corresponding to (1) leading zeros in the
   * input gather or (2) samples for which NMO stretch is excessive. For all
   * other samples, non-zero amplitudes are simply the inverse of warp
   * stretch.
   * @param st uniform time sampling.
   * @param sx offset sampling; need not be uniform.
   * @param vnmo array[nt] of NMO velocities.
   * @param f array[nx][nt] for input gather f(t,x).
   * @return array {t,a} of times and amplitudes.
   */
  public float[][] timesAndAmplitudes(
    Sampling sft, Sampling sgt, float[] g, float[] h)
  {
    float[] t = getTimes(sft,sgt,sft,h);
    float[] a = getAmplitudes(sft,g,_smax,t);
    return new float[][]{t,a};
  }

  /**
   * Returns arrays of times and amplitudes for warp correction.
   * Sets to zero any amplitudes corresponding to (1) leading zeros in the
   * input gather or (2) samples for which NMO stretch is excessive. For all
   * other samples, non-zero amplitudes are simply the inverse of warp
   * stretch.
   * @param sft uniform time sampling.
   * @param sht uniform time sampling.
   * @param sx offset sampling; need not be uniform.
   * @param f array[nx][nt] for input gather f(t,x).
   * @param h array[nx][nt] for shifts.
   * @return array {t,a} of times and amplitudes.
   */
  public float[][][] timesAndAmplitudes(
    Sampling sft, Sampling sht, Sampling sx, float[][] f, float[][] h)
  {
    float[][] t = getTimes(sft,sht,sx,h);
    float[][] a = getAmplitudes(sft,sx,f,_smax,t);
    return new float[][][]{t,a};
  }

  /**
   * Applies this correction with specified times and amplitudes.
   * @param st uniform time sampling.
   * @param t array[nt] of times t(u).
   * @param a array[nt] of amplitudes a(u).
   * @param f array[nt] for input gather f(t).
   * @return array[nt] for output warped trace g(u).
   */
  public float[] apply(Sampling sft, float[] t, float[] a, float[] f) {
    int nt = sft.getCount();
    float[] g = new float[nt];
    float gi = 0.0f;

    for (int it=0; it<nt; ++it) {
      gi = _si.interpolate(sft,f,t[it]);
      g[it] = gi*a[it];
    }
    return g;
  }

  /**
   * Applies this correction for specified times and amplitudes.
   * @param sft uniform time sampling.
   * @param sx uniform time sampling.
   * @param t array[nt] of times t(u,x).
   * @param a array[nt] of amplitudes a(u,x).
   * @param f array[nt] for input gather f(t,x).
   * @return array[nt] for output NMO-corrected gather g(u,x).
   */
  public float[][] apply(Sampling sft, Sampling sx, 
    float[][] t, float[][] a, float[][] f) 
  {
    int nt = f.length;
    double dt = sft.getDelta();
    double ft = sft.getFirst();
    int nx = sx.getCount();
    float[][] g = new float[nx][nt];
    float gi = 0.0f;

    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        gi = _si.interpolate(sft,f[ix],t[ix][ix]);
        g[ix][it] = gi*a[ix][it];
      }
    }
    return g;
  }

  /**
   * For each offset, counts the number of leading zeros in a trace.
   * @param f array[nx][nt] in which to count leading zeros.
   * @return array[nx] counts of leading zeros, one for each trace.
   */
  public static int[] countLeadingZeros(float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    int[] ifnz = new int[nx];
    for (int ix=0; ix<nx; ++ix)
      ifnz[ix] = countLeadingZeros(f[ix]);
    return ifnz;
  }

  /**
   * For each time sample, counts the number of non-zero values.
   * Omits only leading zeros from the counts. In other words, for every
   * offset, any zero values that occur after the first non-zero value are
   * included in the returned counts. These counts are typically used to
   * normalize a stack over offsets, after NMO correction.
   * @param f array[nx][nt] in which to count non-zero values.
   * @return array[nt] of counts, one for each time sample.
   */
  public static int[] countNonZero(float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    int[] nnz = new int[nt];
    for (int ix=0; ix<nx; ++ix) {
      int nz = countLeadingZeros(f[ix]);
      for (int it=nz; it<nt; ++it)
        nnz[it] += 1.0f;
    }
    return nnz;
  }

  //////////////////////////////////////////////////////////////
  // private

  private SincInterp _si = SincInterp.fromErrorAndFrequency(0.01,0.40);
  private float _smax = 0.1f*Float.MAX_VALUE;

  private static int countLeadingZeros(float[] f) {
    int n = f.length;
    int nz = 0;
    for (int i=0; i<n && f[i]==0.0f; ++i)
      ++nz;
    return nz;
  }

  /**
   * Calculates the warped times from the existing 
   * sampling and the given shifts.
   * @param sft the sampling warping to.
   * @param sht the sampling warping from.
   * @param h shifts
   */
  private static float[] getTimes(Sampling sft, Sampling sgt, 
    Sampling sht, float[] h) {
    int nt = sht.getCount();
    float[] twarped = new float[nt];
    for (int it=0; it<nt; ++it) {
      float t = (float)sgt.getValue(it);
      twarped[it] = t+h[it]*(float)(sgt.getDelta());
    }
    return twarped;
  }

  /**
   * Calculates the warped times from the existing 
   * sampling and the given shifts.
   * @param sft the sampling warping to.
   * @param sht the sampling warping from.
   * @param h shifts
   */
  private static float[][] getTimes(Sampling sft, Sampling sht, 
    Sampling sx, float[][] h) 
  {
    int nt = sft.getCount();
    int nx = sx.getCount();
    float[][] twarped = new float[nx][nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        float t = (float)sft.getValue(it);
        twarped[ix][it] = t+h[ix][it]*(float)(sft.getDelta());
      }
    }
    return twarped;
  }

  private static float[] getTimes(Sampling st, float awarp) {
    int nt = st.getCount();
    float[] twarped = new float[nt];
    for (int it=0; it<nt; ++it) {
      float t = (float)st.getValue(it);
      twarped[it] = awarp*t;
    }
    return twarped;
  }

  private static float[] getAmplitudes(
    Sampling sgt, float[] g, float smax, float[] t) 
  {
    int nt = sgt.getCount();
    float dt = (float)sgt.getDelta();
    float ft = (float)sgt.getFirst();
    float odt = 1.0f/dt;
    float dtmin = dt/smax;
    float[] a = new float[nt];

    // Time of first non-zero input sample.
    int nz = countLeadingZeros(g);
    float tnz = ft+nz*dt;

    // Number of leading zeros in output. A leading output sample is zero 
    // if either (1) the corresponding input samples and all prior input 
    // samples are zero, or (2) stretch or squeeze would exceed the maximum.
    nz = 0;
    if (t[0]<tnz || t[1]-t[0]<dtmin)
      ++nz;
    for (int it=1; it<nt; ++it) {
      if (t[it]<tnz || t[it]-t[it-1]<dtmin)
        ++nz;
    }

    // Compute only the non-zero amplitudes. These amplitudes are simply the
    // inverse of stretch or squeeze.
    if (nz==0) {
      a[0] = (t[1]-t[0])*odt;
      ++nz;
    }
    for (int it=nz; it<nt; ++it)
      a[it] = (t[it]-t[it-1])*odt;

    return a;
  }
  private static float[][] getAmplitudes(
    Sampling sft, Sampling sx, float[][] f, float smax, float[][] t) 
  {
    int nt = sft.getCount();
    int nx = sx.getCount();
    float dt = (float)sft.getDelta();
    float ft = (float)sft.getFirst();
    float odt = 1.0f/dt;
    float dtmin = dt/smax;
    float[][] a = new float[nx][nt];

    for (int ix=0; ix<nx; ++ix) {
      // Time of first non-zero input sample.
      int nz = countLeadingZeros(f[ix]);
      float tnz = ft+nz*dt;

      // Number of leading zeros in output. A leading output sample is zero 
      // if either (1) the corresponding input samples and all prior input 
      // samples are zero, or (2) stretch or squeeze would exceed the maximum.
      nz = 0;
      if (t[ix][0]<tnz || t[ix][1]-t[ix][0]<dtmin)
        ++nz;
      for (int it=1; it<nt; ++it) {
        if (t[ix][it]<tnz || t[ix][it]-t[ix][it-1]<dtmin)
          ++nz;
      }

      // Compute only the non-zero amplitudes. These amplitudes are simply the
      // inverse of stretch or squeeze.
      if (nz==0) {
        a[ix][0] = (t[ix][1]-t[ix][0])*odt;
        ++nz;
      }
      for (int it=nz; it<nt; ++it)
        a[ix][it] = (t[ix][it]-t[ix][it-1])*odt;
    }

    return a;
  }
  private static float[] getAmplitudes(
    Sampling st, float[] f, float smax, float awarp) 
  {
    int nt = st.getCount();
    float dt = (float)st.getDelta();
    float ft = (float)st.getFirst();
    float odt = 1.0f/dt;
    float dtmin = dt/smax;
    float[] a = new float[nt];

    for (int it=0; it<nt; ++it)
      a[it] = awarp;

    return a;
  }

}

