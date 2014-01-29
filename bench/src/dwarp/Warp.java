/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dwarp;

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
    Sampling st, float[] shifts, float[] f)
  {
    float[] t = getTimes(st,shifts);
    float[] a = getAmplitudes(st,f,_smax,t);
    return new float[][]{t,a};
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
    Sampling st, float awarp, float[] f)
  {
    float[] t = getTimes(st,awarp);
    float[] a = getAmplitudes(st,f,_smax,awarp);
    return new float[][]{t,a};
  }

  /**
   * Applies this correction for specified times and amplitudes.
   * @param st uniform time sampling.
   * @param t array[nt] of times t(u,x).
   * @param a array[nt] of amplitudes a(u,x).
   * @param f array[nt] for input gather f(t,x).
   * @return array[nt] for output NMO-corrected gather g(u,x).
   */
  public float[] apply(Sampling st, float[] t, float[] a, float[] f) {
    int nt = f.length;
    double dt = st.getDelta();
    double ft = st.getFirst();
    float[] g = new float[nt];
    float gi = 0.0f;

    System.out.println("hi");
    for (int it=0; it<nt; ++it) {
      gi = _si.interpolate(st,f,t[it]);
      g[it] = gi*a[it];
    }
    return g;
  }


  /**
   * Applies this correction for specified shifts.
   * @param st uniform time sampling.
   * @param shifts array[nt] of shifts.
   * @param f array[nt] for input trace.
   * @return array[nt] Shift-corrected output trace.
   */
  public float[] apply(Sampling st, float[] shifts, float[] f) {
    float[][] ta = timesAndAmplitudes(st,shifts,f);
    float[] t = ta[0];
    float[] a = ta[1];
    return apply(st,t,a,f);
  }

  /**
   * Applies this correction for specified shifts.
   * @param st uniform time sampling.
   * @param shifts array[nt] of shifts.
   * @param f array[nt] for input trace.
   * @return array[nt] Shift-corrected output trace.
   */
  public float[] apply(Sampling st, float awarp, float[] f) {
    float[][] ta = timesAndAmplitudes(st,awarp,f);
    
    float[] t = ta[0];
    float[] a = ta[1];
    return apply(st,t,a,f);
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

  private SincInterp _si = SincInterp.fromErrorAndFrequency(0.01,0.45);
  private float _smax = 0.1f*Float.MAX_VALUE;

  private static int countLeadingZeros(float[] f) {
    int n = f.length;
    int nz = 0;
    for (int i=0; i<n && f[i]==0.0f; ++i)
      ++nz;
    return nz;
  }

  private static float[] getTimes(Sampling st, float[] shifts) {
    int nt = st.getCount();
    float[] twarped = new float[nt];
    for (int it=0; it<nt; ++it) {
      float t = (float)st.getValue(it);
      twarped[it] = t+shifts[it];
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
    Sampling st, float[] f, float smax, float[] t) 
  {
    int nt = st.getCount();
    float dt = (float)st.getDelta();
    float ft = (float)st.getFirst();
    float odt = 1.0f/dt;
    float dtmin = dt/smax;
    float[] a = new float[nt];

    // Time of first non-zero input sample.
    int nz = countLeadingZeros(f);
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

