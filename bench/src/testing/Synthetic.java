/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package testing;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;




/**
 * Different sythetics are contained within this class
 * to be used when testing.
 * @author Chris Graziano, CWP
 * @version 2013.12.18
 */ 
public class Synthetic {

  public static float[][] SyntheticCMPGather() {
    System.out.println("NotFinished");
    return null;
  }

  public static float[][] makeCMPReflections(float vel, int nref, 
      Sampling st, Sampling sx) {
    int nt = st.getCount();
    int nx = sx.getCount();
    float dt = (float)st.getDelta();
    float dx = (float)sx.getDelta();
    float ft = (float)st.getFirst();
    float fx = (float)sx.getFirst();
    float[][] p = new float[nt][nx];
    float[] rt = add(ft,mul((nt-1)*dt,randfloat(nref)));//reflection times
    float[] ra = sub(mul(2.0f,randfloat(nref)),1.0f);//reflection amplitudes
    SincInterp si = SincInterp.fromErrorAndFrequency(0.01, 0.45);
    float x = 0;
    float x2dv2 = 0;
    for (int ix=0; ix<nx; ++ix) {
      x = (float)sx.getValue(ix);//value of offset at the ix index
      x2dv2 = x*x/(vel*vel);
      for (int ir=0; ir<nref; ++ir) {
        float rti = rt[ir];
        float rai = ra[ir];
        rti = sqrt(rti*rti+x2dv2);
        si.accumulate(rti,rai,nt,dt,ft,p[ix]);
      }
    }
    return p;
  }

  /**
   * Creates trace with reflection coefficients
   */
  public static float[] makeTrace(Sampling st, float mt, int nref) {
    int nt = st.getCount();
    float dt = (float)st.getDelta();
    float ft = (float)st.getFirst();
    float[] p = new float[nt];
    float[] rt = add(ft+(mt-ft),mul((nt-1)*dt,randfloat(nref)));//reflection times
    float[] ra = sub(mul(2.0f,randfloat(nref)),1.0f);//reflection amplitudes
    float[] trace = new float[nt];
    SincInterp si = SincInterp.fromErrorAndFrequency(0.01, 0.45);
    for (int ir=0; ir<nref; ++ir) {
      float rti = rt[ir];
      float rai = ra[ir];
      si.accumulate(rti,rai,nt,dt,ft,p);
    }
    return p;
  }

  /**
   * Creates trace with reflection coefficients
   */
  public static float[] makeTraceTRef(Sampling st, float mt, float[] tref,
    float[] aref) {
    int nref = tref.length;
    int nt = st.getCount();
    float dt = (float)st.getDelta();
    float ft = (float)st.getFirst();
    float[] p = new float[nt];
    float[] rt = tref;
    float[] ra = aref;
    float[] trace = new float[nt];
    SincInterp si = SincInterp.fromErrorAndFrequency(0.01, 0.45);
    for (int ir=0; ir<nref; ++ir) {
      float rti = rt[ir];
      float rai = ra[ir];
      si.accumulate(rti,rai,nt,dt,ft,p);
    }
    return p;
  }

}
