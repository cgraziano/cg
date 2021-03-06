/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dsiganaly;
import static edu.mines.jtk.util.ArrayMath.*;
import static edu.mines.jtk.dsp.Conv.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.SimplePlot;
import linalgebra.ToeplitzRecursion;
 


/**
 * For testing and understanding!
 * Constructs various objective functions for analysis.
 * To this date, the objective functions constructed are:
 *    The PEF Objective Function for a trace
 *    The PEF Objective Function for a gather
 *    The Semblance Objective Function for a gather
 * @author Chris Graziano, CWP
 * @version 2013.12.18
 */
public class ObjectiveFunction {

  /**
   * Given a trace and a range of a1 and a2 inverse wavelet coeffiecients
   * (prediction error operator coefficients) will compute
   * the objective function at different points of a1 and a2.
   * a0 is assumed to be 1.
   * @param x input trace
   * @param a1 first inverse coefficient
   * @param a2 second inverse coefficient
   */
  public static float[][] PEFTrace(float[] x, 
      Sampling sa1, Sampling sa2) {
    int nt = x.length;
    int na1 = sa1.getCount();
    float da1 = (float) (sa1.getDelta());
    float fa1 = (float) (sa1.getFirst());
    int na2 = sa2.getCount();
    float da2 = (float) (sa2.getDelta());
    float fa2 = (float) (sa2.getFirst());

    float[][] obf = new float[na2][na1];

    float a1 = 0.0f;
    float a2 = 0.0f;
    int na = 3;
    float[] a = new float[na];
    a[0] = 1.0f;
    float sumE = 0.0f;
    float[] y = new float[nt];
    float yt = 0.0f;

    for (int ia1=0; ia1<na1; ++ia1) { 
      a1 = fa1+da1*ia1;
      a[1] = a1;
      for (int ia2=0; ia2<na2; ++ia2) {
        a2 = fa2+da2*ia2;
        a[2] = a2;
        conv(na,0,a,nt,0,x,nt,0,y);
        for (int t=0; t<nt; ++t) {
          yt = y[t]; 
          sumE += yt*yt;
        }
        obf[ia2][ia1] = sumE;
        sumE = 0.0f;
      }
    }
    
    return obf;
  } 

  /**
   * Given a CMP gather and a range of a1 and a2 inverse 
   * wavelet coeffiecients (prediction error operator coefficients), 
   * will compute the objective function at different points of a1 and a2.
   * a0 is assumed to be 1.
   * @param f input gather 
   * @param a1 first inverse coefficient
   * @param a2 second inverse coefficient
   */
  public static float[][] PEFGather(float[][] f, 
      Sampling sa1, Sampling sa2, Sampling sx, Sampling st) {
    int nx = f.length;
    int nt = f[0].length;
    int na1 = sa1.getCount();
    float da1 = (float) (sa1.getDelta());
    float fa1 = (float) (sa1.getFirst());
    int na2 = sa2.getCount();
    float da2 = (float) (sa2.getDelta());
    float fa2 = (float) (sa2.getFirst());

    float[][] obf = new float[na2][na1];
    float[][] fa = new float[nx][nt];//Convolution of a's with gather
    float[] fax = new float[nt];//fa summed across the offsets
    float[] fat = new float[nx];//fa summed across the times 

    float a1 = 0.0f;
    float a2 = 0.0f;
    int na = 3;
    float[] a = new float[na];
    a[0] = 1.0f;
    float sumE = 0.0f;
    float[] y = new float[nt];
    float yt = 0.0f;
    float[][] nmofa = new float[nx][nt];
    for (int ia1=0; ia1<na1; ++ia1) { 
      a1 = fa1+da1*ia1;
      a[1] = a1;
      for (int ia2=0; ia2<na2; ++ia2) {
        a2 = fa2+da2*ia2;
        a[2] = a2;
        for (int x=0; x<nx; ++x) {
          conv(na,0,a,nt,0,f[x],nt,0,fa[x]);
        }

        //fa = wn.applyNmo(st,sx,vnmos,fa);

        fat = new float[nx];
        fax = new float[nt];
        for (int x=0; x<nx; ++x)
          add(fax,mul(fa[x],fa[x]),fax);
        fax = div(fax,nx);
        for (int t=0; t<nt; ++t)
          sumE += fax[t];
        sumE /= nt;
        obf[ia2][ia1] = sumE;
        sumE = 0.0f;
      }
    }
    
    return obf;
  } 

  /**
   * Given a CMP gather and a range of a1 and a2 inverse 
   * wavelet coeffiecients (prediction error operator coefficients), 
   * will compute the (1-semblance) objective function at 
   * different points of a1 and a2.
   * a0 is assumed to be 1.
   * Minimizing (1-(<<g>x^2>t/<<g^2>x>t)),sum over x = <g>x
   * @param f input gather 
   * @param a1 first inverse coefficient
   * @param a2 second inverse coefficient
   */
  public static float[][] SemblanceGather(float[][] f, 
      Sampling sa1, Sampling sa2, Sampling sx, Sampling st) {
    int nx = f.length;
    int nt = f[0].length;
    int na1 = sa1.getCount();
    float da1 = (float) (sa1.getDelta());
    float fa1 = (float) (sa1.getFirst());
    int na2 = sa2.getCount();
    float da2 = (float) (sa2.getDelta());
    float fa2 = (float) (sa2.getFirst());

    float[][] obf = new float[na2][na1];
    float[][] fa = new float[nx][nt];//Convolution of a's with gather
    float[] fax = new float[nt];//fa summed across the offsets
    float[] fat = new float[nx];//fa summed across the times 

    float a1 = 0.0f;
    float a2 = 0.0f;
    int na = 3;
    float[] a = new float[na];
    a[0] = 1.0f;
    float num = 0.0f;
    float den = 0.0f;
    float sumE = 0.0f;
    float[] y = new float[nt];
    float yt = 0.0f;
    float odnt = 1.0f/nt;
    for (int ia1=0; ia1<na1; ++ia1) { 
      a1 = fa1+da1*ia1;
      a[1] = a1;
      for (int ia2=0; ia2<na2; ++ia2) {
        a2 = fa2+da2*ia2;
        a[2] = a2;
        for (int x=0; x<nx; ++x) {
          conv(na,0,a,nt,0,f[x],nt,0,fa[x]);
        }
        //fa = wn.applyNmo(st,sx,vnmos,fa);

        obf[ia2][ia1] = 1.0f-Semblance.smartSemblance(fa);
      }
    }
    
    return obf;
  } 


  
}


