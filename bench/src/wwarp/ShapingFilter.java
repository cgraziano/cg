/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package wwarp;

import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A shaping filter h, such that convolution h * x ~ y for specified x and y.
 * The filter h shapes (approximately) a specified input sequence x into a
 * specified output sequence y. The approximation to minimize the sum of
 * squared differences between the sequences h*x and y.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.25
 */
public class ShapingFilter {

  public static float[] design(
    int itmin, int itmax,
    int nh, int kh, 
    int nx, int kx, float[] x, 
    int ny, int ky, float[] y) 
  {
    float[] cxx = new float[nh];
    float[] cxy = new float[nh];
    int mt = 1+itmax-itmin;
    x = copy(mt,itmin,x);
    y = copy(mt,itmin,y);
    xcor(mt,kx,x,mt,kx,x,nh, 0,cxx);
    //cxx[0] *= 1.0001;
    xcor(mt,kx,x,mt,ky,y,nh,kh,cxy);
    float cxxMax = max(cxx);
    float cxyMax = max(cxy);
    float max = max(cxxMax,cxyMax);

    cxx = div(cxx,max);
    cxy = div(cxy,max);
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cxx);
    return stm.solve(cxy);
  }

  public static float[] design(
    int itmin, int itmax,
    int nh, int kh, 
    int nx1, int kx, float[][] x, 
    int ny1, int ky, float[][] y) 
  {
    int nx2 = x.length;
    float[] cxx = new float[nh];
    float[] cxy = new float[nh];
    float[] cxxSum = new float[nh];
    float[] cxySum = new float[nh];
    int mt = 1+itmax-itmin;
    x = copy(mt,nx2,itmin,0,x);
    y = copy(mt,nx2,itmin,0,y);
    for (int ix=0; ix<nx2; ++ix) {
      xcor(mt,kx,x[ix],mt,kx,x[ix],nh, 0,cxx);
      //cxx[0] *= 1.0001;
      xcor(mt,kx,x[ix],mt,ky,y[ix],nh,kh,cxy);
      cxxSum = add(cxxSum,cxx);
      cxySum = add(cxySum,cxy);
    }
    float cxxMax = max(cxxSum);
    float cxyMax = max(cxySum);
    float max = max(cxxMax,cxyMax);

    cxxSum = div(cxxSum,max);
    cxySum = div(cxySum,max);
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cxxSum);
    return stm.solve(cxySum);
  }

  public static float[] design(
    int itmin, int itmax,
    int nh, int kh, 
    int nx1, int kx, float[][][] x, 
    int ny1, int ky, float[][][] y) 
  {
    int nx2 = x[0].length;
    int nx3 = x.length;
    float[] cxx = new float[nh];
    float[] cxy = new float[nh];
    float[] cxxSum = new float[nh];
    float[] cxySum = new float[nh];
    int mt = 1+itmax-itmin;
    x = copy(mt,nx2,nx3,itmin,0,0,x);
    y = copy(mt,nx2,nx3,itmin,0,0,y);
    for (int ix3=0; ix3<nx3; ++ix3) {
      for (int ix2=0; ix2<nx2; ++ix2) {
        xcor(mt,kx,x[ix3][ix2],mt,kx,x[ix3][ix2],nh, 0,cxx);
        //cxx[0] *= 1.0001;
        xcor(mt,kx,x[ix3][ix2],mt,ky,y[ix3][ix2],nh,kh,cxy);
        cxxSum = add(cxxSum,cxx);
        cxySum = add(cxySum,cxy);
      }
    }
    float cxxMax = max(cxxSum);
    float cxyMax = max(cxySum);
    float max = max(cxxMax,cxyMax);

    cxxSum = div(cxxSum,max);
    cxySum = div(cxySum,max);
    SymmetricToeplitzFMatrix stm = new SymmetricToeplitzFMatrix(cxxSum);
    return stm.solve(cxySum);
  }

}
