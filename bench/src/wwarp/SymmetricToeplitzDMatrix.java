/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package wwarp;

import edu.mines.jtk.lapack.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A symmetric Toeplitz matrix is a square matrix specified by one row.
 * Elements of a Toeplitz matrix are Aij = a[i-j]. In other words, all
 * elements on each diagonal of a Toeplitz matrix are equal. This class
 * therefore enables solution of equations of the form
 * <pre><code>
 *  |a[0]    a[1]    a[2]    a[3]| |x[0]|     |b[0]|
 *  |a[1]    a[0]    a[1]    a[2]| |x[1]|  =  |b[1]|
 *  |a[2]    a[1]    a[0]    a[1]| |x[2]|     |b[2]|
 *  |a[3]    a[2]    a[1]    a[0]| |x[3]|     |b[3]|
 * </code></pre>
 * @author Dave Hale, Colorado School of Mines
 * @version 2013.08.25
 */
public class SymmetricToeplitzDMatrix {

  /**
   * Stores the top row of the symmetric Toeplitz matrix
   * @param a array of elements for the top row of the matrix.
   */
  public SymmetricToeplitzDMatrix(DMatrix a) {
    System.out.println(a.toString());
    int ncols = a.getN();
    _a = new double[ncols];
    for (int ic=0; ic<ncols; ++ic) {
      _a[ic] = a.get(0,ic);
    }
    dump(_a);
    _t = new double[ncols];
  }

  /**
   * Solves this symmetric Toeplitz system for specified right-hand-side.
   * @param b input array (Matrix) containing the right-hand-side column vector.
   * @return array containing the left-hand-side solution vector.
   */
  public double[] solve(DMatrix b) {
    int nrows = b.getM();
    double[] c = new double[nrows];
    for (int ir=0; ir<nrows; ++ir) {
      c[ir] = b.get(ir,0);
    }
    double[] x = new double[nrows];
    solve(c,x);
    return x;
  }

  /**
   * Solves this symmetric Toeplitz system for specified right-hand-side.
   * @param b input array containing the right-hand-side column vector.
   * @param x output array containing the left-hand-side solution vector.
   */
  public void solve(double[] b, double[] x) {
    solve(_t,_a,b,x);
  }

  /**
   * Solves a symmetric Toeplitz system for specified right-hand-side.
   * @param a input array of elements for top row of matrix.
   * @param b input array containing the right-hand-side column vector.
   * @param x output array containing the left-hand-side vector of unknowns.
   */
  public static void solve(double[] a, double[] b, double[] x) {
    double[] t = new double[a.length];
    solve(t,a,b,x);
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private double[] _a; // top row of matrix
  private double[] _t; // work array

  // Solution of Ax = b using a work array t and Levinson recursion.
  private static void solve(double[] t, double[] a, double[] b, double[] x) {
    int n = a.length;
    t[0] = 1.0;
    double v = a[0];
    x[0] = b[0]/a[0];
    for (int i=1; i<n; ++i) {
      t[i] = 0.0;
      x[i] = 0.0;
      double e = 0.0;
      for (int j=0; j<i; ++j)
        e += t[j]*a[i-j];
      double c = e/v;
      v -= c*e;
      for (int j=0; j<=i/2; ++j) {
        double timj = t[i-j]-c*t[j];
        t[j] -= c*t[i-j];
        t[i-j] = timj;
      }
      double w = 0.0;
      for (int j=0; j<i; ++j)
        w += x[j]*a[i-j];
      c = (w-b[i])/v;
      for (int j=0; j<=i; ++j)
        x[j] -= c*t[i-j];
    }
  }
}


