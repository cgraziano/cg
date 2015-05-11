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
import edu.mines.jtk.opt.BrentMinFinder;
import java.awt.Color;
import static edu.mines.jtk.dsp.Conv.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Code to run: bcg && jr wwarp/TestAlpha1D
 *
 * @author Chris Graziano, Colorado School of Mines
 * @version 2015.03.29
 */
 public class TestAlpha1D {
   public static void main(String[] arg) {
     float[] coefficientsf = {10.0f,10.0f,0.0f,-4.0f,1.0f};
     float[] coefficientsfprime = {10.0f,0.0f,-12.0f,4.0f};
     TestAlpha1D ta = new TestAlpha1D(coefficientsf,coefficientsfprime);
     ta.getMinXGN(-10.0f,50);
   }

   /**
    * Constructs a TestAlpha test.
    * @param f the function you want to find the minimum of.
    * @param fprime the 1st derivative of f.
    */
   public TestAlpha1D(float[] coefficientsf, float[] coefficientsfprime) {
     Function f = new Function(coefficientsf);
     Function fprime = new Function(coefficientsfprime);
     _f = f;
     _fprime = fprime;
   }

   /**
    * Finds the minimum of a function using the Gauss-Newton method.
    * @param x0 The first guess for the Gauss-Newton method.
    * @param n The number of iterations the Gauss-Newton method should do.
    */
   public float getMinXGN(float x0, int n) {
     float x = x0;
     float delta = 0.0f;
     float steplength = 0.0f;
     float f = 0.0f;
     for (int i=0; i<n; ++i) {
       trace("starting x = "+x);
       delta = getDelta(x);
       trace("delta = "+delta);
       ErrorFunction ef = new ErrorFunction();
       ef.setParameters(x,delta,_f);
       BrentMinFinder bf = new BrentMinFinder(ef);
       steplength = (float) bf.findMin(0.0,1.0,0.0);
       x = x + delta*steplength;
       f = _f.getValue(x);
       trace("ff = "+(f*f));
       trace("steplength = "+steplength);
       trace("updated x = "+x);
     }
     return x;
   }

//////////////////////////////////////////////////////////////////////////////
//////////////////Private///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
  private float[] _coefficients;
  private Function _f,_fprime;

  private float getDelta(float x) {
    return -_f.getValue(x)/_fprime.getValue(x);
  }

  private static void trace(String s) {
    System.out.println(s);
  }

  private class ErrorFunction implements BrentMinFinder.Function {
    public void setParameters(float x, float delta, Function f) {
      _x = x;
      _delta = delta;
      _f = f;
    }
    public double evaluate(double steplength) {
      float x = _x + (float) steplength*_delta;
      double f = (double) _f.getValue(x);
      return f*f;
    }
    private float _delta, _x;
    private Function _f;
  }

  private class Function {
    /**
     * Constructs a polynomial with set coefficients.
     * @param coefficients An array of coefficients starting with the coefficient 
                            corresponding to x^0 all the way to x^n. Every coefficient must
                            be set.
    */
    public Function(float[] coefficients) {
      _coefficients = coefficients;
      _ncoefficients = coefficients.length;
    }

    /**
     * Returns the polynomial's value at x. 
     */
    public float getValue(float x) {
      float sum = 0.0f;
      for (int i=0; i<_ncoefficients; ++i) {
        sum = sum + _coefficients[i]*pow(x,i);
      }
      return sum;

    }
    private float[] _coefficients;
    private int _ncoefficients;
  }


 }
