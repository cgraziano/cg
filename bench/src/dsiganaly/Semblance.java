/****************************************************************************
Copyright (c) 2013, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package dsiganaly;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.*;


/**
 * Calculates the semblance of the 2D gather.
 * Two methods have been proposed, one where the traces of the gather are
 * averaged together and one where the traces of the gather are averaged 
 * together, but does not allow zeros to be included in this average.
 */
public class Semblance {

  /**
   * Calculates the the semblance of the given 2D gather.
   * The averages used in this semblance will include all of the values
   * in the gather, including zero values.
   * @param f gather
   */
  public static float semblance(float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    float[] fat = new float[nx];
    float[] fax = new float[nt];
    float[] fax2 = new float[nt];
    float[][] f2 = new float[nx][nt];
    float num = 0;
    float den = 0;

    //numerator <<g>x^2>t <> denotes average
    for (int x=0; x<nx; ++x)
      fax = add(fax,f[x]);
    fax = div(fax,nx);
    fax2 = mul(fax,fax);
    for (int t=0; t<nt; ++t)
      num += fax2[t];
    num /= nt;

    //zero fat and fax arrays
    fat = new float[nx];
    fax = new float[nt];

    //<<g^2>x>t denominator 
    f2 = mul(f,f);
    for (int x=0; x<nx; ++x)
      fax = add(fax,f2[x]);
    fax = div(fax,nx);
    for (int t=0; t<nt; ++t)
      den += fax[t];
    den /= nt;

    return num/den;
  }

  /**
   * Calculates the semblance of the given 2D gather. 
   * The averages used in this semblance will include all of the values
   * in the gather, EXCEPT groups of zeros that are next to each other
   * at a particular time in the gather. This type of pattern would occur
   * when a mute has been applied to the gather. To ensure the semblance
   * value is not corrupted with zeros in the averages, the zeros are not
   * included.
   * @param f gather
   */
  public static float smartSemblance(float[][] f) {
    int nx = f.length;
    int nt = f[0].length;
    float[] fat = new float[nx];
    float[] fax = new float[nt];
    float[] fax2 = new float[nt];
    float[][] f2 = new float[nx][nt];
    float num = 0;
    float den = 0;
    float cutoff = .000001f;

    //numerator <<g>x^2>t <> denotes average
    float[][] hit = new float[nx][nt];
    int[] countnx = new int[nt];
    for (int t=0; t<nt; ++t) {
      for (int x=0; x<nx; ++x) {
        if (-cutoff > f[x][t] || f[x][t] > cutoff) {
          hit[x][t] = 1.0f;
          fax[t] = fax[t] + f[x][t];
          countnx[t] += 1;
        }
      }
    }
    System.out.println("countnx[0] = "+countnx[0]);
    SimplePlot sp = new SimplePlot();
    sp.addPixels(hit);
    sp.addColorBar();

    for (int t=0; t<nt; ++t) {
      if (countnx[t] == 0) {
        fax[t] = 0;
      }
      else
        fax[t] = fax[t]/countnx[t];
    }
    
      
    fax2 = mul(fax,fax);
    for (int t=0; t<nt; ++t)
      num += fax2[t];
    num /= nt;

    //zero fat and fax arrays
    fat = new float[nx];
    fax = new float[nt];

    //<<g^2>x>t denominator 
    f2 = mul(f,f);
    System.out.println("f2 = "+f2[0][0]);
    countnx = new int[nt];
    for (int t=0; t<nt; ++t) {
      for (int x=0; x<nx; ++x) {
        if (-cutoff*cutoff > f2[x][t] || f2[x][t] > cutoff*cutoff) {
          fax[t] = fax[t] + f2[x][t];
          countnx[t] += 1;
        }
      }
      
    }
    System.out.println("countnx[0] = "+countnx[0]);

    for (int t=0; t<nt; ++t) {
      if (countnx[t] == 0) {
        fax[t] = 0;
        System.out.println("hi");
      }
        
      else
        fax[t] = fax[t]/countnx[t];
    }
    System.out.println("fax[0] = "+fax[0]);
    
    for (int t=0; t<nt; ++t) {
      den += fax[t];
    }
    den /= nt;
    System.out.println("num = "+num);
    System.out.println("den = "+den);

    return num/den;

  }
}
