package warp;

import static edu.mines.jtk.util.ArrayMath.*;

public class SlopeMax {

  public static double[] cosineTaper(int n1, double rMin, double rMax) {
    float twoPi = 2.0f*FLT_PI;
    float f = 1.0f/(4.0f*n1);
    double[] r = new double[n1]; 
    for (int i=0; i<n1; i++)
      r[i] = cos(twoPi*f*i);
    normalize(r,rMin,rMax);
    for (int i=0; i<n1; i++)
      assert r[i]<=rMax && r[i]>=rMin;
    return r;
  }
  
  private static void normalize(double[] f, double nMin, double nMax) {
    final int n1 = f.length;
    double fMin = min(f);
    double fMax = max(f);
    double range = fMax-fMin;
    double nrange = nMax-nMin;
    for (int i1=0; i1<n1; ++i1) {
      double value = f[i1];
      f[i1] = nrange*(value-fMin)/range + nMin;
    }
  }
  
}
