package util;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.BilinearInterpolator2;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.util.Check;

public class L1C2Interpolator {

  public L1C2Interpolator(float[][] x1, float[] x2, float[][] y) {
    int n2 = x2.length;
    Check.argument(n2==x1.length,"x1 array specified for all n2");
    for (int i2=0; i2<n2; i2++)
      Check.argument(isMonotonic(x1[i2]), "array x1 is monotonic");
    Check.argument(isMonotonic(x2), "array x2 is monotonic");
    _x1 = copy(x1);
    _x2 = copy(x2);
    _n1 = _x1[0].length;
    _n2 = _x2.length;
    makeCoefficients(_x1,_x2,y);
  }
  
  public float[][] interpolate(Sampling s1, Sampling s2) {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    float[][] y = new float[n2][n1];
    int n2x1 = _x2.length;
    int[][] k1a = new int[n2x1][];
    for (int i2=0; i2<n2x1; i2++)
      k1a[i2] = makeIndices(s1,_x1[i2]);
    int[] k2a = makeIndices(s2,_x2);
    for (int i2=0; i2<n2; ++i2) {
      float x2 = (float)s2.getValue(i2);
      int k21 = k2a[i2];
      int k20 = k21-1;
      int k22 = k21+1;
      int k23 = k21+2;
      int[] k2 = new int[]{k20,k21,k22,k23};
      for (int i1=0; i1<n1; ++i1) {
        float x1 = (float)s1.getValue(i1);
        int k11 = k1a[k2[1]][i1];
        int k10 = (k2[0]> -1)?k1a[k2[0]][i1]:k11;
        int k12 = (k2[2]<_n2)?k1a[k2[2]][i1]:k11;
        int k13 = (k2[3]<_n2)?k1a[k2[3]][i1]:k12;
        int[] k1 = new int[]{k10,k11,k12,k13};
        y[i2][i1] = interpolate(x1,x2,k1,k2);
      }
    }
    return y;
  }

  private int _n1,_n2;
  private float[][] _x1;
  private float[] _x2;
  private float[][] _a00,_a10; // interpolation coefficients
  
  private static int index(float x, float[] xs, int i) {
    i = binarySearch(xs,x,i);
    if (i<0) 
      i = (i<-1)?-2-i:0;
    if (i>=xs.length-1)
      i = xs.length-2;
    return i;
  }

  private float interpolate(float x1, float x2, int[] k1, int[] k2) {
    float[] y = new float[4];
    float[] x = new float[4];
    for (int ik=0; ik<4; ik++) {
      int ik1 = k1[ik];
      int ik2 = k2[ik];
      if (ik2<0) {
        x[ik] = _x2[0]-(_x2[1]-_x2[0]);
        ik2 = k2[ik+1];
        ik1 = k1[ik+1];
      } else if (ik2==_n2-1) {
        x[ik] = _x2[_n2-1];
        ik2 = k2[ik-1];
        ik1 = k1[ik-1];
      } else if (ik2>_n2-1) {
        x[ik] = _x2[_n2-1]+(_x2[_n2-1]-_x2[_n2-2]);
        ik2 = k2[ik-2];
        ik1 = k1[ik-2];
      } else
        x[ik] = _x2[ik2];
      float d1 = x1-_x1[ik2][ik1];
      assert(d1>=0):"x1="+x1+", x2="+x2+", ik2="+ik2+", ik1="+ik1+", ik="+ik;
      y[ik] = _a00[ik2][ik1]+d1*_a10[ik2][ik1];
    }
    CubicInterpolator ci = new CubicInterpolator(x,y);
    return ci.interpolate(x2);
  }
  
  private static int[] makeIndices(Sampling si, float[] xs) {
    int n = si.getCount();
    int[] ki = new int[n];
    ki[0] = index((float)si.getValue(0),xs,0);
    for (int i=1; i<n; ++i)
      ki[i] = index((float)si.getValue(i),xs,ki[i-1]);
    return ki;
  }

  private void makeCoefficients(float[][] x1, float[] x2, float[][] y) {
    int n1 = x1[0].length;
    int n2 = x2.length;
    _a00 = new float[n2-1][n1-1];
    _a10 = new float[n2-1][n1-1];
    for (int k2=0; k2<n2-1; ++k2) {
      for (int k1=0; k1<n1-1; ++k1) {
        float d1 = x1[k2][k1+1]-x1[k2][k1];
        float y00 = y[k2  ][k1  ];
        float y10 = y[k2  ][k1+1];
        _a00[k2][k1] = y00;
        _a10[k2][k1] = (y10-y00)/d1;
      }
    }
  }
  
  public static void main(String[] args) {
    float[] x2 = new float[]{0.0f,4.0f,9.0f,12.0f,18.0f,23.0f,30.0f};
    float[][] x1 = new float[][]{
        {0.0f,12.0f,20.0f,25.0f},
        {0.0f,10.0f,18.0f,25.0f},
        {0.0f, 8.0f,19.0f,25.0f},
        {0.0f,10.0f,21.0f,25.0f},
        {0.0f,11.0f,20.0f,25.0f},
        {0.0f,13.0f,17.0f,25.0f},
        {0.0f,15.0f,19.0f,25.0f}};
    int n2 = 7;
    int n1 = 4;
    float[][] y = new float[n2][n1];
    float[][] x2copy = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        y[i2][i1] = x2[i2]+x1[i2][i1];
        x2copy[i2][i1] = x2[i2];
      }
    }
    SimplePlot.asPixels(y);
    Sampling s1 = new Sampling(26,1.0,0.0);
    Sampling s2 = new Sampling(31,1.0,0.0);
    L1C2Interpolator lci = new L1C2Interpolator(x1,x2,y);
    float[][] ylci = lci.interpolate(s1,s2);
    Viewer vlci = new Viewer(ylci);
    vlci.setTitle("Linear/BilinearInterpolator");
    vlci.addPoints2(x1,x2copy);
    vlci.setColorModel1(ColorMap.JET);
    vlci.setSize(900,900);
    vlci.show();
    
    BilinearInterpolator2 bli = new BilinearInterpolator2(x1,x2,y);
    float[][] ybli = bli.interpolate(s1,s2);
    Viewer vbli = new Viewer(ybli);
    vbli.setTitle("BilinearInterpolator");
    vbli.addPoints2(x1,x2copy);
    vbli.setColorModel1(ColorMap.JET);
    vbli.setSize(900,900);
    vbli.show();
  }
}
