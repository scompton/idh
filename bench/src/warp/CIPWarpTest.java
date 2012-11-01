package warp;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.DynamicWarping;

public class CIPWarpTest {

  public static float[][] test1() {
    int N = 100;
    float[] f = zerofloat(N);
    float[] w = zerofloat(N);
    for (int n=0; n<N; n++) {
      if (n==50)
        f[n] = 100;
    }
    return new float[][]{f,w};
  }
  
  public static float[][] findShifts(
      float[][] f, float[][]g, DynamicWarping dw) {
    int n1 = f[0].length;
    int n2 = f[1].length;
    float[][] u = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      u[i2] = dw.findShifts(f[i2],g[i2]);
    }
    return u;
  }
  
  public static float[][] applyShifts(
      float[][] u, float[][] g, DynamicWarping dw) {
    int n1 = g[0].length;
    int n2 = g[1].length;
    float[][] h = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      dw.applyShifts(u[i2],g[i2],h[i2]);
    }
    return h;
  }
  
}
