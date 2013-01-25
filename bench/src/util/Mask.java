package util;

import edu.mines.jtk.util.Check;

public class Mask {

  public static float[][] getMask(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] m = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        if (f[i2][i1]==0.0f)
          m[i2][i1] = 0.0f;
        else
          m[i2][i1] = 1.0f;
      }
    }
    return m;
  }
  
  public static void applyMask(float[][] f, float[][] m) {
    int n2 = f.length;
    int n1 = f[0].length;
    Check.argument(n2==m.length,"f.length==m.length");
    Check.argument(n1==m[0].length,"f[0].length==m[0].length");
    for (int i2=0; i2<n2; i2++)
      for (int i1=0; i1<n1; i1++)
        f[i2][i1] = f[i2][i1]*m[i2][i1];
  }
  
}
