package warp;

import static edu.mines.jtk.util.ArrayMath.*;

public class DynamicWarpingCIPs {

  public DynamicWarpingCIPs() {
    
  }
  
  public static float[][] computeErrors(float[] f, float[] g) {
    int n1 = f.length;
    int nl = n1;
    float[][] e = new float[n1][nl];
    for (int i1=0; i1<n1; n1++) {
      for (int il=0; il<nl; il++) {
        e[i1][il] = pow(abs(f[i1]-g[il]),2);
      }
    }
    return e;
  }
  
}
