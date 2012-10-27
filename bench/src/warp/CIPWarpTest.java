package warp;

import static edu.mines.jtk.util.ArrayMath.*;

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
  
}
