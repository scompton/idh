package util;

import static edu.mines.jtk.util.ArrayMath.*;

public class FakeDataSimple {

  public static float[] makeSimpleTrace(
    int n1, double fpeak, int rickerIndex) 
  {
    float[] trc = zerofloat(n1);
    int ih = (int) (sqrt(3.0/2.0) / (PI*fpeak));
    ih *= 2;
    System.out.println(ih);
    int nh = 1+2*ih;
    float[] h = new float[nh];
    for (int jh=0;jh<nh;jh++)
      h[jh] = ricker(fpeak, jh-ih);
    copy(nh, 0, h, rickerIndex-ih, trc);
    return trc;
  }
  
  /**
   * Returns the amplitude A of a ricker wavelet with a given
   * peak frequency and a given time 
   * @param fpeak
   * @param time
   * @return the amplitude of the ricker wavelet
   */
  private static float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }
}
