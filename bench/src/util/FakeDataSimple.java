package util;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.Random;

import edu.mines.jtk.dsp.Conv;
import edu.mines.jtk.mosaic.SimplePlot;

public class FakeDataSimple {

  public static float[] makeSimpleTrace(
    int n1, double fpeak, int rickerIndex) 
  {
    float[] trc = zerofloat(n1);
    int ih = (int) (sqrt(3.0/2.0) / (PI*fpeak));
    ih *= 2;
    int nh = 1+2*ih;
    float[] h = new float[nh];
    for (int jh=0;jh<nh;jh++)
      h[jh] = ricker(fpeak, jh-ih);
    copy(nh, 0, h, rickerIndex-ih, trc);
    return trc;
  }
  
  public static float[] makeRandomTrace(int n1, double fpeak, long seed) {
    Random r = new Random(seed);
    float[] f = pow(mul(2.0f, sub(randfloat(r, n1), 0.5f)), 15.0f);
//    SimplePlot.asPoints(f);
    int ih = (int)(3.0/fpeak);
    int nh = 1+ih*2;
    float[] h = new float[nh];
    for (int jh=0; jh<nh; jh++)
      h[jh] = ricker(fpeak, jh-ih);
//    SimplePlot.asPoints(f);
    float[] g = new float[n1];
    Conv.conv(nh, -ih, h, n1, 0, f, n1, 0, g);
//    SimplePlot.asPoints(g);
    return g;
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
  
  public static void main(String[] args) {
    makeRandomTrace(1001, 0.05, 1L);
  }
}
