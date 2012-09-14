package stke;

import edu.mines.jtk.dsp.RecursiveGaussianFilter;
import static edu.mines.jtk.util.ArrayMath.*;

public class Correlate {
  
  private RecursiveGaussianFilter _rgf1; 
  private RecursiveGaussianFilter _rgf2; 
  private RecursiveGaussianFilter _rgf3; 
  
  public Correlate(double sigma1, double sigma2, double sigma3) {
    _rgf1 = new RecursiveGaussianFilter(sigma1);
    _rgf2 = new RecursiveGaussianFilter(sigma2);
    _rgf3 = new RecursiveGaussianFilter(sigma3);
  }
  
  public float[][][] compute(float[][][] g, float[][][] s) {
    int n1 = g[0][0].length;
    int n2 = g[0].length;
    int n3 = g.length;
    float[][][] gg = new float[n3][n2][n1];
    float[][][] ss = new float[n3][n2][n1];
    float[][][] gs = new float[n3][n2][n1];
    float[][][] c  = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          gg[i3][i2][i1] = g[i3][i2][i1]*g[i3][i2][i1];
          ss[i3][i2][i1] = s[i3][i2][i1]*s[i3][i2][i1];
          gs[i3][i2][i1] = g[i3][i2][i1]*s[i3][i2][i1];
        }
      }
    }
    _rgf1.apply0XX(gg,gg);
    _rgf2.applyX0X(gg,gg);
    _rgf3.applyXX0(gg,gg);
    _rgf1.apply0XX(ss,ss);
    _rgf2.applyX0X(ss,ss);
    _rgf3.applyXX0(ss,ss);
    _rgf1.apply0XX(gs,gs);
    _rgf2.applyX0X(gs,gs);
    _rgf3.applyXX0(gs,gs);
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          c[i3][i2][i1] = gs[i3][i2][i1]/sqrt(gg[i3][i2][i1]*ss[i3][i2][i1]);
        }
      }
    }
    return c;
  }
}