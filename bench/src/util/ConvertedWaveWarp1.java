package util;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.SincInterpolator;

public class ConvertedWaveWarp1 {

  public ConvertedWaveWarp1(double umax, int n, boolean uS) {
    _n = n;
    _umax = umax;
    _umin = 0.0;
    _scale = (uS)?0.5:3.0;
  }
  
  public float[] warp(float[] f) {
    SincInterpolator si = new SincInterpolator();
    si.setUniform(_n,1.0,0.0,f);
    float[] g = new float[_n];
    for (int i=0; i<_n; ++i) {
      double x = i-u(i);
//      System.out.println(x);
      g[i] = si.interpolate(x);
    }
    return g;
  }
  
  public float[] getShifts(float[] f) {
    float[] u = new float[_n];
    for (int i=0; i<_n; ++i) {
      double x = u(i);
      u[i] = (float)x;
//      u[i] = si.interpolate(x);
    }
    return u;
  }
  
  private int _n;
  private double _umax;
  private double _umin;
  private double _scale;
  
  private double u(double x) {
    if (x-_umin<0.0)
      return 0.0;
    double u = _scale*pow(x-_umin,0.5);
    return (u<_umax)?u:_umax;
  }
  
}
