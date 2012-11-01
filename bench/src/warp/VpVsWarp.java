package warp;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.util.Check;

public class VpVsWarp {

  public VpVsWarp(float a, float b, Sampling s) {
    Check.argument(a>=-1.0 & a<=1.0,"a>=-1.0 & a<=1.0");
    _a = a;
    _b = b;
    _s = s;
  }
  
  public float[] getVpVs() {
    int n1 = _s.getCount();
    double ft = _s.getFirst();
    double dt = _s.getDelta();
    float[] vpvs = new float[n1];
    for (int i1=0; i1<n1; i1++) {
      double t = ft+i1*dt;
      vpvs[i1] = 2.0f + _a*(float)cos(2.0*PI*_b*t);
    }
    return vpvs;
  }
  
  public float getAverage(float[] vpvs) {
    int n1 = vpvs.length;
    return sum(vpvs)/n1;
  }
  
  public float[] warp(float[] f, float[] vpvs) {
    int n1 = f.length;
    float[] u = getU(vpvs);
    SincInterpolator si = new SincInterpolator();
    float[] g = new float[n1];
    si.setUniform(n1,1.0,0.0,f);
    si.interpolate(n1,u,g);
    return g;
  }
  
//  public float[] getU(float[] vpvs) {
//    int n1 = vpvs.length;
//    double ft = _s.getFirst();
//    double dt = _s.getDelta();
//    float[] u = new float[n1];
//    float s = _a/(2.0f*_b);
////    float s = 1.0f/(2.0f*_b);
//    for (int i1=0; i1<n1; i1++) {
//      double t = ft+i1*dt;
//      u[i1] = (float)(0.5*t + s*sin(2.0*PI*_b*t));
////      u[i1] = (float)(_b*t+_a*sin(2.0*PI*_b*t))*s;
//    }
//    return u;
//  }
  
  public float[] getU(float[] vpvs) {
    int n1 = vpvs.length;
    float dt = (float)_s.getDelta();
    dt = 1.0f;
    float[] u = new float[n1];
    u[0] = 0.5f*(vpvs[0]*dt-1.0f);
    for (int i1=1; i1<n1; i1++) {
      u[i1] = u[i1-1]+(0.5f*(vpvs[i1]*dt-1.0f));
    }
    return mul(u,0.5f);
  }
  
  public float[] computeVpVs(float[] u) {
    int n1 = u.length;
    double dt = _s.getDelta();
    dt = 1.0f;
    float[] vpvs = new float[n1];
    vpvs[0] = 1.0f;
    for (int i1=1; i1<n1; i1++) {
      vpvs[i1] = (float)(1.0 + 2.0*((u[i1]-u[i1-1])/dt));
    }
    return vpvs;
  }
  
  private float _a;
  private float _b;
  private Sampling _s;

}
