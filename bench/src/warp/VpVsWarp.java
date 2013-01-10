package warp;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.util.Check;

public class VpVsWarp {

  /**
   * Computes Vp/Vs values at every time sample. The values
   * correspond to a scaled cosine function with a constant
   * frequency b.
   * @param a scale factor between -1.0 and 1.0.
   * @param b frequency in Hz.
   * @param s Sampling.
   * @return synthetic Vp/Vs.
   */
  public static float[] getVpVs(float a, float b, Sampling s) {
    Check.argument(a>=-1.0 & a<=1.0,"a>=-1.0 & a<=1.0");
    int n1 = s.getCount();
    double ft = s.getFirst();
    double dt = s.getDelta();
    float[] vpvs = new float[n1];
    for (int i1=0; i1<n1; i1++) {
      double t = ft+i1*dt;
      vpvs[i1] = 2.0f + a*(float)cos(2.0*PI*b*t);
    }
    return vpvs;
  }
  
  /**
   * Computes the shifts u, from the vpvs values.
   * @param vpvs array of computed vpvs values.
   * @param s Sampling.
   * @return shift from synthetic Vp/Vs.
   */
  public static float[] getU(float[] vpvs, Sampling s) {
    int n1 = vpvs.length;
    float dt = (float)s.getDelta();
    float[] u = new float[n1];
//    u[0] = (vpvs[0]-1.0f)/2.0f;
    u[0] = 0.0f;
    for (int i1=1; i1<n1; i1++) {
//      u[i1] = u[i1-1] + (vpvs[i1]-1.0f)/2.0f*dt;
      u[i1] = u[i1-1] + (vpvs[i1]-1.0f)/2.0f;
    }
    return u;
  }

  // Comparison for getU method. The results of this should be the same.
  public static float[] computeU(float a, float b, Sampling s) {
    Check.argument(a>=-1.0 & a<=1.0,"a>=-1.0 & a<=1.0");
    int n1 = s.getCount();
    double ft = s.getFirst();
    double dt = s.getDelta();
    float[] u = new float[n1];
    float d = (float)(4.0*PI*b);
    for (int i1=0; i1<n1; i1++) {
      double t = ft+i1*dt;
      u[i1] = (float)(a*sin(2.0*PI*b*t)/d + t/2.0f);
    }
    return u;
  }

  public static float getAverage(float[] vpvs) {
    int n1 = vpvs.length;
    return sum(vpvs)/n1;
  }

  public static float[] warp(float[] ps, float[] vpvs, Sampling s) {
    int nps = ps.length;
    float vpvsAvg = getAverage(vpvs);
    int npp = (int)ceil(2.0f*nps/(vpvsAvg+1.0f));
    System.out.println("Average Vp/Vs: "+vpvsAvg+", npp="+npp+", nps="+nps);
    float[] u = getU(vpvs,s);
    SincInterpolator si = new SincInterpolator();
    si.setUniform(nps,1.0,0.0,ps);
    SincInterpolator siu = new SincInterpolator();
    siu.setUniform(u.length,1.0,0.0,u);
    float[] pp = new float[npp];
    for (int i=0; i<npp; ++i) {
      double y = i;
//      double x = y-uy(y,siu);
//      ps[i] = si.interpolate(x);
      pp[i] = si.interpolate(y+u[i]);
    }
    return pp;
  }
  
  public static double uy(double y, SincInterpolator si) {
    double uy = 0.0;
    double up;
    do {
      up = uy;
      uy = si.interpolate(y-uy);
    } while (abs(uy-up)>0.0001);
    return uy;
  }
  
//  public static float[] warp(float[] f, float[] vpvs) {
//    int n1 = f.length;
//    float[] u = getU(vpvs);
//    SincInterpolator si = new SincInterpolator();
//    float[] g = new float[n1];
//    si.setUniform(n1,1.0,0.0,f);
//    si.interpolate(n1,u,g);
//    return g;
//  }
  
}
