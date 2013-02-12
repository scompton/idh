package warp;

import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.interp.CubicInterpolator;
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
    SincInterp si = new SincInterp();
    float[] pp = new float[npp];
    for (int i=0; i<npp; ++i) {
      double y = i;
      pp[i] = si.interpolate(nps,1.0,0.0,ps,y+u[i]);
    }
    return pp;
  }
  
  public static float[][] shift2(int n2, float[] x, double scale, double f) {
    int n1 = x.length;
    float[][] s = new float[n2][n1];
    float[] r = new float[n1];
    ramp(0.0f,1.0f,r);
    CubicInterpolator ci = 
        new CubicInterpolator(CubicInterpolator.Method.LINEAR,r,x);
    for (int i2=0; i2<n2; i2++) {
      float s2 = (float)cos(2.0*PI*f*i2);
      float s2s = (float)(s2*scale);
      assert(abs(s2s)>=0 && abs(s2s)<=scale):"|s2s|="+abs(s2s)+", scale="+scale;
      for (int i1=0; i1<n1; i1++) {
        s[i2][i1] = ci.interpolate(s2s+i1);
      }
    }
    return s;
  }
  
  public static float[][][] shift3(
      int n3, float[][] x, double scale, double f)
  {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][][] s = new float[n3][n2][n1];
    float[] r = new float[n1];
    ramp(0.0f,1.0f,r);
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        CubicInterpolator ci = 
            new CubicInterpolator(CubicInterpolator.Method.LINEAR,r,x[i2]);
        float s2 = (float)sin(2.0*PI*f*i2);
        for (int i1=0; i1<n1; i1++) {
          s[i3][i2][i1] = ci.interpolate((float)(s2*scale)+i1);
        }
      }
    }
    return s;
  }
  
//  public static double uy(double y, SincInterpolator si) {
//    double uy = 0.0;
//    double up;
//    do {
//      up = uy;
//      uy = si.interpolate(y-uy);
//    } while (abs(uy-up)>0.0001);
//    return uy;
//  }
  
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
