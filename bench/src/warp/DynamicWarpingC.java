package warp;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.interp.BicubicInterpolator2;
import edu.mines.jtk.interp.BilinearInterpolator2;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Smooth dynamic warping for converted wave images.
 */
public class DynamicWarpingC {

  public enum Interp {
    LINEAR,
    MONOTONIC,
    SPLINE
  }

  /**
   * Constructor for dynamic warping of a PP and PS trace. The 
   * PS trace is warped to the PP trace. Note that the PP and
   * PS array do not need to be the same length. The maximum
   * shift is determined by the given {@code vpvsAvg}.
   * @param pp 
   * @param ps
   * @param vpvsAvg
   */
  public DynamicWarpingC(float[] pp, float[] ps, float vpvsAvg) {
    Check.argument(vpvsAvg>=1,"vpvsAvg>=1");
    _scale = 2.0f/(vpvsAvg+1.0f);
    int n1pp = pp.length;
    int n1ps = ps.length;
    _fr = 1;
    int[] el = computeErrorLengths(n1pp,n1ps,0);
    _nel = el[1];
    _ne1 = el[0];
    _n1pp = n1pp;
    _n1ps = n1ps;
    _pp1 = pp;
    _ps1 = ps;
    print("PP/PS Traces="+_n2+", PP Sample Length="+_n1pp+
        ", PP Sample Length for Warping="+_ne1+", PS Sample Length="+_n1ps+
        ", Number of Lags="+_nel);
    _si = new SincInterp();
  }

  /**
   * Constructor for dynamic warping of 2D PP and PS images. The 
   * PS traces are warped to the PP traces. Note that the number of 
   * PP and PS traces must be the same, but the traces do not need
   * to be the same length. The maximum shift is determined by the
   * given {@code vpvsAvg}.
   * @param pp 
   * @param ps
   * @param vpvsAvg
   */
  public DynamicWarpingC(float[][] pp, float[][] ps, float vpvsAvg) {
    Check.argument(vpvsAvg>=1,"vpvsAvg>=1");
    Check.argument(pp.length==ps.length,"pp.length==ps.length");
    _scale = 2.0f/(vpvsAvg+1.0f);
    int n1pp = pp[0].length;
    int n1ps = ps[0].length;
    _fr = 1;
    int[] el = computeErrorLengths(n1pp,n1ps,0);
    _nel = el[1];
    _ne1 = el[0];
    _n1pp = n1pp;
    _n1ps = n1ps;
    _n2 = pp.length;
    _pp2 = pp;
    _ps2 = ps;
    print("PP/PS Traces="+_n2+", PP Sample Length="+_n1pp+
        ", PP Sample Length for Warping="+_ne1+", PS Sample Length="+_n1ps+
        ", Number of Lags="+_nel);
    _si = new SincInterp();
  }

  /**
   * Constructor for dynamic warping of 3D PP and PS images. The 
   * PS traces are warped to the PP traces. Note that the number of 
   * PP and PS traces must be the same, but the traces do not need 
   * to be the same length. The maximum shift is determined by the
   * given {@code vpvsAvg}. 
   * @param pp 
   * @param ps
   * @param vpvsAvg
   */
  public DynamicWarpingC(float[][][] pp, float[][][] ps, float vpvsAvg) {
    Check.argument(vpvsAvg>=1,"vpvsAvg>=1");
    Check.argument(pp.length==ps.length,"pp.length==ps.length");
    Check.argument(pp[0].length==ps[0].length,"pp[0].length==ps[0].length");
    _scale = 2.0f/(vpvsAvg+1.0f);
    int n1pp = pp[0][0].length;
    int n1ps = ps[0][0].length;
    _fr = 1;
    int[] el = computeErrorLengths(n1pp,n1ps,0);
    _nel = el[1];
    _ne1 = el[0];
    _n1pp = n1pp;
    _n1ps = n1ps;
    _n2 = pp[0].length;
    _n3 = pp.length;
    _pp3 = pp;
    _ps3 = ps;
    print("PP/PS Ensembles="+_n3+", PP/PS Traces="+_n2+
        ", PP Sample Length="+_n1pp+", PP Sample Length for Warping="+_ne1+
        ", PS Sample Length="+_n1ps+", Number of Lags="+_nel);
    _si = new SincInterp();
  }

  /**
   * Returns the number of PP samples used for alignment errors.
   * This length is computed from the average Vp/Vs given in the 
   * constructor. 
   * @return the number of PP samples used for alignment errors.
   */
  public int getPPErrorLength() {
    return _ne1;
  }

  /**
   * Returns the number of lags. The number of lags is computed from
   * the average Vp/Vs given in the constructor. It is the difference
   * between the scaled PP trace and PS trace length.
   * @return the number of lags.
   */
  public int getNumberOfLags() {
    return _nel;
  }

  /**
   * Find shifts for 1D traces. For the input trace {@code PP[ne1]}, shifts 
   * {@code u} are computed on a subsampled grid such that the computed shifts 
   * are {@code u[ng1]}. This length matches the length of array {@code g1} 
   * and the indices of the subsampled shifts are specified by contents of 
   * this array. The contents of the {@code g1} array are rounded to the 
   * nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input trace, that is 
   * {@code ui[ne1]}.
   * </p>
   * NOTE: The length {@code ne1} is computed from the average Vp/Vs ratio 
   * given in the constructor. The value can be returned from the method 
   * {@link #getPPErrorLength()}.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension. Note for coarse grid 
   *  intervals dg1, if r1max*dg1<1 then the maximum slope is zero.
   * @param g1 array of size [ng1] specifying first dimension sparse grid
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @return shifts for 1D traces.
   */
  public float[] findShifts(
      final float r1min, final float r1max, final float[] g1, Interp interp1) 
  {
    int ng1 = g1.length;
    final int[] g1i = new int[ng1];
    for (int ig1=0; ig1<ng1; ig1++) {
      g1i[ig1] = (int)(g1[ig1]+0.5f);
    }
    float[][] e = computeErrors();
    float[][][] dm = accumulateForwardSparse(e,r1min,r1max,g1i);
    float[] u = backtrackReverse(dm[0],dm[1]);
    return interpolate(g1,u,interp1);
  }

  /**
   * Find shifts for 2D images. For the input image {@code PP[n2][ne1]}, shifts 
   * {@code u} are computed on a subsampled grid such that the computed shifts 
   * are {@code u[ng2][ng1]}. These lengths match the length of arrays 
   * {@code g1,g2} and the indices of the subsampled shifts are specified by 
   * contents of these arrays. The contents of the {@code g1} array are rounded
   * to the nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is 
   * {@code ui[n2][ne1]}.
   * </p>
   * NOTE: The fast dimension length {@code ne1} is computed from the average
   * Vp/Vs ratio given in the constructor. The value can be returned from the
   * method {@link #getPPErrorLength()}.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension. Note for coarse grid 
   *  intervals dg1, if r1max*dg1<1 then the maximum slope is zero.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param r2min minimum slope in second dimension.
   * @param r2max maximum slope in second dimension. Note for coarse grid 
   *  intervals dg2, if r2max*dg2<1 then the maximum slope is zero.
   * @param g2 array of size [ng2] specifying second dimension sparse grid 
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp2 interpolation method for the i2 (slow) dimension.
   * @return shifts for 2D images.
   */
  public float[][] findShifts(
      final float r1min, final float r1max, final float[][] g1,
      final float r2min, final float r2max, final int[]   g2,
      Interp interp1, Interp interp2)
  {
    Check.argument(_n2==g1.length,"_n2==g1.length");
    int ng1 = g1[0].length;
    final int[][] g1i = new int[_n2][ng1];
    for (int i2=0; i2<_n2; i2++) {
      for (int i1=0; i1<ng1; i1++) {
        g1i[i2][i1] = (int)(g1[i2][i1]+0.5f);
      }
    }
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][] es1 = smoothErrors1(_pp2,_ps2,r1min,r1max,g1i);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);
    
    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][] es = smoothErrors2(es1,r2min,r2max,g2);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es);
    
    final int ng2 = es.length;
    final float[][] u = new float[ng2][];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(es[i2],g1i[g2[i2]],r1min,r1max);
      u[i2] = backtrackReverse(dm[0],dm[1]);
    }});
    for (int i2=0; i2<ng2; i2++) {
      for (int i1=1; i1<g1[0].length; i1++) {
        float n = u[i2][i1] - u[i2][i1-1];
        float d = g1i[g2[i2]][i1] - g1i[g2[i2]][i1-1];
        float r = n/d;
        assert r>=r1min && r<=r1max:"n="+n+", d="+d+", r="+r;
      }
    }
    return interpolate(g1,g2,u,interp1,interp2);
  }

  /**
   * Find shifts for 3D images. For the input image {@code PP[n3][n2][ne1]},
   * shifts {@code u} are computed on a subsampled grid such that the computed
   * shifts are {@code u[ng3][ng2][ng1]}. These lengths match the length of
   * arrays {@code g1,g2,g3} and the indices of the subsampled shifts are 
   * specified by contents of these arrays. The contents of the {@code g1} 
   * array are rounded to the nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is 
   * {@code ui[n3][n2][ne1]}.
   * </p>
   * NOTE: The fast dimension length {@code ne1} is computed from the average
   * Vp/Vs ratio given in the constructor. The value can be returned from the
   * method {@link #getPPErrorLength()}.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension. Note for coarse grid 
   *  intervals dg1, if r1max*dg1<1 then the maximum slope is zero.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param r2min minimum slope in second dimension.
   * @param r2max maximum slope in second dimension. Note for coarse grid 
   *  intervals dg2, if r2max*dg2<1 then the maximum slope is zero.
   * @param g2 array of size [ng2] specifying second dimension sparse grid 
   *  locations.
   * @param r3min minimum slope in third dimension.
   * @param r3max maximum slope in third dimension. Note for coarse grid 
   *  intervals dg3, if r3max*dg3<1 then the maximum slope is zero.
   * @param g3 array of size [ng3] specifying third dimension sparse grid 
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp23 interpolation method for the i2 (middle) and i3
   *  (slow) dimensions.
   * @return shifts for 3D images.
   */
  public float[][][] findShifts(
      final float r1min, final float r1max, final float[][][] g1,
      final float r2min, final float r2max, final int[] g2,
      final float r3min, final float r3max, final int[] g3,
      Interp interp1, Interp interp23)
  {
    Check.argument(_n3==g1.length,"_n3==g1.length");
    Check.argument(_n2==g1[0].length,"_n2==g1[0].length");
    int ng1 = g1[0][0].length;
    final int[][][] g1i = new int[_n3][_n2][ng1];
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        for (int i1=0; i1<ng1; i1++) {
          g1i[i3][i2][i1] = (int)(g1[i3][i2][i1]+0.5f);
        }
      }
    }
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][][] es1 = smoothErrors1(_pp3,_ps3,r1min,r1max,g1i);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);

    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][][] es12 = smoothErrors2(es1,r2min,r2max,g2);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es12);
    
    s.restart();
    print("Smoothing 3rd dimension...");
    final float[][][][] es = smoothErrors3(es12,r3min,r3max,g3);
    normalizeErrors(es);
    print("Finished 3rd dimension smoothing in "+s.time()+" seconds");
    
    final int ng2 = es[0].length;
    final int ng3 = es.length;
    final float[][][] u = new float[ng3][ng2][];
    Parallel.loop(ng3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<ng2; i2++) {
        float[][][] dm = 
            accumulateForward(es[i3][i2],g1i[g3[i3]][g2[i2]],r1min,r1max);
        u[i3][i2] = backtrackReverse(dm[0],dm[1]);  
      }
    }});
    for (int i3=0; i3<ng3; i3++) {
      for (int i2=0; i2<ng2; i2++) {
        for (int i1=1; i1<ng1; i1++) {
          float n = u[i3][i2][i1] - u[i3][i2][i1-1];
          float d = g1i[g3[i3]][g2[i2]][i1] - g1i[g3[i3]][g2[i2]][i1-1];
          float r = n/d;
          assert r>=r1min && r<=r1max:"n="+n+", d="+d+", r="+r;
        }
      }
    }
    return interpolate(g1,g2,g3,u,interp1,interp23);
  }

  /**
   * Applies the shifts {@code u} to the PS trace.
   * @param u the shifts that warp the PS trace to the PP trace.
   * @return the warped PS trace.
   */
  public float[] applyShifts(float[] u) {
    int nu = u.length;
    int num = nu-1;
    float[] h = new float[_n1pp];
    for (int iu=0; iu<nu; ++iu) {
      h[iu] = _si.interpolate(_n1ps,1.0,0.0,_ps1,iu+u[iu]);
    }
    for (int i1=nu; i1<_n1pp; ++i1) {
      h[i1] = _si.interpolate(_n1ps,1.0,0.0,_ps1,i1+u[num]);
    }
    return h;
  }
  
  /**
   * Applies the shifts {@code u} to the PS image.
   * @param u the shifts that warp the PS image to the PP image.
   * @return the warped PS image.
   */
  public float[][] applyShifts(float[][] u) {
    final int n1u = u[0].length;
    final int n1um = n1u-1;
    final float[][] uf = u;
    final float[][] hf = new float[_n2][_n1pp];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1u; ++i1) {
        hf[i2][i1] = _si.interpolate(_n1ps,1.0,0.0,_ps2[i2],i1+uf[i2][i1]);
      }
      for (int i1=n1u; i1<_n1pp; ++i1) {
        hf[i2][i1] = _si.interpolate(_n1ps,1.0,0.0,_ps2[i2],i1+uf[i2][n1um]);
      }
    }});
    return hf;
  }
  
  /**
   * Applies the shifts {@code u} to the PS image.
   * @param u the shifts that warp the PS image to the PP image.
   * @return the warped PS image.
   */
  public float[][][] applyShifts(float[][][] u) {
    final int n1u = u[0][0].length;
    final int n1um = n1u-1;
    final float[][][] uf = u;
    final float[][][] hf = new float[_n3][_n2][_n1pp];
    Parallel.loop(_n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; i2++) {
        for (int i1=0; i1<n1u; i1++) {
          hf[i3][i2][i1] = _si.interpolate(_n1ps,1.0,0.0,_ps3[i3][i2],
              i1+uf[i3][i2][i1]);
        }
        for (int i1=n1u; i1<_n1pp; i1++) {
          hf[i3][i2][i1] = _si.interpolate(_n1ps,1.0,0.0,_ps3[i3][i2],
              i1+uf[i3][i2][n1um]);
        }  
      }
    }});
    return hf;
  }
  
  /**
   * Computes the NRMS value after applying shifts {@code u}
   * to the PS trace. This is the RMS of the difference divided
   * by the average RMS of the PP and PS trace.
   * @param u the shifts to warp the PS trace to the PP trace.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public float computeNrms(float[] u) {
    float[] g = applyShifts(u);
    float[] fs = new float[_ne1];
    float[] gs = new float[_ne1];
    for (int i1=0; i1<_ne1; i1++) {
      fs[i1] = _pp1[i1];
      gs[i1] = g[i1];
    }
    float[] d = sub(gs,fs);
    float frms = sqrt(sum(mul(fs,fs))/(_ne1));  
    float grms = sqrt(sum(mul(gs,gs))/(_ne1));
    float drms = sqrt(sum(mul( d, d))/(_ne1));
    float nrms = (2.0f*drms)/(frms+grms);
    return nrms;
  }
  
  /**
   * Computes the NRMS value after applying shifts {@code u}
   * to the PS image. This is the RMS of the difference divided
   * by the average RMS of the PP and PS image.
   * @param u the shifts to warp the PS image to the PP image.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public float computeNrms(float[][] u) {
    float[][] g = applyShifts(u);
    float[][] fs = new float[_n2][_ne1];
    float[][] gs = new float[_n2][_ne1];
    for (int i2=0; i2<_n2; i2++) {
      for (int i1=0; i1<_ne1; i1++) {
        fs[i2][i1] = _pp2[i2][i1];
        gs[i2][i1] = g[i2][i1];
      }
    }
    float[][] d = sub(gs,fs);
    float frms = sqrt(sum(mul(fs,fs))/(_ne1*_n2));  
    float grms = sqrt(sum(mul(gs,gs))/(_ne1*_n2));
    float drms = sqrt(sum(mul( d, d))/(_ne1*_n2));
    float nrms = (2.0f*drms)/(frms+grms);
    return nrms;
  }
  
  /**
   * Computes the NRMS value after applying shifts {@code u}
   * to the PS image. This is the RMS of the difference divided
   * by the average RMS of the PP and PS image.
   * @param u the shifts to warp the PS image to the PP image.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public float computeNrms(float[][][] u) {
    float[][][] g = applyShifts(u);
    float[][][] fs = new float[_n3][_n2][_ne1];
    float[][][] gs = new float[_n3][_n2][_ne1];
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        for (int i1=0; i1<_ne1; i1++) {
          fs[i3][i2][i1] = _pp3[i3][i2][i1];
          gs[i3][i2][i1] = g[i3][i2][i1];
        }
      }
    }
    float[][][] d = sub(gs,fs);
    float frms = sqrt(sum(mul(fs,fs))/(_ne1*_n2*_n3));  
    float grms = sqrt(sum(mul(gs,gs))/(_ne1*_n2*_n3));
    float drms = sqrt(sum(mul( d, d))/(_ne1*_n2*_n3));
    float nrms = (2.0f*drms)/(frms+grms);
    return nrms;
  }
  
  /**
   * Computes an array of VpVs ratios an array of shifts u.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values. 
   */
  public static float[] vpvs(float[] u) {
    int n = u.length;
    int nm1 = n-1;
    float[] vpvs = new float[n];
    vpvs[ 0 ] = 1.0f + 2.0f*(u[ 1 ]-u[  0  ]); // at i1=0, forward difference
    vpvs[nm1] = 1.0f + 2.0f*(u[nm1]-u[nm1-1]); // at i1=nm1, backward difference
    for (int i1=1; i1<nm1; ++i1)
      vpvs[i1] = 1.0f + (u[i1+1]-u[i1-1]);
    return vpvs;
  }
  
  /**
   * Computes an array of VpVs ratios an array of shifts u.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values. 
   */
  public static float[][] vpvs(float[][] u) {
    int n1 = u[0].length;
    int n2 = u.length;
    int n1m1 = n1-1;
    float[][] vpvs = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      // At i1=0, forward difference. At i1=nm1, backward difference.
      vpvs[i2][  0 ] = 1.0f + 2.0f*(u[i2][  1 ]-u[i2][   0  ]); 
      vpvs[i2][n1m1] = 1.0f + 2.0f*(u[i2][n1m1]-u[i2][n1m1-1]);
      for (int i1=1; i1<n1m1; ++i1)
        vpvs[i2][i1] = 1.0f+(u[i2][i1+1]-u[i2][i1-1]);
    }
    return vpvs;
  }
  
  /**
   * Computes an array of VpVs ratios an array of shifts u.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values. 
   */
  public static float[][][] vpvs(float[][][] u) {
    int n1 = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    int n1m1 = n1-1;
    float[][][] vpvs = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        // At i1=0, forward difference. At i1=nm1, backward difference.
        vpvs[i3][i2][  0 ] = 1.0f + 2.0f*(u[i3][i2][  1 ]-u[i3][i2][   0  ]); 
        vpvs[i3][i2][n1m1] = 1.0f + 2.0f*(u[i3][i2][n1m1]-u[i3][i2][n1m1-1]);
        for (int i1=1; i1<n1m1; ++i1)
          vpvs[i3][i2][i1] = 1.0f+(u[i3][i2][i1+1]-u[i3][i2][i1-1]);
      }  
    }
    return vpvs;
  }
  
  public static float[] extrapolate(float[] f, int n1) {
    int n1f = f.length;
    float v = f[n1f-1];
    float[] ef = new float[n1];
    for (int i1=0; i1<n1f; i1++) {
      ef[i1] = f[i1];
    }
    for (int i1=n1f; i1<n1; i1++) {
      ef[i1] = v;
    }
    return ef;
  }
  
  public static float[][] extrapolate(float[][] f, int n1) {
    int n1f = f[0].length;
    int n2 = f.length;
    float[][] ef = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      float v = f[i2][n1f-1];
      for (int i1=0; i1<n1f; i1++) {
        ef[i2][i1] = f[i2][i1];
      }
      for (int i1=n1f; i1<n1; i1++) {
        ef[i2][i1] = v;
      }
    }
    return ef;
  }
  
  public static float[][][] extrapolate(float[][][] f, int n1) {
    int n1f = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] ef = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        float v = f[i3][i2][n1f-1];
        for (int i1=0; i1<n1f; i1++) {
          ef[i3][i2][i1] = f[i3][i2][i1];
        }
        for (int i1=n1f; i1<n1; i1++) {
          ef[i3][i2][i1] = v;
        }
      }  
    }
    return ef;
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // for research and atypical applications
  
  /**
   * Find shifts for 2D images from averaged alignment errors. These
   * shifts are constant for all traces.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param g1 array of subsampled indices.  
   * @return shifts for 2D images from averaged alignment errors.
   */
  public float[] findShifts2(float r1min, float r1max, int[] g1) {
    final float[][] e = computeErrorsSum2();
    final float[][][] dm = accumulateForwardSparse(e,r1min,r1max,g1);
    return backtrackReverse(dm[0],dm[1]);
  }

  public float[] findShifts3(float r1min, float r1max, int[] g1) {
    float[][] e3Avg = computeErrorsSum3();
    final float[][][] dm = accumulateForwardSparse(e3Avg,r1min,r1max,g1);
    return backtrackReverse(dm[0],dm[1]);
  }

  public float[][] applyShifts(
      float[] u1, float[] uS, float[] ps1, float[] ps2)
  {
    int n1 = ps1.length;
    float[] ps1w = new float[n1]; 
    float[] ps2w = new float[n1];
    int nu = u1.length;
    int num = nu-1;
    for (int iu=0; iu<nu; ++iu) {
      ps1w[iu] = _si.interpolate(n1,1.0,0.0,ps1,iu+u1[iu]);
      ps2w[iu] = _si.interpolate(n1,1.0,0.0,ps2,iu+u1[iu]+uS[iu]);
    }
    for (int i1=nu; i1<n1; ++i1) {
      ps1w[i1] = _si.interpolate(n1,1.0,0.0,ps1,i1+u1[num]);
      ps2w[i1] = _si.interpolate(n1,1.0,0.0,ps2,i1+u1[num]+uS[num]);
    }
    return new float[][]{ps1w, ps2w};
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u
   * using a backward difference approximation.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values. 
   */
  public static float[] vpvsBd(float[] u) {
    int n = u.length;
    float[] vpvs = new float[n];
    vpvs[0] = 1.0f + 2.0f*(u[1]-u[0]); // at i1=0, forward difference
    for (int i1=1; i1<n; ++i1)
      vpvs[i1] = 1.0f + 2.0f*(u[i1]-u[i1-1]);
    return vpvs;
  }

  /**
   * Compute alignment errors for 1D traces.
   * @return alignment errors for 1D traces.
   */
  public float[][] computeErrors() {
    float[][] e = new float[_ne1][_nel];
    computeErrors(_pp1,_ps1,e);
    normalizeErrors(e);
    return e;
  }

  /**
   * Compute alignment errors for 2D traces.
   * @return alignment errors for 2D traces.
   */
  public float[][][] computeErrors2() {
    final float[][][] e = new float[_n2][_ne1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrors(_pp2[i2],_ps2[i2],e[i2]);  
    }});
    return e;
  }

  /**
   * Returns alignment errors for 2D images. Alignment errors
   * at every trace are the sum of 2*w nearby traces. 
   * @param w the width to extend left and right from each
   *  trace 
   * @return Alignment errors for 2D image with Alignment errors
   *  at every trace are the sum of 2*w nearby traces.
   */
  public float[][][] computeErrors2Near(final int w) {
    final float[][][] e = new float[_n2][_ne1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      e[i2] = computeErrorsNear(i2,w);
    }});
    return e;
  }

  /**
   * Compute summed alignment errors for 2D images.
   * @return summed alignment errors for 2D images. 
   */
  public float[][] computeErrorsSum2() {
    float[][] e = Parallel.reduce(_n2,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(_pp2[i2],_ps2[i2],e);
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    return e;
  }

  /**
   * Compute summed alignment errors for 3D images.
   * @return summed alignment errors for 3D images. 
   */
  public float[][] computeErrorsSum3() {
    float[][] e = Parallel.reduce(_n2*_n3,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i23) {
      int i2 = i23%_n2;
      int i3 = i23/_n2;
      float[][] e = new float[_ne1][_nel];
      computeErrors(_pp3[i3][i2],_ps3[i3][i2],e);
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    return e;
  }

  public float[][][] accumulateForward(
      float[][] e, int[] g, float rmin, float rmax)
  {
    float[][] d = new float[e.length][e[0].length];
    float[][] m = new float[e.length][e[0].length];
    accumulateFromSparse(1,e,d,m,g,rmin,rmax);
    return new float[][][]{d,m};
  }

  public static float[][][] accumulateForwardSparse(
      float[][] e, float rmin, float rmax, int[] g)
  {
    int nl = e[0].length;
    int ng = g.length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    accumulateSparse(1,rmin,rmax,g,e,d,m);
    return new float[][][]{d,m};
  }
  
  public static float[][][] accumulateReverseSparse(
      float[][] e, float rmin, float rmax, int[] g)
  {
    int nl = e[0].length;
    int ng = g.length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    accumulateSparse(-1,rmin,rmax,g,e,d,m);
    return new float[][][]{d,m};
  }
  
  public float[] backtrackReverse(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrack(-1,_shifts1,d,m,u);
    return u;
  }
  
  public float[] backtrackForward(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrack(1,_shifts1,d,m,u);
    return u;
  }

  public static float[][][] shiftVolume(
      float[] u1, float[] uS, int fr, float[][][] e) {
    int n1 = e.length;
    int nl1 = e[0].length;
    int nlS = e[0][0].length;
    Check.argument(n1==u1.length, "n1==u1.length");
    Check.argument(n1==uS.length, "n1==uS.length");
    float[][][] sv = zerofloat(nlS,nl1,n1);
    for (int i1=0; i1<n1; ++i1) {
      int il1 = (int)u1[i1];
      int ilS = (int)uS[i1]*fr;
      sv[i1][il1][ilS] = 1.0f;
    }
    return sv;
  }

  /**
   * Normalizes values to be in range [0,1].
   * @param e input/output array.
   */
  public static void normalizeErrors(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float emin = e[0][0];
    float emax = e[0][0];
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        float ei = e[i1][il];
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }
    shiftAndScale(emin,emax,e);
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][] e) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] ef = e;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          float ei = ef[i2][i1][il];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }
  
  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][][] e) {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    MinMax mm = Parallel.reduce(n3,new Parallel.ReduceInt<MinMax>() {
      public MinMax compute(int i3) {
        float emin =  Float.MAX_VALUE;
        float emax = -Float.MAX_VALUE;
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            for (int il=0; il<nl; ++il) {
              float ei = ef[i3][i2][i1][il];
              if (ei<emin) emin = ei;
              if (ei>emax) emax = ei;
            }
          }  
        }
        return new MinMax(emin,emax);
      }
      public MinMax combine(MinMax mm1, MinMax mm2) {
        return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
      }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

  public static void normalize(
      float[][][] f, final float nmin, final float nmax) 
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] ff = f;
    final float vmin = min(f);
    final float vmax = max(f);
    final float range = vmax-vmin;
    final float nrange = nmax-nmin;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float vi = ff[i3][i2][i1];
          ff[i3][i2][i1] = nrange*(vi-vmin)/range + nmin;
        }
      }
    }});
  }
  
  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n1][nl] of errors.
   * @return transposed array[nl][n1] of errors.
   */
  public static float[][] transposeLag(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] t = new float[nl][n1];
    for (int il=0; il<nl; ++il) {
      for (int i1=0; i1<n1; ++i1) {
        t[il][i1] = e[i1][il];
      }
    }
    return t;
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n2][n1][nl] of errors.
   * @return transposed array[nl][n2][n1] of errors.
   */
  public static float[][][] transposeLag(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[nl][n2][n1];
    for (int il=0; il<nl; ++il) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[il][i2][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }
  
  public static float[][][] transposeLag12(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[n2][nl][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int il=0; il<nl; ++il) {
        for (int i1=0; i1<n1; ++i1) {
          t[i2][il][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }
  
  public static float[][][] transposeLag23(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[n1][n2][nl];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        t[i1][i2] = e[i2][i1];
      }
    }
    return t;
  }
  
  public static float[][][][] transposeLag12(float[][][][] e) {
    int nl = e[0][0][0].length;
    int n1 = e[0][0].length;
    int n2 = e[0].length;
    int n3 = e.length;
    float[][][][] t = new float[n3][n2][nl][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          for (int i1=0; i1<n1; ++i1) {
            t[i3][i2][il][i1] = e[i3][i2][i1][il];
          }
        }
      }  
    }
    return t;
  }
  
  public static float[][][][] getEmptySparseErrors(int n3) {
    return new float[n3][][][];
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nel; // number of lags
  private int _ne1; // number of alignment error samples
  private int _n1pp; // number of pp samples before scaling
  private int _n1ps; // number of pp samples before scaling
  private int _n2; // number of traces
  private int _n3; // number of ensembles
  private float[] _pp1; // pp trace
  private float[] _ps1; // ps trace
  private float[][] _pp2; // pp traces
  private float[][] _ps2; // ps traces
  private float[][][] _pp3; // pp traces
  private float[][][] _ps3; // ps traces
  private int _fr; // fractional shift factor (1.0/_fr is the shift interval)
  private Sampling _shifts1; // sampling of pp to ps1 shift values
  private Sampling _shiftsS; // sampling of ps1 to ps2 shift values
  private float _scale; // computed scaling for pp traces
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private SincInterp _si; // for warping with non-integer shifts
  
  /**
   * Interpolates subsampled shifts u[ng1] to uniformly sampled shifts
   * ui[_ne1].
   * @param g1 sparse grid indices.
   * @param u sparse shifts.
   * @param doLinear {@code true} for linear interpolation,
   *  {@code false} for cubic.
   * @return the interpolated shifts.
   */
  private float[] interpolate(float[] g1, float[] u, Interp interp1) {
    int ng1 = g1.length;
    float[] ui = new float[_ne1];
    CubicInterpolator.Method m1;
    switch (interp1) {
      case LINEAR:    m1 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m1 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m1 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }
    CubicInterpolator ci = new CubicInterpolator(m1,ng1,g1,u);
    for (int i=0; i<_ne1; i++)
      ui[i] = ci.interpolate(i);
    return ui;
  }
  
  /**
   * Interpolates subsampled shifts u[ng2][ng1] to uniformly sampled 
   * shifts ui[_n2][_ne1]. The locations of the subsampled shifts u 
   * are defined by the input arrays g1 and g2. The g2 array is
   * assumed to contain integer indices consistent for all n1 locations. 
   * That is, g2 indices are the same at every i1, for i1=0,...,n1-1.
   * </p>
   * The g1 coordinate array defines the subsampled locations of the
   * fastest dimension of shifts u, for all n2 indices.
   * </p>
   * The interpolation is done in two passes. First interpolation
   * is done in the second dimension where subsampling must be 
   * regular. The second pass does interpolation in the first
   * dimension where subsampling may be irregular.
   * @param g1 first dimension sparse grid coordinates.
   * @param g2 second dimension sparse grid indices.
   * @param u sparse shifts.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp2 interpolation method for the i2 (slow) dimension.
   * @return the interpolated shifts ui[_n2][_ne1].
   */
  private float[][] interpolate(
      float[][] g1, int[] g2, float[][] u, Interp interp1, Interp interp2)
  {
    int ng1 = g1[0].length;
    int ng2 = g2.length;
    float[] g2f = new float[ng2];
    for (int ig2=0; ig2<ng2; ig2++)
      g2f[ig2] = g2[ig2];
    
    float[][] ui = new float[_n2][_ne1];
    CubicInterpolator.Method m2;
    switch (interp2) {
      case LINEAR:    m2 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m2 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m2 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp2.toString()+" is not a recognized interpolation method.");
    }
    CubicInterpolator.Method m1;
    switch (interp1) {
      case LINEAR:    m1 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m1 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m1 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }
    
    // interpolate in the second dimension.
    float[] u2 = new float[ng2];
    float[][] ui2 = new float[_n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i2=0; i2<ng2; i2++)
        u2[i2] = u[i2][i1];
      CubicInterpolator ciu = new CubicInterpolator(m2,g2f,u2);
      for (int i2=0; i2<_n2; i2++)
        ui2[i2][i1] = ciu.interpolate(i2);
    }

    // interpolate in the first dimension.
    for (int i2=0; i2<_n2; i2++) {
      CubicInterpolator ci = new CubicInterpolator(m1,g1[i2],ui2[i2]);
      for (int i1=0; i1<_ne1; i1++)
        ui[i2][i1] = ci.interpolate(i1);
    }
    
    return ui;
  }
  
  /**
   * Interpolates subsampled shifts u[ng3][ng2][ng1] to uniformly sampled 
   * shifts ui[_n3][_n2][_ne1]. The locations of the subsampled shifts u are 
   * defined by the input arrays g1, g2, and g3. The indices in the g3 array 
   * must be the same at every i2, for i2=0,...,n2-1 and i1, for i1=0,...n1-1. 
   * The indices in the g2 array must be the same at every i3, for 
   * i3=0,...,n3-1, and i1, for i1=0,...,n1-1.
   * </p>
   * The g1 coordinate array defines the subsampled locations of the
   * fastest dimension of shifts u, for all n3 and n2 indices.
   * </p>
   * The interpolation is done in two passes. First interpolation
   * is done in the second and third dimension where subsampling
   * must be regular. The second pass does interpolation in the 
   * first dimension where subsampling may be irregular.
   * @param g1 first dimension sparse grid coordinates.
   * @param g2 second dimension sparse grid indices.
   * @param g3 second dimension sparse grid indices.
   * @param u sparse shifts.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp23 interpolation method for the i2 (middle) and i3
   *  (slow) dimensions.
   * @return the interpolated shifts ui[_n3][_n2][_ne1].
   */
  private float[][][] interpolate(
      float[][][] g1, int[] g2, int[] g3, float[][][] u, 
      Interp interp1, Interp interp23) 
  {
    int ng1 = g1[0][0].length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    Check.argument(
        ng1==u[0][0].length,"ng1==u[0][0].length: "+ng1+"=="+u[0][0].length);
    Check.argument(ng2==u[0].length,"ng2==u[0].length: "+ng2+"=="+u[0].length);
    Check.argument(ng3==u.length,"ng3==u.length: "+ng3+"=="+u.length);
    float[] g2f = new float[ng2];
    float[] g3f = new float[ng3];
    for (int ig2=0; ig2<ng2; ig2++)
      g2f[ig2] = g2[ig2];
    for (int ig3=0; ig3<ng3; ig3++)
      g3f[ig3] = g3[ig3];
    
    CubicInterpolator.Method m1;
    switch (interp1) {
      case LINEAR:    m1 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m1 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m1 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }
    
    BicubicInterpolator2.Method m23 = null;
    boolean doLinear;
    switch (interp23) {
      case LINEAR:    doLinear = true; 
                      break;
      case MONOTONIC: doLinear = false; 
                      m23 = BicubicInterpolator2.Method.MONOTONIC; 
                      break;
      case SPLINE:    doLinear = false; 
                      m23 = BicubicInterpolator2.Method.SPLINE; 
                      break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }

    // interpolate the second and third dimension.
    float[][] u23 = new float[ng3][ng2];
    float[][][] ui23 = new float[_n3][_n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i3=0; i3<ng3; i3++)
        for (int i2=0; i2<ng2; i2++)
          u23[i3][i2] = u[i3][i2][i1];
      if (doLinear) {
        BilinearInterpolator2 bli = new BilinearInterpolator2(g2f,g3f,u23);
        for (int i3=0; i3<_n3; i3++)
          for (int i2=0; i2<_n2; i2++)
            ui23[i3][i2][i1] = bli.interpolate(i2,i3);
      } else {
        BicubicInterpolator2 bci = 
            new BicubicInterpolator2(m23,m23,g2f,g3f,u23);
        for (int i3=0; i3<_n3; i3++)
          for (int i2=0; i2<_n2; i2++)
            ui23[i3][i2][i1] = bci.interpolate(i2,i3);
      }
    }

    float[][][] ui = new float[_n3][_n2][_ne1];
    // interpolate the first dimension
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        CubicInterpolator ci = 
            new CubicInterpolator(m1,g1[i3][i2],ui23[i3][i2]);
        for (int i1=0; i1<_ne1; i1++)
          ui[i3][i2][i1] = ci.interpolate(i1);
      }
    }
    return ui;
  }
  
  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

//  private void computeErrors(
//      float[] pp, float[] ps1, float[] ps2, float[][][] e)
//  {
//    int n1 = pp.length;
//    float e1m=0.0f,e2m=0.0f,e3m=0.0f; // last computed values for extrapolation
//
//    // Notes for indexing:
//    // 0 <= il1 < _nl1, where il1 is index for lag between pp and ps1
//    // 0 <= ilS < _nlS, where ilS is index for lag between ps1 and ps2    
//    // 0 <= i1 < n1, where i1 is index for sequence pp
//    // 0 <= j1 < n1, where j1 is index for sequence ps1
//    // 0 <= jS < n1, where jS is index for sequence ps2
//    // j1 = i1+_l1min, where _l1min is the minimum shift u1 where u2 = u1 + uS
//    // jS = j1+_lSmin, where _lSmin is the minimum shift uS where u2 = u1 + uS
//
//    // Compute errors for pp, ps1, and ps2.
//    for (int i1=0; i1<n1; ++i1) {
//      for (int il1=0,j1=i1+_l1min; il1<_nl1; ++il1,++j1) {
//        float e1;
//        if (j1<n1) {
//          e1 = error(pp[i1],ps1[j1]);
//          e1m = e1;
//        } else
//          e1 = e1m;
//        for (int ilS=0,jS=j1+_lSmin; ilS<_nlS; ++ilS,++jS) {
//          float e2,e3;
//          if (jS<n1) {
//            e2 = error(pp [i1],ps2[jS]);
//            e3 = error(ps1[j1],ps2[jS]);
//            e2m = e2;
//            e3m = e3;
//          } else {
//            e2 = e2m;
//            e3 = e3m;
//          }
//          e[i1][il1][ilS] = e1+e2+e3;
//        }
//      }
//    }
//  }

  /**
   * Compute the error array lengths [n1][nl]. n1 is scaled by the
   * vpvs ratio set in the constructor and nl is difference between
   * the given n1 and the scaled n1. In this way, we do not compute
   * errors outside of the possible bounds determined by the average
   * vpvs scalar. In other words, the scalar determines the maximum
   * possible shift. Alternatively a maxShift parameter can be set
   * directly if maxShift > 0.
   * @param npp
   * @param maxShift
   * @return
   */
  private int[] computeErrorLengths(int npp, int nps, int maxShift) {
    int nppMax = (int)(nps*_scale);
    int nl = nps-nppMax;
    print("nppMax="+nppMax+", nl="+nl);
    if (maxShift!=0) {
      nl = maxShift;
      nppMax = nps-nl;
    }
    nppMax = (nppMax>npp)?npp:nppMax;
    nl = (nl-1)*_fr+1;
    _shifts1 = new Sampling(nl,1.0/_fr,0.0);
    return new int[]{nppMax,nl};
  }
  
  /**
   * Compute the error array lengths [n1][nl1][nlS]. nl1 is the
   * difference between the given n1 and the scaled n1, where the
   * scalar is determined by the vpvs ratio set in the constructor.
   * nlS is set as a percentage of nl1 by the nl1P parameter. n1 is
   * then the difference the given n1 and nl1 and nlS.
   * <p>
   * In this way, we do not compute errors outside of the possible 
   * bounds determined by the average vpvs scalar. In other words,
   * the scalar determines the maximum possible shift.
   * @param n1
   * @param nl1P
   * @return
   */
  private int[] computeErrorLengths(int n1, float nl1P) {
    int nl1 = n1-(int)(n1*_scale);
    int nlS = (int)(nl1*nl1P);
    int ppN1Max = n1-nl1-nlS+2;
    nlS = (nlS-1)*_fr+1;
    System.out.println("ppN1Max="+ppN1Max+", nl1="+nl1+", nlS="+nlS);
    _shifts1 = new Sampling(nl1,1.0,0.0);
    _shiftsS = new Sampling(nlS,1.0/_fr,0.0);
    return new int[]{ppN1Max,nl1,nlS};
  }
  
  /**
   * Computes alignment errors for {@code f} and {@code g}. Note
   * that values of the {@code e} array are not replaced, but
   * added to. 
   * @param f
   * @param g
   * @param e
   */
  private void computeErrors(float[] f, float[] g, float[][] e) {
    int n1max = e.length;
    int nl1 = e[0].length;
    for (int i1=0; i1<n1max; ++i1) {
      for (int il1=0,j1=i1*_fr; il1<nl1; ++il1,++j1) {
        e[i1][il1] += error(f[i1],g[j1]);
      }
    }
  }
  
  /**
   * Returns alignment errors at index i2, from a sum of alignment
   * errors at nearby traces.
   * @param i2 the center index used for summing nearby alignment
   *  errors.
   * @param w the width to extend left and right from i2. 
   * @return summed alignment errors for index i2.
   */
  private float[][] computeErrorsNear(final int i2, final int w) {
    final float[][] e = new float[_ne1][_nel];
    final int n2m1 = _n2-1;
    int i2min = max(0,i2-w);
    int i2max = min(n2m1,i2+w);
    for (int j=i2min; j<=i2max; j++)
      computeErrors(_pp2[j],_ps2[j],e);
    return e;
  }
  
  private void computeErrors(
      float[] pp, float[] ps1, float[] ps2, float[][][] e)
  {
    int n1max = e.length;
    int nl1 = e[0].length;
    int nlS = e[0][0].length;

    // Compute errors for pp, ps1, and ps2.
    for (int i1=0; i1<n1max; ++i1) {
      for (int il1=0,j1=i1; il1<nl1; ++il1,++j1) {
        float e1 = error(pp[i1],ps1[j1]);
        for (int ilS=0,jS=j1*_fr; ilS<nlS; ++ilS,++jS) {
          float e2 = error(pp [i1],ps2[jS]);
          float e3 = error(ps1[j1],ps2[jS]);
          e[i1][il1][ilS] = e1+e2+e3;
        }
      }
    }
  }
  
  /**
   * Returns smooth alignment errors on the sparse grid defined
   * by the indices of g. 
   * @param e 2D array of alignment errors.
   * @param rmin minimum slope.
   * @param rmax maximum slope.
   * @param g sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [g.length][e[0].length].
   */
  public static float[][] smoothErrors(
      float[][] e, float rmin, float rmax, int[] g)      
  {
    int ng = g.length;
    int nel = e[0].length;
    float[][] ef = new float[ng][nel];
    float[][] er = new float[ng][nel];
    float[][] es = new float[ng][nel];
    accumulateSparse( 1,rmin,rmax,g,e,ef,null);
    accumulateSparse(-1,rmin,rmax,g,e,er,null);
    float scale = 1.0f/e.length;
    for (int i1=0; i1<ng; i1++) {
      for (int il=0; il<nel; il++) {
        es[i1][il] = scale*(ef[i1][il]+er[i1][il]-e[g[i1]][il]);
      }
    }
    return es;
  }
  
  private float[][][] smoothErrors1W(final float[][] pp, final float[][] ps,
      final float r1min, final float r1max, final int[] g1, final int w)
  {
    final int ng1 = g1.length;
    final float[][][] es1 = new float[_n2][ng1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = computeErrorsNear(i2,w);
      es1[i2] = smoothErrors(e,r1min,r1max,g1);
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second dimension. Alignment errors are
   * computed on the fly.
   * @param pp the PP image.
   * @param ps the PS image.
   * @param r1min minimum slope in the first dimension.
   * @param r1max maximum slope in the first dimension.
   * @param g1 first dimension sparse grid indices specified
   *  for all _n2 (size[_n2][])
   * @return smoothed alignment errors with size 
   *  [_n2][g1.length][_nel].
   */
  private float[][][] smoothErrors1(
      final float[][] pp, final float[][] ps,
      final float r1min, final float r1max, final int[][] g1)
  {
    final int ng1 = g1.length;
    final float[][][] es1 = new float[_n2][ng1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(pp[i2],ps[i2],e);
      es1[i2] = smoothErrors(e,r1min,r1max,g1[i2]);
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second and third dimension. Alignment
   * errors are computed on the fly.
   * @param pp the PP image.
   * @param ps the PS image.
   * @param r1min minimum slope in the first dimension.
   * @param r1max maximum slope in the first dimension.
   * @param g1 first dimension sparse grid indices specified
   *  for all _n3 and _n2 (size[_n3][_n2][]).
   * @return smoothed alignment errors with size 
   *  [_n3][_n2][g1.length][_nel].
   */
  private float[][][][] smoothErrors1(
      final float[][][] pp, final float[][][] ps,
      final float r1min, final float r1max, final int[][][] g1)
  {
    final int ng1 = g1.length;
    final float[][][][] es1 = new float[_n3][_n2][ng1][_nel];
    Parallel.loop(_n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; i2++) {
        float[][] e = new float[_ne1][_nel];
        computeErrors(pp[i3][i2],ps[i3][i2],e);
        es1[i3][i2] = smoothErrors(e,r1min,r1max,g1[i3][i2]);  
      }
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the second dimension.
   * Returned errors are sparse in the second dimension, and
   * unchanged in the first dimension.
   * @param e alignment errors.
   * @param r2min minimum slope in the second dimension.
   * @param r2max maximum slope in the second dimension.
   * @param g2 second dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [g2.length][e[0].length][e[0][0].length].
   */
  private static float[][][] smoothErrors2(
      final float[][][] e, 
      final float r2min, final float r2max, final int[] g2)
  {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final int ng2 = g2.length;
    final float[][][] es = new float[ng2][n1][nl]; // smoothed errors
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][]  e2 = new float[n2][nl]; // errors at index i1
      for (int i2=0; i2<n2; ++i2)
        e2[i2] = e[i2][i1];
      float[][] es2 = smoothErrors(e2,r2min,r2max,g2);
      for (int i2=0; i2<ng2; i2++)
        es[i2][i1] = es2[i2];
    }});
    return es;
  }

  /**
   * Returns alignment errors smoothed in the second dimension.
   * Returned errors are sparse in the second dimension, and
   * unchanged in the first dimension and third dimension.
   * @param e alignment errors.
   * @param r2min minimum slope in the second dimension.
   * @param r2max maximum slope in the second dimension.
   * @param g2 second dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [e.length][g2.length][e[0][0].length][e[0][0][0].length].
   */
  private static float[][][][] smoothErrors2(
      final float[][][][] e, 
      final float r2min, final float r2max, final int[] g2)
  {
    final int n3 = e.length;
    final float[][][][] es = new float[n3][][][]; // smoothed errors
    for (int i3=0; i3<n3; i3++)
      es[i3] = smoothErrors2(e[i3],r2min,r2max,g2);
    return es;
  }

  /**
   * Returns alignment errors smoothed in the third dimension.
   * Returned errors are sparse in the third dimension, and
   * unchanged in the first and second dimension.
   * @param e alignment errors.
   * @param r3min minimum slope in the third dimension.
   * @param r3max maximum slope in the third dimension.
   * @param g3 third dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [g3.length][e[0].length][e[0][0].length][e[0][0][0].length].
   */
  private static float[][][][] smoothErrors3(
      final float[][][][] e, 
      final float r3min, final float r3max, final int[] g3)
  {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final int ng3 = g3.length;
    final float[][][][] es = new float[ng3][n2][n1][nl]; // smoothed errors
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      for (int i2=0; i2<n2; i2++) {
        float[][]  e3 = new float[n3][nl]; // smooth errors at index i1,i2
        for (int i3=0; i3<n3; i3++)
          e3[i3] = e[i3][i2][i1];
        float[][] es3 = smoothErrors(e3,r3min,r3max,g3);
        for (int i3=0; i3<ng3; i3++)
          es[i3][i2][i1] = es3[i3];
      }
    }});
    return es;
  }
  
  /**
   * Accumulation for 3D alignment Errors
   * @param dir
   * @param b1
   * @param bS
   * @param e
   * @param d
   */
  private static void accumulate(
      int dir, int b1, int bS, float[][][] e, float[][][] d)
  {
    int nl1 = e[0].length;
    int nl1m1 = nl1-1;
    int nlS = e[0][0].length;
    int nlSm1 = nlS-1;
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    for (int il1=0; il1<nl1; ++il1)
      for (int ilS=0; ilS<nlS; ++ilS)
        d[ib][il1][ilS] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb1 = max(0,min(nim1,ii-is*b1));
      int jbS = max(0,min(nim1,ii-is*bS));
      for (int il1=0; il1<nl1; ++il1) {
        int il1c = il1+ic;
        il1c = (il1c==-1)?0:(il1c==nl1)?nl1m1:il1c;
        for (int ilS=0; ilS<nlS; ++ilS) {
          int ilSc = ilS+ic;
          ilSc = (ilSc==-1)?0:(ilSc==nlS)?nlSm1:ilSc;
          float dc1 = d[jb1][il1c][ilS ];
          float dcS = d[jbS][il1 ][ilSc];
          float dc1S= d[jbS][il1c][ilSc];
          float di  = d[ji ][il1 ][ilS ];
          for (int kb1=ji; kb1!=jb1; kb1-=is) {
            dc1 += e[kb1][il1c][ilS ];
          }
          for (int kbS=ji; kbS!=jbS; kbS-=is) {
            dcS += e[kbS][il1 ][ilSc];
            dc1S+= e[kbS][il1c][ilSc];
          }
          d[ii][il1][ilS] = min4(dc1,dcS,dc1S,di)+e[ii][il1][ilS];
        }
      }
    }
  }
  
  private static void accumulateSparse(
      int dir, float rmin, float rmax, int[] g,
      float[][] e, float[][] d, float[][] m)
  {
    int nl   = e[0].length;
    int ng   = g.length;
    int ngm1 = ng-1;
    int ibg = (dir>0)?0:ngm1; // beginning index
    int ieg = (dir>0)?ng:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int isp = ibg; // sparse grid index
    int ie = g[isp]; // error index
    
    // Initialize accumulation values
    for (int il=0; il<nl; ++il)
      d[isp][il] = e[ie][il];
    isp += is;
    // Loop over all sparse grid points.
    for (; isp!=ieg; isp+=is) {
      int ispm1 = isp-is; // previous sparse grid index
      ie = g[isp]; // new error index
      int je = g[ispm1]; // min error index for interpolation.
      int dg = ie-je; // sparse grid delta
      int kmin, kmax;
      if (dg>0) {
        kmin = (int) ceil(-rmax*dg);
        kmax = (int)floor(-rmin*dg);
      } else {
        kmin = (int) ceil(-rmin*dg);
        kmax = (int)floor(-rmax*dg);
      }
      assert kmin<=kmax : "kmin="+kmin+", kmax="+kmax;
      float[] dm = new float[nl];
      fill(Float.MAX_VALUE,d[isp]);
      // loop over all slope indices
      for (int k=kmin; k<=kmax; k++) {
        int ils = max(0,-k);
        int ile = min(nl,nl-k);
        for (int il=ils; il<ile; il++)
          dm[il] = d[ispm1][il+k] + e[ie][il];
        float r = (float)k/(float)dg; // slope
        if (r==0) { // zero slope, no interpolation necessary
          for (int x=je+is; x!=ie; x+=is)
            for (int il=ils; il<ile; il++)
              dm[il] += e[x][il]; 
        } else { // linearly interpolate
          for (int x=je+is; x!=ie; x+=is) {
            float ky = r*(ie-x);
            int k1 = (int)ky;
            if (ky<0.0f) --k1;
            int k2 = k1+1;
            float w1 = k2-ky;
            float w2 = 1.0f-w1;
            for (int il=ils; il<ile; il++)
              dm[il] += w1*e[x][k1+il]+w2*e[x][k2+il];
          }  
        }
        // update previous errors and record moves.
        for (int il=ils; il<ile; il++) {
          if (dm[il]<d[isp][il]) {
            d[isp][il] = dm[il];
            if (m!=null)
              m[ispm1][il] = k;
          }
        }
      }
    }
  }
  
  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   * @param m output array[ni][nl] of recorded moves.
   * @param g
   * @param rmin
   * @param rmax
   */
  private static void accumulateFromSparse(
      int dir, float[][] e, float[][] d, float[][] m, int[] g,
      float rmin, float rmax)
  {
    int nl = e[0].length;
    int nlm1 = nl-1;
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    float dmax = ni;
    int ii=ib;
    for (int il=0; il<nl; ++il)
      d[ii][il] = e[ii][il];
    ii+=is;
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      int iemax = g[ii];
      int iemin = g[ji];
      int dg = abs(iemax-iemin); // sparse grid delta
      int kmin = (int)ceil( rmin*dg);
      int kmax = (int)floor(rmax*dg);
      for (int il=0; il<nl; ++il) {
        float dm = dmax;
        int mi = 0;
        for (int k=kmin; k<=kmax; k++) {
          int rk = k*ic;
          int ik = il+rk;
          if (ik<0 || ik>nlm1)
            continue;
          float dc = d[ji][ik];
          if (dc<dm) {
            dm = dc;
            mi = rk;
          }
        }
        m[ji][il] = mi;
        d[ii][il] = dm+e[ii][il];
      }
    }
  }
  
  private static void backtrack(
      int dir, Sampling shifts, float[][] d, float[][] m, float[] u) 
    {
      int nl = d[0].length;
      int ni = d.length;
      int nlm1 = nl-1;
      int nim1 = ni-1;
      int ib = (dir>0)?0:nim1; // begining index
      int ie = (dir>0)?nim1:0; // end index
      int is = (dir>0)?1:-1;   // stride
      int ii = ib;
      // Set initial lag for the case that all errors at ii are equal.
      int il = (dir>0)?0:nlm1; 
      float dl = d[ii][il]; // Current accumulated error value.
          
      // Find minimum lag value(dl) and index(il) at trace index ii.
      for (int jl=0; jl<nl; ++jl) {
        if (d[ii][jl]<dl) {
          dl = d[ii][jl];
          il = jl;
        }
      }
      u[ii] = (float)shifts.getValue(il);
      while (ii!=ie) {
        ii += is;
        il += (int)m[ii][il];
        u[ii] = (float)shifts.getValue(il);
      }
    }
  
  private static void backtrack(
      int dir, int b1, int bS, Sampling shifts1, Sampling shiftsS,
      float[][][] d, float[][][] e, float[] u1, float[] uS) 
    {
      float ob1 = 1.0f/b1;
      float obS = 1.0f/bS;
      int nl1 = d[0].length;
      int nlS = d[0][0].length;
      int ni = d.length;
      int nl1m1 = nl1-1;
      int nlSm1 = nlS-1;
      int nim1 = ni-1;
      int ib = (dir>0)?0:nim1;
      int ie = (dir>0)?nim1:0;
      int is = (dir>0)?1:-1;
      int ic = (dir>0)?1:-1;
      int ii = ib;
      int il1 = (dir>0)?0:nl1m1;
      int ilS = (dir>0)?0:nlSm1;
      float dl = d[ii][il1][ilS];
      for (int jl1=0; jl1<nl1; ++jl1) {
        for (int jlS=0; jlS<nlS; ++jlS) {
          if (d[ii][jl1][jlS]<dl) {
            dl = d[ii][jl1][jlS];
            il1 = jl1;
            ilS = jlS;
          }
        }
      }
      u1[ii] = (float)shifts1.getValue(il1);
      uS[ii] = (float)shiftsS.getValue(ilS);
      
      // Notes for backtracking:
      // ii, the trace index of the current minimum value
      // ji, the next trace index determined by the direction
      // jb1, the next trace for shift 1 that satisfies strain limits b1 
      // jbS, the next trace for shift S that satisfies strain limits bS
      // il1, the lag index for shift 1 of the current minimum value
      // ilc1, the lag index for shift 1 of a possible min value that
      //       satisfies the constraint ic
      // ilS, the lag index for shift S of the current minimum value
      // ilcS, the lag index for shift S of a possible min value that 
      //       satisifies the constraint ic
      // dc1, dcS, dc1S, di, values on the possible minimum paths 
      while (ii!=ie) {
        int ji  = max(0,min(nim1,ii+is));
        int jb1 = max(0,min(nim1,ii+is*b1));
        int jbS = max(0,min(nim1,ii+is*bS));
        int ilc1 = il1+ic;
        int ilcS = ilS+ic;
        ilc1 = (ilc1==-1)?0:(ilc1==nl1)?nl1m1:ilc1;
        ilcS = (ilcS==-1)?0:(ilcS==nlS)?nlSm1:ilcS;
        float dc1 = d[jb1][ilc1][ilS ];
        float dcS = d[jbS][il1 ][ilcS];
        float dc1S= d[jbS][ilc1][ilcS];
        float di  = d[ji ][il1 ][ilS ];
        for (int kb1=ji; kb1!=jb1; kb1+=is) {
          dc1 += e[kb1][ilc1][ilS ];
        }
        for (int kbS=ji; kbS!=jbS; kbS+=is) {
          dcS += e[kbS][il1 ][ilcS];
          dc1S+= e[kbS][ilc1][ilcS];
        }
        dl = min4(dc1,dcS,dc1S,di);
        if (dl!=di) {
          if (dl==dc1S) {
            il1 = ilc1;
            ilS = ilcS;
          } else if (dl==dc1) {
            il1 = ilc1;
          } else if (dl==dcS) {
            ilS = ilcS;
          }
        }
        ii += is;
        u1[ii] = (float)shifts1.getValue(il1);
        uS[ii] = (float)shiftsS.getValue(ilS);
        
        // Adjust shifts for strain limits, if applicable.
        float du1 = 0.0f;
        float duS = 0.0f;
        if (il1==ilc1 && ilS==ilcS) {
          du1 = (u1[ii]-u1[ii-is])*obS; // u1 must satisfy the strain
          duS = (uS[ii]-uS[ii-is])*obS; // limits imposed by uS.
          u1[ii] = u1[ii-is]+du1;
          uS[ii] = uS[ii-is]+duS;
          for (int kbS=ji; kbS!=jbS; kbS+=is) {
            ii += is;
            u1[ii] = u1[ii-is]+du1;
            uS[ii] = uS[ii-is]+duS;
          }
        } else if (il1==ilc1) {
          du1 = (u1[ii]-u1[ii-is])*ob1;
          u1[ii] = u1[ii-is]+du1;
          for (int kb1=ji; kb1!=jb1; kb1+=is) {
            ii += is;
            u1[ii] = u1[ii-is]+du1;
            uS[ii] = uS[ii-is];
          }
        } else if (ilS==ilcS) {
          duS = (uS[ii]-uS[ii-is])*obS;
          uS[ii] = uS[ii-is]+duS;
          for (int kbS=ji; kbS!=jbS; kbS+=is) {
            ii += is;
            u1[ii] = u1[ii-is];
            uS[ii] = uS[ii-is]+duS;
          }
        }
      }
    }
  
  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    int nl = e[0].length;
    int n1 = e.length;
    float eshift = emin;
    float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        e[i1][il] = (e[i1][il]-eshift)*escale;
      }
    }
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] ef = e;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          ef[i2][i1][il] = (ef[i2][i1][il]-eshift)*escale;
        }
      }
    }});
  }
  
  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][][] ef = e;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          for (int il=0; il<nl; ++il) {
            ef[i3][i2][i1][il] = (ef[i3][i2][i1][il]-eshift)*escale;
          }
        }  
      }
    }});
  }

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }
  
  private static float min4(float a, float b, float c, float d) {
    float min = a;
    if (b<=min)
      min = b;
    if (c<=min)
      min = c;
    if (d<=min)
      min = d;
    return min;
  }

  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }

  private static void print(String s) {
    System.out.println(s);
  }
}
