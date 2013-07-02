package warp;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.dsp.SincInterp.Extrapolation;
import edu.mines.jtk.interp.BicubicInterpolator2;
import edu.mines.jtk.interp.BilinearInterpolator2;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Smooth dynamic warping for converted wave images.
 */
public class DynamicWarpingC {

  /**
   * Constants for interpolation.
   */
  public enum Interp {
    LINEAR,
    MONOTONIC,
    SPLINE
  }

  /**
   * Constructor for smooth dynamic warping of converted wave images. Required
   * parameters are minimum and maximum lag. The number of lags is nl = 
   * {@code lMax}-{@code lMin}+1. For warping PS to PP images the typical 
   * minimum lag is 0, but negative values may need to be allowed due to statics
   * corrections applied to the data. The maximum lag can be computed from the
   * {@link #computeMaxLag(int, float)} method.
   * @param lMin the minimum lag.
   * @param lMax the maximum lag.
   */
  public DynamicWarpingC(int lMin, int lMax) {
    Check.argument(lMin<=lMax,"lMin<=lMax");
    _nl = lMax-lMin+1;
    _sl1 = new Sampling(_nl,1.0,lMin);
    _r1Min =  0.0;
    _r1Max =  1.0;
    _r2Min = -1.0;
    _r2Max =  1.0;
    _r3Min = -1.0;
    _r3Max =  1.0;
    _si = new SincInterp();
    _si.setExtrapolation(Extrapolation.CONSTANT);
  }
  
  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   */
  public void setStrainLimits(double r1Min, double r1Max) {
    setStrainLimits(r1Min,r1Max,-1.0,1.0,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   * @param r2Min lower bound on strain in 2nd dimension.
   * @param r2Max upper bound on strain in 2nd dimension.
   */
  public void setStrainLimits(
    double r1Min, double r1Max,
    double r2Min, double r2Max)
  {
    setStrainLimits(r1Min,r1Max,r2Min,r2Max,-1.0,1.0);
  }

  /**
   * Sets bounds on strains for this dynamic warping.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1Min lower bound on strain in 1st dimension.
   * @param r1Max upper bound on strain in 1st dimension.
   * @param r2Min lower bound on strain in 2nd dimension.
   * @param r2Max upper bound on strain in 2nd dimension.
   * @param r3Min lower bound on strain in 3rd dimension.
   * @param r3Max upper bound on strain in 3rd dimension.
   */
  public void setStrainLimits(
    double r1Min, double r1Max,
    double r2Min, double r2Max,
    double r3Min, double r3Max)
  {
    _r1Min = r1Min; _r1Max = r1Max;
    _r2Min = r2Min; _r2Max = r2Max;
    _r3Min = r3Min; _r3Max = r3Max;
  }
  
  /**
   * Returns the number of lags.
   * @return the number of lags.
   */
  public int getNumberOfLags() {
    return _nl;
  }
  
  /**
   * Find shifts for 1D traces. For the input trace {@code f[n1]}, shifts 
   * {@code u} are computed on a subsampled grid such that the computed shifts 
   * are {@code u[ng1]}. This length matches the length of array {@code g1} 
   * and the indices of the subsampled shifts are specified by contents of 
   * this array. The contents of the {@code g1} array are rounded to the 
   * nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input trace, that is 
   * {@code ui[n1]}.
   * @param f the PP trace.
   * @param g the PS trace.
   * @param g1 array of size [ng1] specifying first dimension sparse grid
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @return shifts for 1D traces.
   */
  public float[] findShifts(float[] f, float[] g, float[] g1, Interp interp1) {
    int n1 = f.length;
    int ng1 = g1.length;
    final int[] g1i = new int[ng1];
    for (int ig1=0; ig1<ng1; ig1++) {
      g1i[ig1] = (int)(g1[ig1]+0.5f);
    }
    float[][] e = computeErrors(f,g);
    float[][][] dm = accumulateForwardSparse(e,_r1Min,_r1Max,g1i);
    float[] u = backtrackReverse(dm[0],dm[1]);
    return interpolate(n1,g1,u,interp1);
  }

  /**
   * Find shifts for 2D images. For the input image {@code f[n2][n1]}, shifts 
   * {@code u} are computed on a subsampled grid such that the computed shifts 
   * are {@code u[ng2][ng1]}. These lengths match the length of arrays 
   * {@code g1,g2} and the indices of the subsampled shifts are specified by 
   * contents of these arrays. The contents of the {@code g1} array are rounded
   * to the nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is 
   * {@code ui[n2][n1]}.
   * @param f the PP traces.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid 
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp2 interpolation method for the i2 (slow) dimension.
   * @return shifts for 2D images.
   */
  public float[][] findShifts(
      float[][] f, float[][] g, final float[][] g1, final int[] g2,
      Interp interp1, Interp interp2)
  {
    int n2 = f.length;
    int n1 = f[0].length;
    Check.argument(n2==g1.length,"_n2==g1.length");
    int ng1 = g1[0].length;
    final int[][] g1i = new int[n2][ng1];
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<ng1; i1++) {
        g1i[i2][i1] = (int)(g1[i2][i1]+0.5f);
      }
    }
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][] es1 = smoothErrors1(f,g,_r1Min,_r1Max,g1i);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);
    
    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][] es = smoothErrors2(es1,_r2Min,_r2Max,g2);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es);
    
    final int ng2 = es.length;
    final float[][] u = new float[ng2][];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(es[i2],g1i[g2[i2]],_r1Min,_r1Max);
      u[i2] = backtrackReverse(dm[0],dm[1]);
    }});
    for (int i2=0; i2<ng2; i2++) {
      for (int i1=1; i1<g1[0].length; i1++) {
        float n = u[i2][i1] - u[i2][i1-1];
        float d = g1i[g2[i2]][i1] - g1i[g2[i2]][i1-1];
        float r = n/d;
        assert r>=_r1Min && r<=_r1Max:"n="+n+", d="+d+", r="+r;
      }
    }
    return interpolate(n1,n2,g1,g2,u,interp1,interp2);
  }

  /**
   * Find shifts for 3D images. For the input image {@code f[n3][n2][n1]},
   * shifts {@code u} are computed on a subsampled grid such that the computed
   * shifts are {@code u[ng3][ng2][ng1]}. These lengths match the length of
   * arrays {@code g1,g2,g3} and the indices of the subsampled shifts are 
   * specified by contents of these arrays. The contents of the {@code g1} 
   * array are rounded to the nearest integer.
   * </p>
   * The sparsely computed shifts are interpolated back to the fine grid such
   * that the returned shifts match the size of the input image, that is 
   * {@code ui[n3][n2][ne1]}.
   * @param f the PP traces.
   * @param g the PS traces.
   * @param g1 array of size [n2][ng1] specifying first dimension sparse grid
   *  locations for all n2.
   * @param g2 array of size [ng2] specifying second dimension sparse grid 
   *  locations.
   * @param g3 array of size [ng3] specifying third dimension sparse grid 
   *  locations.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp23 interpolation method for the i2 (middle) and i3
   *  (slow) dimensions.
   * @return shifts for 3D images.
   */
  public float[][][] findShifts(
      float[][][] f, float[][][] g, 
      final float[][][] g1, final int[] g2, final int[] g3,
      Interp interp1, Interp interp23)
  {
    int n3 = f.length;
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    Check.argument(n3==g1.length,"_n3==g1.length");
    Check.argument(n2==g1[0].length,"_n2==g1[0].length");
    int ng1 = g1[0][0].length;
    final int[][][] g1i = new int[n3][n2][ng1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<ng1; i1++) {
          g1i[i3][i2][i1] = (int)(g1[i3][i2][i1]+0.5f);
        }
      }
    }
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][][] es1 = smoothErrors1(f,g,_r1Min,_r1Max,g1i);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);

    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][][] es12 = smoothErrors2(es1,_r2Min,_r2Max,g2);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es12);
    
    s.restart();
    print("Smoothing 3rd dimension...");
    final float[][][][] es = smoothErrors3(es12,_r3Min,_r3Max,g3);
    normalizeErrors(es);
    print("Finished 3rd dimension smoothing in "+s.time()+" seconds");
    
    final int ng2 = es[0].length;
    final int ng3 = es.length;
    final float[][][] u = new float[ng3][ng2][];
    Parallel.loop(ng3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<ng2; i2++) {
        float[][][] dm = 
            accumulateForward(es[i3][i2],g1i[g3[i3]][g2[i2]],_r1Min,_r1Max);
        u[i3][i2] = backtrackReverse(dm[0],dm[1]);  
      }
    }});
    for (int i3=0; i3<ng3; i3++) {
      for (int i2=0; i2<ng2; i2++) {
        for (int i1=1; i1<ng1; i1++) {
          float n = u[i3][i2][i1] - u[i3][i2][i1-1];
          float d = g1i[g3[i3]][g2[i2]][i1] - g1i[g3[i3]][g2[i2]][i1-1];
          float r = n/d;
          assert r>=_r1Min && r<=_r1Max:"n="+n+", d="+d+", r="+r;
        }
      }
    }
    return interpolate(n1,n2,n3,g1,g2,g3,u,interp1,interp23);
  }

  /**
   * Applies the shifts {@code u} to the PS trace {@code g}.
   * @param n1f the length of PP trace. This is the length of the returned 
   *  warped trace.
   * @param g the PS trace to be warped.
   * @param u the shifts that warp the PS trace to the PP trace.
   * @return the warped PS trace.
   */
  public float[] applyShifts(int n1f, float[] g, float[] u) {
    int n1g = g.length;
    int nu = u.length;
    int num = nu-1;
    float[] h = new float[n1f];
    for (int iu=0; iu<nu; ++iu) {
      h[iu] = _si.interpolate(n1g,1.0,0.0,g,iu+u[iu]);
    }
    for (int i1=nu; i1<n1f; ++i1) {
      h[i1] = _si.interpolate(n1g,1.0,0.0,g,i1+u[num]);
    }
    return h;
  }
  
  /**
   * Applies the shifts {@code u} to the PS image {@code g}.
   * @param n1f the length of PP traces. This is the length of the returned 
   *  warped traces.
   * @param g the PS trace to be warped.
   * @param u the shifts that warp the PS image to the PP image.
   * @return the warped PS image.
   */
  public float[][] applyShifts(
      final int n1f, final float[][] g, final float[][] u) 
  {
    final int n2 = g.length;
    final int n1g = g[0].length;
    final int n1u = u[0].length;
    final int n1um = n1u-1;
    final float[][] hf = new float[n2][n1f];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1u; ++i1) {
        hf[i2][i1] = _si.interpolate(n1g,1.0,0.0,g[i2],i1+u[i2][i1]);
      }
      for (int i1=n1u; i1<n1f; ++i1) {
        hf[i2][i1] = _si.interpolate(n1g,1.0,0.0,g[i2],i1+u[i2][n1um]);
      }
    }});
    return hf;
  }
  
  /**
   * Applies the shifts {@code u} to the PS image {@code g}.
   * @param n1f the length of PP traces. This is the length of the returned 
   *  warped traces.
   * @param g the PS trace to be warped.
   * @param u the shifts that warp the PS image to the PP image.
   * @return the warped PS image.
   */
  public float[][][] applyShifts(
      final int n1f, final float[][][] g, final float[][][] u) 
  {
    final int n3 = g.length;
    final int n2 = g[0].length;
    final int n1g = g[0][0].length;
    final int n1u = u[0][0].length;
    final int n1um = n1u-1;
    final float[][][] hf = new float[n3][n2][n1f];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1u; i1++) {
          hf[i3][i2][i1] = _si.interpolate(n1g,1.0,0.0,g[i3][i2],
              i1+u[i3][i2][i1]);
        }
        for (int i1=n1u; i1<n1f; i1++) {
          hf[i3][i2][i1] = _si.interpolate(n1g,1.0,0.0,g[i3][i2],
              i1+u[i3][i2][n1um]);
        }  
      }
    }});
    return hf;
  }

  /**
   * Compute the maximum length of PP traces from a guess of the average Vp/Vs
   * ratio. This method can be used to find the maximum PP length useful for 
   * warping. Truncating the PP data to this length will improve efficiency.
   * This calculation uses the entire length of the PS trace {@code n1PS} and 
   * determines the maximum sample in the PP trace that could correspond to the 
   * final PS sample. The length is [2.0/({@code vpvsAvg}+1.0)]*{@code n1PS} or 
   * {@code n1PP} if the computed value is larger than {@code n1PP}.
   * @param n1PP the length of PP traces.
   * @param n1PS the length of PS traces.
   * @param vpvsAvg a guess of the average Vp/Vs ratio at the PP sample that 
   *  corresponds with the maximum PS sample. Usually a value of 2.0 is a good 
   *  starting point.
   * @return the maximum length of PP traces useful for warping, based on the 
   *  {@code vpvsAvg} and {@code n1PS}.
   */
  public static int computeMaxLength(int n1PP, int n1PS, float vpvsAvg) {
    float scale = getScale(vpvsAvg);
    int n1PPMax = (int)(n1PS*scale);
    return (n1PPMax>n1PP)?n1PP:n1PPMax;
  }

  /**
   * Computes the maximum lag for warping PS to PP traces from a guess of the 
   * average Vp/Vs ratio. This value can be used to construct a 
   * {@link DynamicWarpingC} instance. This method uses the approach described 
   * in the {@link #computeMaxLength(int, int, float)} method. The maximum lag 
   * is simply the difference between the computed PP trace length useful for 
   * warping and the PS trace length.
   * @param n1PS the length of PS traces.
   * @param vpvsAvg a guess of the average Vp/Vs ratio at the PP sample that 
   *  corresponds with the maximum PS sample. Usually a value of 2.0 is a good 
   *  starting point.
   * @return the maximum lag (in samples) for warping, based on the 
   *  {@code vpvsAvg} and {@code n1PS}. 
   */
  public static int computeMaxLag(int n1PS, float vpvsAvg) {
    float scale = getScale(vpvsAvg);
    int n1PPMax = (int)(n1PS*scale);
    int lMax = n1PS-n1PPMax-1;
    Check.argument(lMax>0,"lMax>0");
    return lMax;
  }

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace 
   * {@code f}. The NRMS value is the RMS of the difference between {@code f} 
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this 
   *  computation.
   * @param f the PP trace.
   * @param h the warped PS trace.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(int n1Max, float[] f, float[] h) {
    int n1f = f.length;
    int n1h = h.length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    float[] fs = copy(n1Max,f);
    float[] hs = copy(n1Max,h);
    float scale = 1.0f/n1Max;
    float[] d = sub(hs,fs);
    float frms = sqrt(sum(mul(fs,fs))*scale);  
    float hrms = sqrt(sum(mul(hs,hs))*scale);
    float drms = sqrt(sum(mul( d, d))*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
  }

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace 
   * {@code f}. The NRMS value is the RMS of the difference between {@code f} 
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this 
   *  computation.
   * @param f the PP traces.
   * @param h the warped PS traces.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(
      final int n1Max, final float[][] f, final float[][] h) 
  {
    int n1f = f[0].length;
    int n1h = h[0].length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    int n2 = f.length;
    float scale = 1.0f/(n1Max*n2);
    float[] rms = Parallel.reduce(n2,new Parallel.ReduceInt<float[]>() {
      @Override
      public float[] compute(int i2) {
        float[] fhdSq = new float[3];
        for (int i1=0; i1<n1Max; i1++) {
          float fv = f[i2][i1];
          float hv = h[i2][i1];
          fhdSq[0] += fv*fv;
          fhdSq[1] += hv*hv;
          fhdSq[2] += (hv-fv)*(hv-fv);
        }
        return fhdSq;
      }

      @Override
      public float[] combine(float[] fhdSq1, float[] fhdSq2) {
        float[] rms = new float[3];
        rms[0] = fhdSq1[0]+fhdSq2[0];
        rms[1] = fhdSq1[1]+fhdSq2[1];
        rms[2] = fhdSq1[2]+fhdSq2[2];
        return rms;
      }
    });
    float frms = sqrt(rms[0]*scale);
    float hrms = sqrt(rms[1]*scale);
    float drms = sqrt(rms[2]*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
  }

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace 
   * {@code f}. The NRMS value is the RMS of the difference between {@code f} 
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this 
   *  computation.
   * @param f the PP traces.
   * @param h the warped PS traces.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(
      final int n1Max, final float[][][] f, final float[][][] h) 
  {
    int n1f = f[0][0].length;
    int n1h = h[0][0].length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n23 = n2*n3;
    float scale = 1.0f/(n1Max*n23);
    float[] rms = Parallel.reduce(n23,new Parallel.ReduceInt<float[]>() {
      @Override
      public float[] compute(int i23) {
        int i2 = i23%n2;
        int i3 = i23/n2;
        float[] fhdSq = new float[3];
        for (int i1=0; i1<n1Max; i1++) {
          float fv = f[i3][i2][i1];
          float hv = h[i3][i2][i1];
          fhdSq[0] += fv*fv;
          fhdSq[1] += hv*hv;
          fhdSq[2] += (hv-fv)*(hv-fv);
        }
        return fhdSq;
      }

      @Override
      public float[] combine(float[] fhdSq1, float[] fhdSq2) {
        float[] rms = new float[3];
        rms[0] = fhdSq1[0]+fhdSq2[0];
        rms[1] = fhdSq1[1]+fhdSq2[1];
        rms[2] = fhdSq1[2]+fhdSq2[2];
        return rms;
      }
    });
    float frms = sqrt(rms[0]*scale);
    float hrms = sqrt(rms[1]*scale);
    float drms = sqrt(rms[2]*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. The relationship is 
   * defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values. 
   */
  public static float[] vpvs(float[] u) {
    int n = u.length;
    int nm1 = n-1;
    float[] vpvs = new float[n];
    vpvs[ 0 ] = 1.0f + 2.0f*(u[ 1 ]-u[  0  ]); // at i1=0, forward difference
    vpvs[nm1] = 1.0f + 2.0f*(u[nm1]-u[nm1-1]); // at i1=nm1, backward difference
    for (int i1=1; i1<nm1; i1++)
      vpvs[i1] = 1.0f + (u[i1+1]-u[i1-1]);
    return vpvs;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. The relationship is 
   * defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values. 
   */
  public static float[][] vpvs(float[][] u) {
    int n2 = u.length;
    float[][] vpvs = new float[n2][];
    for (int i2=0; i2<n2; i2++)
      vpvs[i2] = vpvs(u[i2]);
    return vpvs;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. The relationship is 
   * defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values. 
   */
  public static float[][][] vpvs(float[][][] u) {
    int n3 = u.length;
    float[][][] vpvs = new float[n3][][];
    for (int i3=0; i3<n3; i3++)
      vpvs[i3] = vpvs(u[i3]);
    return vpvs;
  }

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[] extrapolate(int n1, float[] f) {
    int n1f = f.length;
    float v = f[n1f-1];
    float[] ef = new float[n1];
    copy(n1f,f,ef);
    for (int i1=n1f; i1<n1; i1++)
      ef[i1] = v;
    return ef;
  }

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[][] extrapolate(int n1, float[][] f) {
    int n2 = f.length;
    float[][] ef = new float[n2][];
    for (int i2=0; i2<n2; i2++)
      ef[i2] = extrapolate(n1,f[i2]);
    return ef;
  }

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[][][] extrapolate(int n1, float[][][] f) {
    int n3 = f.length;
    float[][][] ef = new float[n3][][];
    for (int i3=0; i3<n3; i3++)
      ef[i3] = extrapolate(n1,f[i3]);
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
  public float[] findShifts2(float[][] f, float[][] g, int[] g1) {
    final float[][] e = computeErrorsSum2(f,g);
    final float[][][] dm = accumulateForwardSparse(e,_r1Min,_r1Max,g1);
    return backtrackReverse(dm[0],dm[1]);
  }

  public float[] findShifts3(float[][][] f, float[][][] g, int[] g1) {
    float[][] e3Avg = computeErrorsSum3(f,g);
    final float[][][] dm = accumulateForwardSparse(e3Avg,_r1Min,_r1Max,g1);
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
  public float[][] computeErrors(float[] f, float[] g) {
    int n1 = f.length;
    float[][] e = new float[n1][_nl];
    computeErrors(f,g,e);
    normalizeErrors(e);
    return e;
  }

  /**
   * Compute alignment errors for 2D traces.
   * @return alignment errors for 2D traces.
   */
  public float[][][] computeErrors2(final float[][] f, final float[][] g) {
    int n2 = f.length;
    int n1 = f[0].length;
    final float[][][] e = new float[n2][n1][_nl];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrors(f[i2],g[i2],e[i2]);  
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
  public float[][][] computeErrors2Near(
      final float[][] f, final float[][] g, final int w) 
  {
    int n2 = f.length;
    int n1 = f[0].length;
    final float[][][] e = new float[n2][n1][_nl];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      e[i2] = computeErrorsNear(f,g,i2,w);
    }});
    return e;
  }

  /**
   * Compute summed alignment errors for 2D images.
   * @return summed alignment errors for 2D images. 
   */
  public float[][] computeErrorsSum2(final float[][] f, final float[][] g) {
    final int n2 = f.length;
    final int n1 = f[0].length;
    float[][] e = Parallel.reduce(n2,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i2) {
      float[][] e = new float[n1][_nl];
      computeErrors(f[i2],g[i2],e);
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
  public float[][] computeErrorsSum3(
      final float[][][] f, final float[][][] g) 
  {
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n1 = f[0][0].length;
    float[][] e = Parallel.reduce(n2*n3,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i23) {
      int i2 = i23%n2;
      int i3 = i23/n2;
      float[][] e = new float[n1][_nl];
      computeErrors(f[i3][i2],g[i3][i2],e);
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    return e;
  }

  public float[][][] accumulateForward(
      float[][] e, int[] g, double rMin, double rMax)
  {
    float[][] d = new float[e.length][e[0].length];
    float[][] m = new float[e.length][e[0].length];
    accumulateFromSparse(1,e,d,m,g,rMin,rMax);
    return new float[][][]{d,m};
  }

  public static float[][][] accumulateForwardSparse(
      float[][] e, double rmin, double rmax, int[] g)
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
    backtrack(-1,_sl1,d,m,u);
    return u;
  }
  
  public float[] backtrackForward(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrack(1,_sl1,d,m,u);
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

  private double _r1Min;
  private double _r1Max;
  private double _r2Min;
  private double _r2Max;
  private double _r3Min;
  private double _r3Max;
  private int _nl; // number of lags
  private Sampling _sl1; // sampling of pp to ps1 lags
  private Sampling _shiftsS; // sampling of ps1 to ps2 shift values
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private SincInterp _si; // for warping with non-integer shifts
  
  /**
   * Computes scale to apply to PS traces. n1PP = scale*n1PS
   * @param vpvsAvg
   * @return a scaler to apply to PS traces.
   */
  private static float getScale(float vpvsAvg) {
    Check.argument(vpvsAvg>=1,"vpvsAvg>=1");
    return 2.0f/(vpvsAvg+1.0f);
  }
  
  
  /**
   * Interpolates subsampled shifts u[ng1] to uniformly sampled shifts
   * ui[_ne1].
   * @param g1 sparse grid indices.
   * @param u sparse shifts.
   * @param doLinear {@code true} for linear interpolation,
   *  {@code false} for cubic.
   * @return the interpolated shifts.
   */
  private float[] interpolate(int n1, float[] g1, float[] u, Interp interp1) {
    int ng1 = g1.length;
    float[] ui = new float[n1];
    CubicInterpolator.Method m1;
    switch (interp1) {
      case LINEAR:    m1 = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m1 = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m1 = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
          interp1.toString()+" is not a recognized interpolation method.");
    }
    CubicInterpolator ci = new CubicInterpolator(m1,ng1,g1,u);
    for (int i=0; i<n1; i++)
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
      int n1, int n2, float[][] g1, int[] g2, float[][] u,
      Interp interp1, Interp interp2)
  {
    int ng1 = g1[0].length;
    int ng2 = g2.length;
    float[] g2f = new float[ng2];
    for (int ig2=0; ig2<ng2; ig2++)
      g2f[ig2] = g2[ig2];
    
    float[][] ui = new float[n2][n1];
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
    float[][] ui2 = new float[n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i2=0; i2<ng2; i2++)
        u2[i2] = u[i2][i1];
      CubicInterpolator ciu = new CubicInterpolator(m2,g2f,u2);
      for (int i2=0; i2<n2; i2++)
        ui2[i2][i1] = ciu.interpolate(i2);
    }

    // interpolate in the first dimension.
    for (int i2=0; i2<n2; i2++) {
      CubicInterpolator ci = new CubicInterpolator(m1,g1[i2],ui2[i2]);
      for (int i1=0; i1<n1; i1++)
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
      int n1, int n2, int n3, float[][][] g1, int[] g2, int[] g3, float[][][] u, 
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
    float[][][] ui23 = new float[n3][n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i3=0; i3<ng3; i3++)
        for (int i2=0; i2<ng2; i2++)
          u23[i3][i2] = u[i3][i2][i1];
      if (doLinear) {
        BilinearInterpolator2 bli = new BilinearInterpolator2(g2f,g3f,u23);
        for (int i3=0; i3<n3; i3++)
          for (int i2=0; i2<n2; i2++)
            ui23[i3][i2][i1] = bli.interpolate(i2,i3);
      } else {
        BicubicInterpolator2 bci = 
            new BicubicInterpolator2(m23,m23,g2f,g3f,u23);
        for (int i3=0; i3<n3; i3++)
          for (int i2=0; i2<n2; i2++)
            ui23[i3][i2][i1] = bci.interpolate(i2,i3);
      }
    }

    float[][][] ui = new float[n3][n2][n1];
    // interpolate the first dimension
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        CubicInterpolator ci = 
            new CubicInterpolator(m1,g1[i3][i2],ui23[i3][i2]);
        for (int i1=0; i1<n1; i1++)
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

//  private void initShifts(int n1pp, int n1ps, float vpvsAvg, int sMin) {
//    _sMin = sMin;
//    float scale = 2.0f/(vpvsAvg+1.0f);
//    int[] el = computeErrorLengths(n1pp,n1ps,0,scale);
//    _nel = el[1]-_sMin;
//    _ne1 = el[0];
//    _n1pp = n1pp;
//    _n1ps = n1ps;
//    _sl1 = new Sampling(_nel,1.0,_sMin);
//  }
  
  /**
   * Computes alignment errors for {@code f} and {@code g}. Note
   * that values of the {@code e} array are not replaced, but
   * added to. 
   * @param f
   * @param g
   * @param e
   */
  private void computeErrors(float[] f, float[] g, float[][] e) {
    int n1Max = e.length;
    int nl = e[0].length;
    int ng = g.length;
    float[] gi = new float[n1Max];
    for (int il=0; il<nl; il++) {
      _si.interpolate(ng,1.0,0.0,g,n1Max,1.0,_sl1.getValue(il),gi);
      for (int i1=0; i1<n1Max; i1++) {
        e[i1][il] += error(f[i1],gi[i1]);
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
  private float[][] computeErrorsNear(float[][] f, float[][] g, int i2, int w) {
    int n2 = f.length;
    int n1 = f[0].length;
    final float[][] e = new float[n1][_nl];
    final int n2m1 = n2-1;
    int i2min = max(0,i2-w);
    int i2max = min(n2m1,i2+w);
    for (int j=i2min; j<=i2max; j++)
      computeErrors(f[j],g[j],e);
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
        for (int ilS=0,jS=j1; ilS<nlS; ++ilS,++jS) {
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
      float[][] e, double rmin, double rmax, int[] g)      
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
      final double r1min, final double r1max, final int[][] g1)
  {
    final int n2 = pp.length;
    final int n1 = pp[0].length;
    final int ng1 = g1.length;
    final float[][][] es1 = new float[n2][ng1][_nl];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[n1][_nl];
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
      final double r1min, final double r1max, final int[][][] g1)
  {
    final int n3 = pp.length;
    final int n2 = pp[0].length;
    final int n1 = pp[0][0].length;
    final int ng1 = g1.length;
    final float[][][][] es1 = new float[n3][n2][ng1][_nl];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        float[][] e = new float[n1][_nl];
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
      final double r2min, final double r2max, final int[] g2)
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
   * @param r2Min minimum slope in the second dimension.
   * @param r2Max maximum slope in the second dimension.
   * @param g2 second dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [e.length][g2.length][e[0][0].length][e[0][0][0].length].
   */
  private static float[][][][] smoothErrors2(
      final float[][][][] e, 
      final double r2Min, final double r2Max, final int[] g2)
  {
    final int n3 = e.length;
    final float[][][][] es = new float[n3][][][]; // smoothed errors
    for (int i3=0; i3<n3; i3++)
      es[i3] = smoothErrors2(e[i3],r2Min,r2Max,g2);
    return es;
  }

  /**
   * Returns alignment errors smoothed in the third dimension.
   * Returned errors are sparse in the third dimension, and
   * unchanged in the first and second dimension.
   * @param e alignment errors.
   * @param r3Min minimum slope in the third dimension.
   * @param r3Max maximum slope in the third dimension.
   * @param g3 third dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [g3.length][e[0].length][e[0][0].length][e[0][0][0].length].
   */
  private static float[][][][] smoothErrors3(
      final float[][][][] e, 
      final double r3Min, final double r3Max, final int[] g3)
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
        float[][] es3 = smoothErrors(e3,r3Min,r3Max,g3);
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
      int dir, double rMin, double rMax, int[] g,
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
        kmin = (int) ceil(-rMax*dg);
        kmax = (int)floor(-rMin*dg);
      } else {
        kmin = (int) ceil(-rMin*dg);
        kmax = (int)floor(-rMax*dg);
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
   * @param rMin
   * @param rMax
   */
  private static void accumulateFromSparse(
      int dir, float[][] e, float[][] d, float[][] m, int[] g,
      double rMin, double rMax)
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
      int kmin = (int)ceil( rMin*dg);
      int kmax = (int)floor(rMax*dg);
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
