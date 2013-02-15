package warp;

import java.util.Map;

import util.Viewer;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.HilbertTransformFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.interp.BicubicInterpolator2;
import edu.mines.jtk.interp.BilinearInterpolator2;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.interp.CubicInterpolator.Method;
import edu.mines.jtk.interp.TricubicInterpolator3;
import edu.mines.jtk.interp.TrilinearInterpolator3;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic warping for converted wave images.
 */
public class DynamicWarpingC {

  /**
   * Constructor that specifies the maximum shift (computed
   * from vpvsAvg) and how finely to interpolate when computing
   * errors. Reasonable defaults are vpvsAvg=2.0f, fr=2, frMax=2.
   * @param vpvsAvg
   * @param fr
   * @param frMax
   */
  public DynamicWarpingC(float vpvsAvg, int fr, int frMax) {
    Check.argument(vpvsAvg>=1,"vpvsAvg>=1");
    Check.argument(fr>0,"fr>0");
    _fr = fr;
    _frMax = frMax;
    _scale = 2.0f/(vpvsAvg+1.0f);
    _si = new SincInterp();
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
    _frMax = 1;
    int[] el = computeErrorLengths(n1pp,n1ps,0);
    _nel = el[1];
    _ne1 = el[0];
    _n1pp = n1pp;
    _n1ps = n1ps;
    _pp1 = pp;
    _ps1 = ps;
    System.out.println("PP/PS Traces="+_n2+", PP Sample Length="+_n1pp+
        ", PP Sample Length for Warping="+_ne1+", PS Sample Length="+_n1ps+
        ", Number of Lags="+_nel);
    _si = new SincInterp();
  }
  
  /**
   * Constructor for dynamic warping of 2D PP and PS images. The 
   * PS traces are warped to the PP traces. Note that the number of 
   * PP and PS traces must be the same, but the number of samples 
   * do not need to be the same length. The maximum shift is 
   * determined by the given {@code vpvsAvg}.
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
    _frMax = 1;
    int[] el = computeErrorLengths(n1pp,n1ps,0);
    _nel = el[1];
    _ne1 = el[0];
    _n1pp = n1pp;
    _n1ps = n1ps;
    _n2 = pp.length;
    _pp2 = pp;
    _ps2 = ps;
    System.out.println("PP/PS Traces="+_n2+", PP Sample Length="+_n1pp+
        ", PP Sample Length for Warping="+_ne1+", PS Sample Length="+_n1ps+
        ", Number of Lags="+_nel);
    _si = new SincInterp();
  }
  
  /**
   * Constructor for dynamic warping of 3D PP and PS images. The 
   * PS traces are warped to the PP traces. Note that the number of 
   * PP and PS traces must be the same, but the number of samples 
   * do not need to be the same length. The maximum shift is 
   * determined by the given {@code vpvsAvg}.
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
    _frMax = 1;
    int[] el = computeErrorLengths(n1pp,n1ps,0);
    _nel = el[1];
    _ne1 = el[0];
    _n1pp = n1pp;
    _n1ps = n1ps;
    _n2 = pp[0].length;
    _n3 = pp.length;
    _pp3 = pp;
    _ps3 = ps;
    System.out.println("PP/PS Ensembles="+_n3+", PP/PS Traces="+_n2+
        ", PP Sample Length="+_n1pp+", PP Sample Length for Warping="+_ne1+
        ", PS Sample Length="+_n1ps+", Number of Lags="+_nel);
_si = new SincInterp();
  }
  
  public void setInterpMethod(CubicInterpolator.Method m) {
    _m1 = m;
  }
  public void setInterpMethod(BicubicInterpolator2.Method m) {
    _m2 = m;
  }
  public void setInterpMethod(TricubicInterpolator3.Method m) {
    _m3 = m;
  }
  
  public int getPPErrorLength() {
    return _ne1;
  }
  
  public int getNumberOfLags() {
    return _nel;
  }
  
  /**
   * Find shifts for 1D traces. Shifts are computed on a sparse
   * grid and then interpolated back to the fine grid. The
   * sparse grid has a minimum interval of ceil(1/dr1). This
   * interval is variable in order to include the first and last
   * indices of the fine grid.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param dr1 sparse grid increment in first dimension. 
   * @return shifts for 1D traces.
   */
  public float[] findShifts(float r1min, float r1max, float dr1) {
    float[][] e = computeErrors();
    int dg1 = (int)ceil(1.0f/dr1);
    int[] g1 = Subsample.subsample(_ne1,dg1);
    dump(g1);
    float[][][] dm = accumulateForwardSparse(e,r1min,r1max,g1);
    float[] u = backtrackReverse(dm[0],dm[1]);
    return interpolateSparseShifts(_ne1,g1,u,_m1);
  }
  
  /**
   * Find shifts for 2D images. Shifts are computed on a sparse
   * grid and then interpolated back to the fine grid. The
   * sparse grid has a minimum interval of ceil(1/dr1) in the
   * first dimension and ceil(1/dr2) in the second dimension. 
   * These intervals are variable in order to include the first 
   * and last indices of the fine grid.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param dr1 sparse grid increment in first dimension. 
   * @param r2max maximum slope in second dimension.
   * @param dr2 sparse grid increment in second dimension. 
   * @return shifts for 2D images.
   */
  public float[][] findShifts(
      final float r1min, final float r1max, final float dr1,
      final float r2min, final float r2max, final float dr2)
  {
    final int dg1 = (int)ceil(1.0f/dr1);
    final int dg2 = (int)ceil(1.0f/dr2);
    final int[] g1 = Subsample.subsample(_ne1,dg1);
    final int[] g2 = Subsample.subsample( _n2,dg2);
    dump(g1); dump(g2);
    Stopwatch s = new Stopwatch();
    s.start();
    final float[][][] es1 = smoothErrors1(_pp2,_ps2,r1min,r1max,g1);
    System.out.println("Finished 1st Dimension Smoothing in "+s.time()+
        " seconds");
    normalizeErrors(es1);
    s.restart();
    final float[][][] es = smoothErrors2(es1,r2min,r2max,g2);
    System.out.println("Finished 2nd Dimension Smoothing in "+s.time()+
        " seconds");
    normalizeErrors(es);
    final int ng2 = es.length;
    final float[][] u = new float[ng2][];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(es[i2],g1,r1min,r1max);
      u[i2] = backtrackReverse(dm[0],dm[1]);
    }});
    return interpolateSparseShifts(_ne1,_n2,g1,g2,u);
  }
  
  public float[][][] findShifts(
      final float r1min, final float r1max, final float dr1, 
      final float r2min, final float r2max, final float dr2,
      final float r3min, final float r3max, final float dr3)
  {
    final int dg1 = (int)ceil(1.0f/dr1);
    final int dg2 = (int)ceil(1.0f/dr2);
    final int dg3 = (int)ceil(1.0f/dr3);
    final int[] g1 = Subsample.subsample(_ne1,dg1);
    final int[] g2 = Subsample.subsample( _n2,dg2);
    final int[] g3 = Subsample.subsample( _n3,dg3);
    dump(g1); dump(g2); dump(g3);
//    Stopwatch s = new Stopwatch();
//    s.start();
//    final float[][][][] es1 = smoothErrors1(_pp3,_ps3,r1min,r1max,g1);
//    System.out.println("Finished 1st Dimension Smoothing in "+s.time()+
//        " seconds");
//    normalizeErrors(es1);
//    s.restart();
//    final float[][][] es = smoothErrors2(es1,r2min,r2max,g2);
//    System.out.println("Finished 2nd Dimension Smoothing in "+s.time()+
//        " seconds");
//    normalizeErrors(es);
    final float[][][][] es = 
        smoothErrorsSparse(r1min,r1max,g1,r3min,r2max,g2,r3min,r3max,g3);
    final int ng3 = es.length;
    final int ng2 = es[0].length;
    final float[][][] u = new float[ng3][ng2][];
    Parallel.loop(ng3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<ng2; i2++) {
        float[][][] dm = accumulateForward(es[i3][i2],g1,r1min,r1max);
        u[i3][i2] = backtrackReverse(dm[0],dm[1]);
      }
    }});
    return interpolateSparseShifts(_ne1,_n2,_n3,g1,g2,g3,u);
  }
  
  /**
   * Find shifts for 1D traces. This method first smooths
   * alignment errors. Shifts are computed on a sparse grid
   * and then are interpolated back to the fine grid. The
   * sparse grid has a minimum interval of ceil(1/dr1). This
   * interval is variable in order to include the first and last
   * indices of the fine grid.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param dr1 sparse grid increment in first dimension. 
   * @return shifts for 1D traces.
   */
  public float[] findShiftsSmooth(float r1min, float r1max, float dr1) {
    int dg1 = (int)ceil(1.0f/dr1);
    int[] g1 = Subsample.subsample(_ne1,dg1);
    dump(g1);
    float[][] e = computeErrors();
    float[][] es = smoothErrorsSparse(e,r1min,r1max,g1);
    normalizeErrors(es);
    float[][][] dm = accumulateForward(es,g1,r1min,r1max);
    float[] u = backtrackReverse(dm[0],dm[1]);
    return interpolateSparseShifts(_ne1,g1,u,_m1);
  }
  
  /**
   * Find shifts for 1D traces. Shifts are computed on a sparse
   * grid and then interpolated back to the fine grid. The
   * sparse grid indices are selected preferentially by the 
   * maximum amplitudes in the pp trace, however the indices 
   * must have a minimum interval of ceil(1/dr1). 
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param dr1 sparse grid increment in first dimension. 
   * @return shifts for 1D traces.
   */
  public float[] findShiftsM(
      float r1min, float r1max, float dr1, boolean useNG)
  {
    float[][] e = computeErrors();
    float[] ppe = getEnvelope(_pp1);
    float[] pp = new float[_ne1];
    for (int i1=0; i1<_ne1; i1++)
      pp[i1] = ppe[i1]; // use envelope of trace
    int[] g1;
    int dg1 = (int)ceil(1.0f/dr1);
    if (useNG) {
      int ng = 1+(_ne1-1)/dg1;
      g1 = Subsample.subsample(pp,dg1,ng);
    } else
      g1 = Subsample.subsample(pp,dg1);
    System.out.println("ng="+g1.length);
    dump(g1);
    float[][][] dm = accumulateForwardSparse(e,r1min,r1max,g1);
    float[] u = backtrackReverse(dm[0],dm[1]);
    return interpolateSparseShifts(_ne1,g1,u,_m1);
  }
  
  /**
   * Find shifts for 2D images. Shifts are computed on a sparse
   * grid and then interpolated back to the fine grid. For the 
   * first dimension, the sparse grid indices are selected 
   * preferentially by the maximum amplitudes in the pp traces.
   * The sparse grid has a minimum interval of ceil(1/dr2) in the 
   * second dimension. This interval is variable in order to 
   * include the first and last indices of the fine grid.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param dr1 sparse grid increment in first dimension. 
   * @param r2max maximum slope in second dimension.
   * @param dr2 sparse grid increment in second dimension. 
   * @return shifts for 2D images.
   */
  public float[][] findShiftsM(
      float[][] ppf, float[][] x1, 
      final float r1min, final float r1max, final float dr1,
      final float r2min, final float r2max, final float dr2)
  {
    final int[][][] grids = getSparseGrid(ppf,x1,dr1);
    final int[][] g1 = grids[0];
    final int dg2 = (int)ceil(1.0f/dr2);
    final int[] g2 = Subsample.subsample(_n2,dg2);
    dump(g2);
    Stopwatch s = new Stopwatch();
    s.start();
    final float[][][] es1 = smoothErrors1(_pp2,_ps2,r1min,r1max,g1);
    System.out.println("Finished 1st Dimension Smoothing in "+s.time()+
        " seconds");
    normalizeErrors(es1);
    s.restart();
    final float[][][] es = smoothErrors2(es1,r2min,r2max,g2);
    System.out.println("Finished 2nd Dimension Smoothing in "+s.time()+
        " seconds");
    normalizeErrors(es);
//    final float[][][] es = 
//        smoothErrorsSparseM(_pp2,_ps2,r1min,r1max,g1,r2min,r2max,g2);
    // Second smoothing iteration
//    smoothSparseErrors(es,r1min,r1max,g1,g2,r2max);
    
    final int ng2 = es.length;
    final float[][] u = new float[ng2][];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(es[i2],g1[i2],r1min,r1max);
      u[i2] = backtrackReverse(dm[0],dm[1]);
    }});
//    Viewer ve2 = new Viewer(transposeLag12(es),Orientation.X1RIGHT_X2UP);
//    ve2.setTitle("Smoothed Errros 2");
//    ve2.setSize(900,600);
//    ve2.setColorModel1(ColorMap.JET);
//    ve2.setClips1(0.0f,0.9f);
//    ve2.addPoints(u);
//    ve2.show();
//    Viewer vet2 = new Viewer(
//        transposeLag12(transposeLag23(es)),Orientation.X1RIGHT_X2UP);
//    vet2.setTitle("Smoothed Errros 2 Transposed");
//    vet2.setSize(900,600);
//    vet2.setColorModel1(ColorMap.JET);
//    vet2.setClips1(0.0f,0.9f);
//    vet2.addPoints(transposeLag(u));
//    vet2.show();
    return interpolateSparseShifts(_ne1,_n2,grids[1][0],g2,u,x1);
  }
  
  /**
   * Find shifts for 3D images. Shifts are computed on a sparse
   * grid and then interpolated back to the fine grid. For the 
   * first dimension, the sparse grid indices are selected 
   * preferentially by the maximum amplitudes in the pp traces.
   * The sparse grid has a minimum interval of ceil(1/dr2) in the 
   * second dimension. and a minimum interval of ceil(1/dr3) in 
   * the third dimension. These intervals are variable in order to 
   * include the first and last indices of the fine grid.
   * @param ppf flattened pp image for computing sparse grid points.
   * @param u1 
   * @param x1 
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param dr1 sparse grid increment in first dimension. 
   * @param r2max maximum slope in second dimension.
   * @param dr2 sparse grid increment in second dimension. 
   * @param r3max maximum slope in third dimension.
   * @param dr3 sparse grid increment in third dimension.
   * @return shifts for 3D images.
   */
  public float[][][] findShiftsM(
      float[][][] ppf, float[][][] x1, 
      final float r1min, final float r1max, final float dr1,
      final float r2min, final float r2max, final float dr2,
      final float r3min, final float r3max, final float dr3)
  {
    int[][][][] grids = getSparseGrid(ppf,x1,dr1);
    final int dg2 = (int)ceil(1.0f/dr2);
    final int dg3 = (int)ceil(1.0f/dr3);
    final int[][][] g1 = grids[0];
    final int[] g2 = Subsample.subsample(_n2,dg2);
    final int[] g3 = Subsample.subsample(_n3,dg3);
    dump(g2); dump(g3);
    final float[][][][] es = smoothErrorsSparseM(
        r1min,r1max,g1,r2min,r2max,g2,r3min,r3max,g3);
    smoothSparseErrors(es,r1min,r1max,g1,g2,g3,r2max,r3max);
    final int ng2 = es[0].length;
    final int ng3 = es.length;
    final float[][][] u = new float[ng3][ng2][];
    Parallel.loop(ng3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<ng2; i2++) {
        float[][][] dm = accumulateForward(es[i3][i2],g1[i3][i2],r1min,r1max);
        u[i3][i2] = backtrackReverse(dm[0],dm[1]);  
      }
    }});
//    interpolateSparseShiftsFlat(_ne1,_n2,_n3,grids[1][0][0],g2,g3,u);
    return interpolateSparseShifts(_ne1,_n2,_n3,grids[1][0][0],g2,g3,u,x1);
  }
  
  /**
   * Find shifts for 2D images from averaged alignment errors. These
   * shifts are constant for all traces.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param dr1 sparse grid increment in first dimension. 
   * @return shifts for 2D images from averaged alignment errors.
   */
  public float[][] findShifts2(float r1min, float r1max, float dr1) {
    final float[][] e = computeErrors2Average();
    final int dg1 = (int)ceil(1.0f/dr1);
    final int[] g1 = Subsample.subsample(_ne1,dg1);
    dump(g1);
    final float[][][] dm = accumulateForwardSparse(e,r1min,r1max,g1);
    final float[] u1 = interpolateSparseShifts(_ne1,g1,
        backtrackReverse(dm[0],dm[1]),_m1);
    final float[][] u = new float[_n2][];
    for (int i2=0; i2<_n2; i2++) {
      u[i2] = copy(u1);
    }
    return u;
  }
  
  public float[][] findShifts2M(float r1min, float r1max, float dr1) {
    final float[][] e = computeErrors2Average();
    final int dg1 = (int)ceil(1.0f/dr1);
    float[] ea = new float[_ne1];
    float[] et;
    for (int i2=0; i2<_n2; i2++) {
      et = getEnvelope(_pp2[i2]);
      for (int i1=0; i1<_ne1; i1++) {
        ea[i1] += et[i1];
      }
    }
    final int[] g1 = Subsample.subsample(ea,dg1);
    dump(g1);
    final float[][][] dm = accumulateForwardSparse(e,r1min,r1max,g1);
    final float[] u1 = interpolateSparseShifts(_ne1,g1,
        backtrackReverse(dm[0],dm[1]),_m1);
    final float[][] u = new float[_n2][];
    for (int i2=0; i2<_n2; i2++) {
      u[i2] = copy(u1);
    }
    return u;
  }
  
  public float[][] findShifts2MFlat(
      float[][] ppf, float[][] x1, float r1min, float r1max, float dr1) 
  {
    final float[][] e = computeErrors2Average(x1);
    int nppf = ppf[0].length;
    int dg1 = (int)ceil(1.0f/dr1);
    float[] ea = new float[nppf];
    float[] et;
    for (int i2=0; i2<_n2; i2++) {
      et = getEnvelope(ppf[i2]);
      for (int i1=0; i1<nppf; i1++)
        ea[i1] += et[i1];
    }
    int[] g1 = Subsample.subsample(ea,dg1);
    dump(g1);
    final float[][][] dm = accumulateForwardSparse(e,r1min,r1max,g1);
    final float[] u1 = interpolateSparseShifts(_ne1,g1,
        backtrackReverse(dm[0],dm[1]),_m1);
    final float[][] u = new float[_n2][];
    for (int i2=0; i2<_n2; i2++) {
      u[i2] = copy(u1);
    }
    return u;
  }
  
  /**
   * Find shifts for 2D images from averaged alignment errors. These
   * shifts are constant for all traces.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param dr1 sparse grid increment in first dimension. 
   * @return shifts for 2D images from averaged alignment errors.
   */
  public float[][][] findShifts3(float r1min, float r1max, float dr1) {
    float[][] e3Avg = computeErrors3Average();
    return findShifts3(e3Avg,r1min,r1max,dr1);
  }
  
  public float[][][] findShifts3(
      float[][] e3Avg, float r1min, float r1max, float dr1) 
  {
    final int dg1 = (int)ceil(1.0f/dr1);
    final int[] g1 = Subsample.subsample(_ne1,dg1);
    dump(g1);
    final float[][][] dm = accumulateForwardSparse(e3Avg,r1min,r1max,g1);
    final float[] u1 = interpolateSparseShifts(_ne1,g1,
        backtrackReverse(dm[0],dm[1]),_m1);
    final float[][][] u = new float[_n3][_n2][];
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        u[i3][i2] = copy(u1);
      }
    }
    return u;
  }
  
  public float[][][] findShifts3M(float r1min, float r1max, float dr1) {
    final float[][] e3Avg = computeErrors3Average();
    return findShifts3M(e3Avg,r1min,r1max,dr1); 
  }
  
  public float[][][] findShifts3M(
      float[][] e3Avg, float r1min, float r1max, float dr1) 
  {
    final float[][] e = computeErrors3Average();
    final int dg1 = (int)ceil(1.0f/dr1);
    float[] ea = new float[_ne1];
    float[] et;
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        et = getEnvelope(_pp3[i3][i2]);
        for (int i1=0; i1<_ne1; i1++) {
          ea[i1] += et[i1];
        }
      }
    }
    final int[] g1 = Subsample.subsample(ea,dg1);
    dump(g1);
    final float[][][] dm = accumulateForwardSparse(e,r1min,r1max,g1);
    final float[] u1 = interpolateSparseShifts(_ne1,g1,
        backtrackReverse(dm[0],dm[1]),_m1);
    final float[][][] u = new float[_n3][_n2][];
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        u[i3][i2] = copy(u1);
      }
    }
    return u;
  }
  
  public float[][][] findShifts3MFlat(
      float[][][] ppf, float[][][] x1, float r1min, float r1max, float dr1)
  {
    float[][] e3AvgF = computeErrors3Average(x1);
    return findShifts3MFlat(e3AvgF,ppf,x1,r1min,r1max,dr1);
  }
  
  public float[][][] findShifts3MFlat(
      float[][] e3AvgF, float[][][] ppf, float[][][] x1,
      float r1min, float r1max, float dr1) 
  {
    int nppf = ppf[0][0].length;
    int dg1 = (int)ceil(1.0f/dr1);
    float[] ea = new float[nppf];
    float[] et;
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        et = getEnvelope(ppf[i3][i2]);
        for (int i1=0; i1<nppf; i1++)
          ea[i1] += et[i1];
      }
    }
    int[] g1 = Subsample.subsample(ea,dg1);
    dump(g1);
    final float[][][] dm = accumulateForwardSparse(e3AvgF,r1min,r1max,g1);
    final float[] u1 = interpolateSparseShifts(_ne1,g1,
        backtrackReverse(dm[0],dm[1]),_m1);
    final float[][][] u = new float[_n3][_n2][];
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        u[i3][i2] = copy(u1);
      }
    }
    return u;
  }
  
  /**
   * Gets a sparse grid for warping in the first dimension
   * of the pp image. This method flattens the PP image using 
   * given mappings, fm. From the flattened image the average of
   * envelopes of all the traces is computed and used to 
   * preferentially select grid locations with the highest
   * amplitudes, that still honor the minimum interval defined
   * by ceil(1/dr1). These grid location correspond to strong
   * horizons in the flattened space. Again using the provided
   * mappings, we can warp the sparse grid back to the
   * unflattened space. This unflattened sparse grid corresponds
   * to the strong horizons in the original PP image.
   * @param fm Mappings to flatten and unflatten the PP image. 
   *  The validity of the mappings is in no way checked by this
   *   method.
   * @param dr1 sparse grid increment in first dimension. 
   * @return an array [2][_n2][_ne1] of coarse grids for 
   *  warping that follows maximum amplitude events in the PP 
   *  image. The [0] array is the coarse grid in the unflattened
   *  space. The [1] array is the coarse grid in the flat space.
   *  This grid is sparse in the first dimension, but fine in 
   *  the second dimension. 
   */
  public int[][][] getSparseGrid(float[][] ppf, float[][] x1, float dr1) {
    int nppf = ppf[0].length;
    int dg1 = (int)ceil(1.0f/dr1);
    float[] ea = new float[nppf];
    float[] et;
    for (int i2=0; i2<_n2; i2++) {
      et = getEnvelope(ppf[i2]);
      for (int i1=0; i1<nppf; i1++)
        ea[i1] += et[i1];
    }
    int[] g1 = Subsample.subsample(ea,dg1);
    int ng1 = g1.length;
    int ng1m1 = ng1-1; 
    int[][] g = new int[_n2][];
    for (int i2=0; i2<_n2; i2++)
      g[i2] = copy(g1);
    int[][] gw = new int[_n2][ng1];
    for (int i2=0; i2<_n2; i2++) {
      gw[i2][0] = 0;
      gw[i2][ng1m1] = _ne1-1;
      for (int i1=1; i1<ng1m1; i1++) {
        int u1 = g[i2][i1];
        gw[i2][i1] = (int)(x1[i2][u1]);
        // Make sure we don't pass the last grid point _ne1-1
        if (gw[i2][i1]>=(_ne1-1)) {
          int d = ng1m1-i1; // This many left will be out-of-bounds
          for (int id=d; id>0; id--)
            gw[i2][ng1m1-id] = gw[i2][ng1m1]-id;
          for (int i1p=1; i1p<=i1; i1p++) { // Adjust previous grid points
            while (gw[i2][i1-i1p]>=gw[i2][i1-i1p+1])
              gw[i2][i1-i1p] = gw[i2][i1-i1p]-1;
            assert(gw[i2][i1-i1p]<gw[i2][i1-i1p+1]):"i2="+i2+", i1="+i1;
          }
          break; // The rest of i1 indices are squeezed in, next trace.
        }
      }
    }
    return new int[][][]{gw,g};
  }

  /**
   * Gets a sparse grid for warping in the first dimension
   * of the pp image. This method uses the flattened PP image, 
   * ppf, to average the envelopes of all the traces and
   * preferentially select grid locations with the highest
   * amplitudes. The horizon interval still honors the minimum 
   * interval defined by ceil(1/dr1). These grid location 
   * correspond to strong horizons in the flattened space and  
   * are then shifted back to the unflattened space. This 
   * unflattened sparse grid corresponds to the strong horizons
   * in the original PP image.
   * @param ppf
   * @param u1
   * @param x1
   * @param dr1 sparse grid increment in first dimension. 
   * @return an array [2][_n3][_n2][_ne1] of coarse grids for 
   *  warping that follows maximum amplitude events in the PP 
   *  image. The [0] array is the coarse grid in the unflattened
   *  space. The [1] array is the coarse grid in the flat space.
   *  This grid is sparse in the first dimension, but fine in 
   *  the second and third dimensions. 
   */
  public int[][][][] getSparseGrid(
      float[][][] ppf, float[][][] x1, float dr1) 
  {
    int nppf = ppf[0][0].length;
    int dg1 = (int)ceil(1.0f/dr1);
    float[] ea = new float[nppf];
    float[] et;
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        et = getEnvelope(ppf[i3][i2]);
        for (int i1=0; i1<nppf; i1++)
          ea[i1] += et[i1];
      }
    }
    int[] g1 = Subsample.subsample(ea,dg1);
    int ng1 = g1.length;
    int ng1m1 = ng1-1; 
    int[][][] g = new int[_n3][_n2][];
    for (int i3=0; i3<_n3; i3++)
      for (int i2=0; i2<_n2; i2++)
        g[i3][i2] = copy(g1);
    int[][][] gw = new int[_n3][_n2][ng1];
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        gw[i3][i2][0] = 0;
        gw[i3][i2][ng1m1] = _ne1-1;
        for (int i1=1; i1<ng1m1; i1++) {
          int u = g[i3][i2][i1];
          gw[i3][i2][i1] = (int)(x1[i3][i2][u]);
          // Make sure we don't pass the last grid point _ne1-1
          if (gw[i3][i2][i1]>=(_ne1-1)) {
            int d = ng1m1-i1; // This many left will be out-of-bounds
            for (int id=d; id>0; id--)
              gw[i3][i2][ng1m1-id] = gw[i3][i2][ng1m1]-id;
            for (int i1p=1; i1p<=i1; i1p++) { // Adjust previous grid points
              while (gw[i3][i2][i1-i1p]>=gw[i3][i2][i1-i1p+1])
                gw[i3][i2][i1-i1p] = gw[i3][i2][i1-i1p]-1;
              assert(gw[i3][i2][i1-i1p]<gw[i3][i2][i1-i1p+1]):
                "i3="+i3+", i2="+i2+", i1="+i1;
            }
            break; // The rest of i1 indices are squeezed in, next trace.
          }
        }
      }
    }
    return new int[][][][]{gw,g};
  }

  public void plotSparseGrid(
      float[][] ppf, int[][] gf, int[][] gw, float dr2) 
  {
    int[] g2 = Subsample.subsample(_n2,(int)ceil(1.0f/dr2));
    int ng1 = gf[0].length;
    int ng2 = g2.length;
    float[][] x1f = new float[ng2][ng1];
    float[][] x2f = new float[ng2][ng1];
    for (int i2=0; i2<ng2; i2++) {
      for (int i1=0; i1<ng1; i1++) {
        x1f[i2][i1] = (float)gf[g2[i2]][i1];
        x2f[i2][i1] = (float)g2[i2];
      }  
    }
    float[][] x1w = new float[ng2][ng1];
    float[][] x2w = new float[ng2][ng1];
    for (int i2=0; i2<ng2; i2++) {
      for (int i1=0; i1<ng1; i1++) {
        x1w[i2][i1] = (float)gw[g2[i2]][i1];
        x2w[i2][i1] = (float)g2[i2];
      }
    }
    Viewer vf = new Viewer(ppf);
    vf.setTitle("Flattened Image - Coarse Grid (dr2="+dr2+")");
    vf.setHLabel("crossline");
    vf.setVLabel("Tau (samples)");
    vf.addColorBar("Amplitude");
    vf.setClips1(-3.0f,3.0f);
    vf.addPoints2(x1f,x2f);
    vf.setHLimits(0,_n2-1);
    vf.setVLimits(0,ppf[0].length-1);
    vf.setSize(600,900);
    vf.show();
    
    Viewer vw = new Viewer(_pp2);
    vw.setTitle("Original Image - Coarse Grid (dr2="+dr2+")");
    vw.setHLabel("crossline");
    vw.setVLabel("time (samples)");
    vw.addColorBar("Amplitude");
    vw.setClips1(-3.0f,3.0f);
    vw.addPoints2(x1w,x2w);
    vw.setHLimits(0,_n2-1);
    vw.setVLimits(0,_ne1-1);
    vw.setSize(600,900);
    vw.show();
  }
  
  public void plotSparseGrid(
      float[][][] ppf, int[][][] gf, int[][][] gw, float dr2) 
  {
    int[] g2 = Subsample.subsample(_n2,(int)ceil(1.0f/dr2));
    int ng1 = gf[0][0].length;
    int ng2 = g2.length;
    float[][][] x1f = new float[_n3][ng2][ng1];
    float[][][] x2f = new float[_n3][ng2][ng1];
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<ng2; i2++) {
        for (int i1=0; i1<ng1; i1++) {
          x1f[i3][i2][i1] = (float)gf[i3][g2[i2]][i1];
          x2f[i3][i2][i1] = (float)g2[i2];
        }  
      }
    }
    float[][][] x1w = new float[_n3][ng2][ng1];
    float[][][] x2w = new float[_n3][ng2][ng1];
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<ng2; i2++) {
        for (int i1=0; i1<ng1; i1++) {
          x1w[i3][i2][i1] = (float)gw[i3][g2[i2]][i1];
          x2w[i3][i2][i1] = (float)g2[i2];
        }
      }
    }
    Viewer vf = new Viewer(ppf);
    vf.setTitle("Flattened Image - Coarse Grid (dr2="+dr2+")");
    vf.setHLabel("crossline");
    vf.setVLabel("Tau (samples)");
    vf.addColorBar("Amplitude");
    vf.setClips1(-3.0f,3.0f);
    vf.addPoints3(x1f,x2f);
    vf.setHLimits(0,_n2-1);
    vf.setVLimits(0,ppf[0][0].length-1);
    vf.setSize(600,900);
    vf.show();
    
    Viewer vw = new Viewer(_pp3);
    vw.setTitle("Original Image - Coarse Grid (dr2="+dr2+")");
    vw.setHLabel("crossline");
    vw.setVLabel("time (samples)");
    vw.addColorBar("Amplitude");
    vw.setClips1(-3.0f,3.0f);
    vw.addPoints3(x1w,x2w);
    vw.setHLimits(0,_n2-1);
    vw.setVLimits(0,_ne1-1);
    vw.setSize(600,900);
    vw.show();
  }
  
  public float[][] getX1X2(float[] u, float dr1, Sampling s1) {
    s1 = (s1==null)?new Sampling(u.length):s1;
    int dg1 = (int)ceil(1.0f/dr1);
    int[] g1 = Subsample.subsample(_ne1,dg1);
    int ng = g1.length;
    float[] x1 = new float[ng];
    float[] x2 = new float[ng];
    for (int ig=0; ig<ng; ig++) {
      x1[ig] = (float)s1.getValue(g1[ig]);
      x2[ig] = (float)u[g1[ig]];
    }
    return new float[][]{x1,x2};
  }
  
  public float[][] getX1X2M(float[] u, float dr1, Sampling s1, boolean useNG) {
    s1 = (s1==null)?new Sampling(u.length):s1;
    float[] ppe = new float[_ne1];
    for (int i1=0; i1<_ne1; i1++)
      ppe[i1] = _pp1[i1];
    int[] g1;
    int dg1 = (int)ceil(1.0f/dr1);
    if (useNG) {
      int ng = 1+(_ne1-1)/dg1;
      g1 = Subsample.subsample(ppe,dg1,ng);
    } else
      g1 = Subsample.subsample(ppe,dg1);
    int ng = g1.length;
    float[] x1 = new float[ng];
    float[] x2 = new float[ng];
    for (int ig=0; ig<ng; ig++) {
      x1[ig] = (float)s1.getValue(g1[ig]);
      x2[ig] = (float)u[g1[ig]];
    }
    return new float[][]{x1,x2};
  }
  
  public float[][] getX1X2AvgM2(
      float[] u, float dr1, Sampling s1, boolean useNG) 
  {
    s1 = (s1==null)?new Sampling(u.length):s1;
    float[] ea = new float[_ne1];
    float[] et;
    for (int i2=0; i2<_n2; i2++) {
      et = getEnvelope(_pp2[i2]);
      for (int i1=0; i1<_ne1; i1++) {
        ea[i1] += et[i1];
      }
    }
    int dg1 = (int)ceil(1.0f/dr1);
    int[] g1;
    if (useNG) {
      int ng = 1+(_ne1-1)/dg1;
      g1 = Subsample.subsample(ea,dg1,ng);
    } else
      g1 = Subsample.subsample(ea,dg1);
    int ng = g1.length;
    float[] x1 = new float[ng];
    float[] x2 = new float[ng];
    for (int ig=0; ig<ng; ig++) {
      x1[ig] = (float)s1.getValue(g1[ig]);
      x2[ig] = (float)u[g1[ig]];
    }
    return new float[][]{x1,x2};
  }
  
  public float[][] getX1X2AvgM3(
      float[] u, float dr1, Sampling s1, boolean useNG) 
  {
    s1 = (s1==null)?new Sampling(u.length):s1;
    float[] ea = new float[_ne1];
    float[] et;
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        et = getEnvelope(_pp3[i3][i2]);
        for (int i1=0; i1<_ne1; i1++) {
          ea[i1] += et[i1];
        }
      }
    }
    int dg1 = (int)ceil(1.0f/dr1);
    int[] g1;
    if (useNG) {
      int ng = 1+(_ne1-1)/dg1;
      g1 = Subsample.subsample(ea,dg1,ng);
    } else
      g1 = Subsample.subsample(ea,dg1);
    int ng = g1.length;
    float[] x1 = new float[ng];
    float[] x2 = new float[ng];
    for (int ig=0; ig<ng; ig++) {
      x1[ig] = (float)s1.getValue(g1[ig]);
      x2[ig] = (float)u[g1[ig]];
    }
    return new float[][]{x1,x2};
  }
  
  public float[][] getX1X2AvgMFlat2(
      float[] u, float[][] ppf, float dr1, Sampling s1, boolean useNG) 
  {
    s1 = (s1==null)?new Sampling(u.length):s1;
    int nppf = ppf[0].length;
    float[] ea = new float[nppf];
    float[] et;
    for (int i2=0; i2<_n2; i2++) {
      et = getEnvelope(ppf[i2]);
      for (int i1=0; i1<nppf; i1++) {
        ea[i1] += et[i1];
      }
    }
    int dg1 = (int)ceil(1.0f/dr1);
    int[] g1;
    if (useNG) {
      int ng = 1+(nppf-1)/dg1;
      g1 = Subsample.subsample(ea,dg1,ng);
    } else
      g1 = Subsample.subsample(ea,dg1);
    int ng = g1.length;
    float[] x1 = new float[ng];
    float[] x2 = new float[ng];
    for (int ig=0; ig<ng; ig++) {
      x1[ig] = (float)s1.getValue(g1[ig]);
      x2[ig] = (float)u[g1[ig]];
    }
    return new float[][]{x1,x2};
  }
  
  public float[][] getX1X2AvgMFlat3(
      float[] u, float[][][] ppf, float dr1, Sampling s1, boolean useNG) 
  {
    s1 = (s1==null)?new Sampling(u.length):s1;
    int nppf = ppf[0][0].length;
    float[] ea = new float[nppf];
    float[] et;
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        et = getEnvelope(ppf[i3][i2]);
        for (int i1=0; i1<nppf; i1++) {
          ea[i1] += et[i1];
        }
      }
    }
    int dg1 = (int)ceil(1.0f/dr1);
    int[] g1;
    if (useNG) {
      int ng = 1+(nppf-1)/dg1;
      g1 = Subsample.subsample(ea,dg1,ng);
    } else
      g1 = Subsample.subsample(ea,dg1);
    int ng = g1.length;
    float[] x1 = new float[ng];
    float[] x2 = new float[ng];
    for (int ig=0; ig<ng; ig++) {
      x1[ig] = (float)s1.getValue(g1[ig]);
      x2[ig] = (float)u[g1[ig]];
    }
    return new float[][]{x1,x2};
  }
  
  public float[][][] getX1X2(float[][] u, float dr1, Sampling s1) {
    s1 = (s1==null)?new Sampling(u.length):s1;
    int dg1 = (int)ceil(1.0f/dr1);
    int[][] g = new int[_n2][];
    int[] g1 = Subsample.subsample(_ne1,dg1);
    for (int i2=0; i2<_n2; i2++)
      g[i2] = copy(g1);
    int ng = g[0].length;
    float[][] x1 = new float[_n2][ng];
    float[][] x2 = new float[_n2][ng];
    for (int i2=0; i2<_n2; i2++) {
      for (int ig=0; ig<ng; ig++) {
        x1[i2][ig] = (float)s1.getValue(g[i2][ig]);
        x2[i2][ig] = (float)u[i2][g[i2][ig]];
      }  
    }
    return new float[][][]{x1,x2};
  }
  
  public float[][][] getX1X2T(float[][] u, float dr2, Sampling s2) {
    s2 = (s2==null)?new Sampling(u[0].length):s2;
    int dg2 = (int)ceil(1.0f/dr2);
    int[][] g = new int[_ne1][];
    int[] g2 = Subsample.subsample(_n2,dg2);
    for (int i1=0; i1<_ne1; i1++)
      g[i1] = copy(g2);
    int ng = g[0].length;
    float[][] x1 = new float[_ne1][ng];
    float[][] x2 = new float[_ne1][ng];
    for (int i1=0; i1<_ne1; i1++) {
      for (int ig=0; ig<ng; ig++) {
        x1[i1][ig] = (float)s2.getValue(g[i1][ig]);
        x2[i1][ig] = (float)u[i1][g[i1][ig]];
      }  
    }
    return new float[][][]{x1,x2};
  }
  
  public float[][][] getX1X2M(
      float[][] u, float[][] ppf, float[][] x1, float dr1, Sampling s1)
  {
    s1 = (s1==null)?new Sampling(u.length):s1;
//    float[][] ppf = fm.flatten(_pp2);
//    int ppfMax = (int)(fm.u1[0][_ne1]);
//    int nppf = ppfMax+1;
    int nppf = ppf[0].length;
    int dg1 = (int)ceil(1.0f/dr1);
    float[] ea = new float[nppf];
    float[] et;
    for (int i2=0; i2<_n2; i2++) {
      et = getEnvelope(ppf[i2]);
      for (int i1=0; i1<nppf; i1++)
        ea[i1] += et[i1];
    }
    int[] g1 = Subsample.subsample(ea,dg1);
    int ng1 = g1.length;
    int ng1m1 = ng1-1; 
    int[][] gw = new int[_n2][ng1];
    for (int i2=0; i2<_n2; i2++) {
      gw[i2][0] = 0;
      gw[i2][ng1m1] = _ne1-1;
      for (int i1=1; i1<ng1m1; i1++) {
        int u1 = g1[i1];
        gw[i2][i1] = (int)(x1[i2][u1]);
        // Make sure we don't pass the last grid point _ne1-1
        if (gw[i2][i1]>=(_ne1-1)) {
          int d = ng1m1-i1; // This many left will be out-of-bounds
          for (int id=d; id>0; id--)
            gw[i2][ng1m1-id] = gw[i2][ng1m1]-id;
          for (int i1p=1; i1p<=i1; i1p++) { // Adjust previous grid points
            while (gw[i2][i1-i1p]>=gw[i2][i1-i1p+1])
              gw[i2][i1-i1p] = gw[i2][i1-i1p]-1;
            assert(gw[i2][i1-i1p]<gw[i2][i1-i1p+1]):"i2="+i2+", i1="+i1;
          }
          break; // The rest of i1 indices are squeezed in, next trace.
        }
      }
    }
    float[][] x1c = new float[_n2][ng1];
    float[][] x2c = new float[_n2][ng1];
    for (int i2=0; i2<_n2; i2++) {
      for (int i1=0; i1<ng1; i1++) {
        x1c[i2][i1] = (float)s1.getValue(gw[i2][i1]);
        x2c[i2][i1] = (float)u[i2][gw[i2][i1]];
      }
    }
    return new float[][][]{x1c,x2c};
  }
  
  public float[] findShiftSums(float[][] u) {
    int nu = u.length;
    float[] sums = new float[nu];
    for (int iu=0; iu<nu; iu++) {
      sums[iu] = sum(u[iu]);
    }
    normalize(sums);
    return sums;
  }
  
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
   * Computes an array of VpVs ratios an array of shifts u.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values. 
   */
  public static float[] vpvs(float[] u) {
    int n = u.length;
    int nm1 = n-1;
    float[] vpvs = new float[n];
    vpvs[ 0 ] = 1.0f + 2*(u[ 1 ]-u[  0  ]); // at i1=0, forward difference
    vpvs[nm1] = 1.0f + 2*(u[nm1]-u[nm1-1]); // at i1=nm1, backward difference
    for (int i1=1; i1<nm1; ++i1) {
      vpvs[i1] = 1.0f + 2*((u[i1+1]-u[i1-1])/2.0f);
    }
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
      vpvs[i2][  0 ] = 1.0f + 2*(u[i2][  1 ]-u[i2][   0  ]); 
      vpvs[i2][n1m1] = 1.0f + 2*(u[i2][n1m1]-u[i2][n1m1-1]);
      for (int i1=1; i1<n1m1; ++i1) {
        vpvs[i2][i1] = 1.0f+2*((u[i2][i1+1]-u[i2][i1-1])/2.0f);
      }
    }
    return vpvs;
  }
  
  /**
   * Computes an array of VpVs ratios an array of shifts u.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values. 
   */
  public static float[][] vpvs(float[][] u, float[][] u1m) {
    int n1 = u[0].length;
    int n2 = u.length;
    int n1m1 = n1-1;
    float[][] vpvs = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      // At i1=0, forward difference. At i1=nm1, backward difference.
      vpvs[i2][  0 ] = 1.0f + 2*(u[i2][  1 ]-u[i2][   0  ]); 
      vpvs[i2][n1m1] = 1.0f + 2*(u[i2][n1m1]-u[i2][n1m1-1]);
      for (int i1=1; i1<n1m1; ++i1) {
        vpvs[i2][i1] = 1.0f+2*((u[i2][i1+1]-u[i2][i1-1])/2.0f);
      }
    }
    float[][] vpvsi = new float[n2][n1];
    float[] r = new float[n1];
    ramp(0.0f,1.0f,r);
    for (int i2=0; i2<n2; i2++) {
      float[] u1 = new float[n1];
      for (int i1=0; i1<n1; i1++)
        u1[i1] = u1m[i2][i1];
      CubicInterpolator civ = new CubicInterpolator(Method.LINEAR,r,vpvs[i2]);
      vpvsi[i2] = civ.interpolate(u1);
    }
    return vpvsi;
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
        vpvs[i3][i2][  0 ] = 1.0f + 2*(u[i3][i2][  1 ]-u[i3][i2][   0  ]); 
        vpvs[i3][i2][n1m1] = 1.0f + 2*(u[i3][i2][n1m1]-u[i3][i2][n1m1-1]);
        for (int i1=1; i1<n1m1; ++i1) {
          vpvs[i3][i2][i1] = 1.0f+2*((u[i3][i2][i1+1]-u[i3][i2][i1-1])/2.0f);
        }
      }  
    }
    return vpvs;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values. 
   */
  public static float[][][] vpvs(float[][][] u, float[][][] u1) {
    int n1 = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    int n1m1 = n1-1;
    float[][][] vpvs = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        // At i1=0, forward difference. At i1=nm1, backward difference.
        vpvs[i3][i2][  0 ] = 1.0f + 2*(u[i3][i2][  1 ]-u[i3][i2][   0  ]); 
        vpvs[i3][i2][n1m1] = 1.0f + 2*(u[i3][i2][n1m1]-u[i3][i2][n1m1-1]);
        for (int i1=1; i1<n1m1; ++i1) {
          vpvs[i3][i2][i1] = 1.0f+2*((u[i3][i2][i1+1]-u[i3][i2][i1-1])/2.0f);
        }
      }  
    }
    float[][][] vpvsi = new float[n3][n2][n1];
    float[] r = new float[n1];
    ramp(0.0f,1.0f,r);
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        float[] u1s = new float[n1];
        for (int i1=0; i1<n1; i1++)
          u1s[i1] = u1[i3][i2][i1];
        CubicInterpolator civ = new CubicInterpolator(r,vpvs[i3][i2]);
        vpvsi[i3][i2] = civ.interpolate(u1s);
      }
    }
    return vpvsi;
  }
  
  public static float[][] unflattenShifts(float[][] uf, float[][] u1m) {
    int n1 = uf[0].length;
    int n2 = uf.length;
    float[][] u = new float[n2][n1];
    float[] r = new float[n1];
    ramp(0.0f,1.0f,r);
    for (int i2=0; i2<n2; i2++) {
      float[] u1s = new float[n1];
      for (int i1=0; i1<n1; i1++)
        u1s[i1] = u1m[i2][i1];
      CubicInterpolator civ = new CubicInterpolator(Method.LINEAR,r,uf[i2]);
      u[i2] = civ.interpolate(u1s);
    }
    return u;
  }
  
  public static float[][][] unflattenShifts(float[][][] uf, float[][][] u1m) {
    int n1 = uf[0][0].length;
    int n2 = uf[0].length;
    int n3 = uf.length;
    float[][][] u = new float[n3][n2][n1];
    float[] r = new float[n1];
    ramp(0.0f,1.0f,r);
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        float[] u1s = new float[n1];
        for (int i1=0; i1<n1; i1++)
          u1s[i1] = u1m[i3][i2][i1];
        CubicInterpolator civ = new CubicInterpolator(r,uf[i3][i2]);
        u[i3][i2] = civ.interpolate(u1s);
      }
    }
    return u;
  }

  ///////////////////////////////////////////////////////////////////////////
  // for research and atypical applications

  public float[][] computeErrors() {
    float[][] e = new float[_ne1][_nel];
    computeErrors(_pp1,_ps1,e);
    normalizeErrors(e);
    return e;
  }
  
  public float[][] computeErrors(float[] pp, float[] ps) {
    int npp = pp.length;
    int nps = ps.length;
    int[] el = computeErrorLengths(npp,nps,0);
    float[][] e = new float[el[0]][el[1]];
    if (_fr>1) {
      int nx = (npp-1)*_fr+1;
      float[] psi = new float[nx];
      _si.interpolate(npp,1.0,0.0,ps,nx,_shifts1.getDelta(),0.0,psi);
      computeErrors(pp,psi,e);
    } else {
      computeErrors(pp,ps,e);
    }
    normalizeErrors(e);
    return e;
  }
  
  /**
   * Compute alignment errors for 2D images.
   * @return Normalized alignment errors.
   */
  public float[][][] computeErrors2() {
    final float[][][] e = new float[_n2][_ne1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrors(_pp2[i2],_ps2[i2],e[i2]);  
    }});
//    normalizeErrors(e);
    return e;
  }
  
  public float[][][][] computeErrors(float[][][] pp, float[][][] ps) {
    final int npp1 = pp[0][0].length;
    final int npp2 = pp[0].length;
    final int npp3 = pp.length;
    final int nps1 = ps[0][0].length;
    final int nps2 = ps[0].length;
    final int nps3 = ps.length;
    Check.argument(npp2==nps2,"npp2==nps2");
    Check.argument(npp3==nps3,"npp3==nps3");
    final float[][][] fpp = pp;
    final float[][][] fps = ps;
    int[] el = computeErrorLengths(npp1,nps1,0);
    final float[][][][] e = new float[npp3][npp2][el[0]][el[1]];
    Parallel.loop(npp3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<npp2; i2++) {
        computeErrors(fpp[i3][i2],fps[i3][i2],e[i3][i2]);    
      }
    }});
    normalizeErrors(e);
    return e;
  }
  
  /**
   * Compute average alignment errors for 2D images.
   * @return Normalized averaged alignment errors for 2D images. 
   */
  public float[][] computeErrors2Average() {
    float[][] e = Parallel.reduce(_n2,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(_pp2[i2],_ps2[i2],e);
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    normalizeErrors(e);
    return e;
  }
  
  public float[][] computeErrors2Average(float[][] x1) {
    final float[][][] e2 = new float[_n2][_ne1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrors(_pp2[i2],_ps2[i2],e2[i2]);  
    }});
    float[][] e = e2[0];
    for (int i2=1; i2<_n2; i2++) {
      for (int i1=0; i1<_ne1; i1++) {
        int x = (int)(x1[i2][i1]+0.5);
        e[i1] = add(e[i1],e2[i2][x]);  
      }
    }
    normalizeErrors(e);
    return e;
  }
  
  /**
   * Compute average alignment errors for 3D images.
   * @return Normalized averaged alignment errors for 3D images. 
   */
  public float[][] computeErrors3Average() {
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
    normalizeErrors(e);
    return e;
  }
  
  public float[][] computeErrors3Average(float[][][] x1) {
    final float[][][][] e3 = new float[_n3][_n2][_ne1][_nel];
    Parallel.loop(_n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; i2++)
        computeErrors(_pp3[i3][i2],_ps3[i3][i2],e3[i3][i2]);  
    }});
    float[][] e = zerofloat(_nel,_ne1);
    for (int i3=0; i3<_n3; i3++) {
      for (int i2=0; i2<_n2; i2++) {
        for (int i1=0; i1<_ne1; i1++) {
          int x = (int)(x1[i3][i2][i1]+0.5);
          e[i1] = add(e[i1],e3[i3][i2][x]);  
        }
      }
    }
    normalizeErrors(e);
    return e;
  }
  
  public float[][] computeErrors1(float[][] pp, float[][] ps1) {
    final int n1 = pp[0].length;
    int n1ps = ps1[0].length;
    int[] el = computeErrorLengths(n1,n1ps,0);
    final int n1M = el[0]; // maximum pp time index
    final int nl1 = el[1]; // number of pp-ps1 lags
    final int n2 = pp.length;
    final float[][] fpp = pp;
    final float[][] fps1 = ps1;
    float[][] e = Parallel.reduce(n2, new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i2) {
      float[][] e = new float[n1M][nl1];
      if (_fr>1) {
        int nx = (n1-1)*_fr+1;
        float[] ps1i = new float[nx];
        _si.interpolate(n1,1.0,0.0,fps1[i2],nx,_shifts1.getDelta(),0.0,ps1i);
        computeErrors(fpp[i2],ps1i,e);
      } else {
        computeErrors(fpp[i2],fps1[i2],e);
      }  
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    normalizeErrors(e);
    return e;
  }
  
  public float[][] computeErrors1(float[][][] pp, float[][][] ps) {
    final int n1pp = pp[0][0].length;
    final int n1ps = ps[0][0].length;
    final int n2 = pp[0].length;
    final int n3 = pp.length;
    int[] el = computeErrorLengths(n1pp,n1ps,0);
    final int n1M = el[0]; // maximum pp time index
    final int nl1 = el[1]; // number of pp-ps1 lags
    final float[][][] fpp = pp;
    final float[][][] fps1 = ps;
    float[][] e = Parallel.reduce(n2*n3,new Parallel.ReduceInt<float[][]>() {
      public float[][] compute(int i23) {
        int i2 = i23%n2;
        int i3 = i23/n2;
        float[][] e = new float[n1M][nl1];
        if (_fr>1) {
          int nx = (n1pp-1)*_fr+1;
          float[] ps1i = new float[nx];
          _si.interpolate(n1pp,1.0,0.0,fps1[i3][i2],
              nx,_shifts1.getDelta(),0.0,ps1i);
          computeErrors(fpp[i3][i2],ps1i,e);
        } else {
          computeErrors(fpp[i3][i2],fps1[i3][i2],e);
        }  
        return e;
      }
      public float[][] combine(float[][] ea, float[][] eb) {
        return add(ea,eb);
      }});
    normalizeErrors(e);
    return e;
  }
  
  public float[][] computeErrors(float[] ps1, float[] ps2, int maxShift) {
    int n1 = ps1.length;
    int[] el = computeErrorLengths(n1,maxShift);
    float[][] e = new float[el[0]][el[1]];
    if (_fr>1) {
      int nx = (n1-1)*_fr+1;
      float[] ps2i = new float[nx];
      _si.interpolate(n1,1.0,0.0,ps2,nx,_shifts1.getDelta(),0.0,ps2i);
      computeErrors(ps1,ps2i,e);
    } else {
      computeErrors(ps1,ps2,e);
    }  
    normalizeErrors(e);
    return e;
  }
  
  public float[][][] computeErrors(
      float[] pp, float[] ps1, float[] ps2, float fsp)
  {
    int n1 = pp.length;
    int[] el = computeErrorLengths(n1,fsp);
    float[][][] e = new float[el[0]][el[1]][el[2]];
    if (_fr>1) {
      int nx = (n1-1)*_fr+1;
      float[] ps2i = new float[nx];
      _si.interpolate(n1,1.0,0.0,ps2,nx,_shiftsS.getDelta(),0.0,ps2i);
      computeErrors(pp,ps1,ps2i,e);
    } else {
      computeErrors(pp,ps1,ps2,e);
    }
    normalizeErrors(e);
    return e;
  }

  public float[][] computeErrors1(float[][] ps1, float[][] ps2, int maxShift) {
    final int n1 = ps1[0].length;
    int[] el = computeErrorLengths(n1,maxShift);
    final int n1M = el[0]; // maximum pp time index
    final int nl1 = el[1]; // number of pp-ps1 lags
    final int n2 = ps1.length;
    final float[][] fpp = ps1;
    final float[][] fps2 = ps2;
    float[][] e = Parallel.reduce(n2, new Parallel.ReduceInt<float[][]>() {
      public float[][] compute(int i2) {
        float[][] e = new float[n1M][nl1];
        if (_fr>1) {
          int nx = (n1-1)*_fr+1;
          float[] ps2i = new float[nx];
          _si.interpolate(n1,1.0,0.0,fps2[i2],nx,_shifts1.getDelta(),0.0,ps2i);
          computeErrors(fpp[i2],ps2i,e);
        } else {
          computeErrors(fpp[i2],fps2[i2],e);
        }
        return e;
      }
      public float[][] combine(float[][] ea, float[][] eb) {
        return add(ea,eb);
      }});
      normalizeErrors(e);
    return e;
  }
  
  public float[][][] computeErrors1(
      float[][] pp, float[][] ps1, float[][] ps2, float fsp)
  {
    final int n1 = pp[0].length;
    int[] el = computeErrorLengths(n1,fsp);
    final int n1M = el[0]; // maximum pp time index
    final int nl1 = el[1]; // number of pp-ps1 lags
    final int nlS = el[2]; // number of ps1-ps2 lags
    final int n2 = pp.length;
    final float[][] fpp = pp;
    final float[][] fps1 = ps1;
    final float[][] fps2 = ps2;
    float[][][] e = Parallel.reduce(n2, new Parallel.ReduceInt<float[][][]>() {
    public float[][][] compute(int i2) {
      float[][][] e = new float[n1M][nl1][nlS];
      if (_fr>1) {
        int nx = (n1-1)*_fr+1;
        float[] ps2i = new float[nx];
        _si.interpolate(n1,1.0,0.0,fps2[i2],nx,_shifts1.getDelta(),0.0,ps2i);
        computeErrors(fpp[i2],fps1[i2],ps2i,e);
      } else {
        computeErrors(fpp[i2],fps1[i2],fps2[i2],e);  
      }
      return e;
    }
    public float[][][] combine(float[][][] ea, float[][][] eb) {
      return add(ea,eb);
    }});
    normalizeErrors(e);
    return e;
  }

  public static float[][] smoothErrors(float[][] e) {
    int n1 = e.length;
    int nl = e[0].length;
    float[][] ef = like(e);
    float[][] er = like(e);
    float[][] es = like(e);
    accumulate( 1,e,ef);
    accumulate(-1,e,er);
    for (int i1=0; i1<n1; i1++) {
      for (int il=0; il<nl; il++) {
        es[i1][il] = ef[i1][il]+er[i1][il]-e[i1][il];
      }
    }
    normalizeErrors(es);
    return es;
  }
  
  /**
   * Returns smooth alignment errors on a sparse grid. 
   * @param e 2D array of alignment errors.
   * @param rmin minimum slope in first dimension.
   * @param rmax maximum slope in first dimension.
   * @param g sparse grid points in first dimension.
   * @return smoothed alignment errors.
   */
  public static float[][] smoothErrorsSparse(
      float[][] e, float rmin, float rmax, int[] g)      
  {
    int ng = g.length;
    int nel = e[0].length;
    float[][] ef = new float[ng][nel];
    float[][] er = new float[ng][nel];
    float[][] es = new float[ng][nel];
    accumulateSparseNew( 1,rmin,rmax,g,e,ef,null);
    accumulateSparseNew(-1,rmin,rmax,g,e,er,null);
    float scale = 1.0f/e.length;
    for (int i1=0; i1<ng; i1++) {
      for (int il=0; il<nel; il++) {
        es[i1][il] = scale*(ef[i1][il]+er[i1][il]-e[g[i1]][il]);
      }
    }
    return es;
  }
  
  public float[][][] smoothErrors1(final float[][] pp, final float[][] ps,
      final float r1min, final float r1max, final int[] g1)
  {
    final int ng1 = g1.length;
    final float[][][] es1 = new float[_n2][ng1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(pp[i2],ps[i2],e);
      es1[i2] = smoothErrorsSparse(e,r1min,r1max,g1);
    }});
    return es1;
  }
  
  public float[][][] smoothErrors1(final float[][] pp, final float[][] ps,
      final float r1min, final float r1max, final int[][] g1)
  {
    final int ng1 = g1.length;
    final float[][][] es1 = new float[_n2][ng1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(pp[i2],ps[i2],e);
      es1[i2] = smoothErrorsSparse(e,r1min,r1max,g1[i2]);
    }});
    return es1;
  }
  
  public float[][][][] smoothErrors1(final float[][][] pp, final float[][][] ps,
      final float r1min, final float r1max, final int[] g1)
  {
    final int ng1 = g1.length;
    final float[][][][] es1 = new float[_n3][_n2][ng1][_nel];
    Parallel.loop(_n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; i2++) {
        float[][] e = new float[_ne1][_nel];
        computeErrors(pp[i3][i2],ps[i3][i2],e);
        es1[i3][i2] = smoothErrorsSparse(e,r1min,r1max,g1);  
      }
    }});
    return es1;
  }
  
  public float[][][] smoothErrors2(
      final float[][][] es1, 
      final float r2min, final float r2max, final int[] g2)
  {
    final int ng1 = es1[0].length;
    final int ng2 = g2.length;
    final float[][][] es = new float[ng2][ng1][_nel]; // smoothed errors
    Parallel.loop(ng1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][]  e2 = new float[_n2][_nel]; // errors at index i1
      for (int i2=0; i2<_n2; ++i2)
        e2[i2] = es1[i2][i1];
      float[][] es2 = smoothErrorsSparse(e2,r2min,r2max,g2);
      for (int i2=0; i2<ng2; i2++)
        es[i2][i1] = es2[i2];
    }});
    return es;
  }

  /**
   * Returns smooth alignment errors for 2D images on a sparse grid in
   * each dimension.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param g1 sparse grid points in first dimension.
   * @param r2max maximum slope in second dimension.
   * @param g2 sparse grid points in second dimension.
   * @return smoothed alignment errors.
   */
  public float[][][] smoothErrorsSparse(
      final float[][] pp, final float[][] ps,
      final float r1min, final float r1max, final int[] g1,
      final float r2min, final float r2max, final int[] g2)
  {
    // Smooth alignment errors in first dimension for each trace.
    final int ng1 = g1.length;
    Stopwatch s = new Stopwatch();
    s.start();
    final float[][][] es1 = new float[_n2][ng1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(pp[i2],ps[i2],e);
      es1[i2] = smoothErrorsSparse(e,r1min,r1max,g1);
    }});
    System.out.println("Finished 1st Dimension Smoothing in "+s.time()+
        " seconds");
    normalizeErrors(es1);
    s.restart();
    
    // Smooth in the second dimension
    final int ng2 = g2.length;
    final float[][][] es = new float[ng2][ng1][_nel]; // smoothed errors
    final Parallel.Unsafe<float[][][]> e2u = 
        new Parallel.Unsafe<float[][][]>();
    Parallel.loop(ng1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][][] ee = e2u.get();
      if (ee==null) e2u.set(ee=new float[3][ng2][_nel]);
      float[][] es2 = ee[0]; // smooth sparse errors at index i1 
      float[][] ef2 = ee[1]; // forward errors
      float[][] er2 = ee[2]; // reverse errors
      float[][]  e2 = new float[_n2][_nel]; // smooth errors at index i1
      for (int i2=0; i2<ng2; ++i2) {
        es2[i2] = es[i2][i1];
        for (int il=0; il<_nel; ++il) {
          ef2[i2][il] = 0.0f;
          er2[i2][il] = 0.0f;
        }
      }
      for (int i2=0; i2<_n2; ++i2) {
        e2[i2] = es1[i2][i1];
      }
      accumulateSparseNew( 1,r2min,r2max,g2,e2,ef2,null);
      accumulateSparseNew(-1,r2min,r2max,g2,e2,er2,null);
      float scale = 1.0f/ng2;
      for (int i2=0; i2<ng2; ++i2) {
        for (int il=0; il<_nel; ++il) {
          es2[i2][il] = scale*(ef2[i2][il]+er2[i2][il]-e2[g2[i2]][il]);
        }
      }
    }});
    normalizeErrors(es);
    System.out.println("Finished 2nd Dimension Smoothing in "+s.time()+
        " seconds");
    return es;
  }
  
  /**
   * Returns smooth alignment errors for 3D images on a sparse grid in
   * each dimension.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param g1 sparse grid points in first dimension.
   * @param r2max maximum slope in second dimension.
   * @param g2 sparse grid points in second dimension.
   * @param r3min minimum slope in third dimension.
   * @param r3max maximum slope in third dimension.
   * @param g3 sparse grid points in third dimension.
   * @return smoothed alignment errors.
   */
  public float[][][][] smoothErrorsSparse(
      final float r1min, final float r1max, final int[] g1,
      final float r2min, final float r2max, final int[] g2,
      final float r3min, final float r3max, final int[] g3)
  {
    final int ng1 = g1.length;
    final int ng2 = g2.length;
    final int ng3 = g3.length;
    final float[][][][] es12 = new float[_n3][ng2][ng1][_nel]; // smoothed errors
    for (int i3=0; i3<_n3; i3++) {
      System.out.println("Starting i3 "+i3);
      es12[i3] = 
          smoothErrorsSparse(_pp3[i3],_ps3[i3],r1min,r1max,g1,r2min,r2max,g2);
    }
    normalizeErrors(es12);

    final float[][][][] es = new float[ng3][ng2][ng1][_nel]; // smoothed errors
    final Parallel.Unsafe<float[][][]> e3u = 
        new Parallel.Unsafe<float[][][]>();
    Parallel.loop(ng1,new Parallel.LoopInt() {
    public void compute(int i1) {
      for (int i2=0; i2<ng2; i2++) {
        float[][][] ee = e3u.get();
        if (ee==null) e3u.set(ee=new float[3][ng3][_nel]);
        float[][] es3 = ee[0]; // smooth sparse errors at index i1 
        float[][] ef3 = ee[1]; // forward errors
        float[][] er3 = ee[2]; // reverse errors
        float[][]  e3 = new float[_n3][_nel]; // smooth errors at index i1
        for (int i3=0; i3<ng3; i3++) {
          es3[i3] = es[i3][i2][i1];
          for (int il=0; il<_nel; ++il) {
            ef3[i3][il] = 0.0f;
            er3[i3][il] = 0.0f;
          }
        }
        for (int i3=0; i3<_n3; i3++) {
          e3[i3] = es12[i3][i2][i1];
        }
        accumulateSparseNew( 1,r3min,r3max,g3,e3,ef3,null);
        accumulateSparseNew(-1,r3min,r3max,g3,e3,er3,null);
        for (int i3=0; i3<ng3; i3++) {
          for (int il=0; il<_nel; il++) {
            es3[i3][il] = ef3[i3][il]+er3[i3][il]-e3[g3[i3]][il];
          }
        }  
      }
    }});
    normalizeErrors(es);
    return es;
  }
  
  /**
   * Returns smooth alignment errors for 2D images on a sparse grid in
   * each dimension.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension.
   * @param g1 sparse grid points in first dimension.
   * @param r2max maximum slope in second dimension.
   * @param g2 sparse grid points in second dimension.
   * @return smoothed alignment errors.
   */
  public float[][][] smoothErrorsSparseM(
      final float[][] pp, final float[][] ps,
      final float r1min, final float r1max, final int[][] g1,
      final float r2min, final float r2max, final int[] g2)
  {
    // Smooth alignment errors in first dimension for each trace.
    Stopwatch s = new Stopwatch();
    s.start();
    int ng1 = g1[0].length;
    final float[][][] es1 = new float[_n2][ng1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(pp[i2],ps[i2],e);
      es1[i2] = smoothErrorsSparse(e,r1min,r1max,g1[i2]);
    }});
    System.out.println("Finished 1st Dimension Smoothing in "+s.time()+
        " seconds");
    normalizeErrors(es1);
    
    s.restart();
    final int ng2 = g2.length;
    final float[][][] es = new float[ng2][ng1][_nel]; // smoothed errors
    final Parallel.Unsafe<float[][][]> e1u = 
        new Parallel.Unsafe<float[][][]>();
    Parallel.loop(ng1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][][] ee = e1u.get();
      if (ee==null) e1u.set(ee=new float[3][ng2][_nel]);
      float[][]   m = ee[0]; // moves
      float[][] es2 = ee[1]; // smooth sparse errors at index i1 
      float[][] ef2 = ee[2]; // forward errors
      float[][] er2 = ee[3]; // reverse errors
      float[][]  e2 = new float[_n2][_nel]; // smooth errors at index i1
      for (int i2=0; i2<ng2; ++i2) {
        es2[i2] = es[i2][i1];
        for (int il=0; il<_nel; ++il) {
          ef2[i2][il] = 0.0f;
          er2[i2][il] = 0.0f;
        }
      }
      for (int i2=0; i2<_n2; ++i2) {
        e2[i2] = es1[i2][i1];
      }
      //TODO remove
//      accumulateSparse2( 1,r2max,g2,e2,ef2,m);
//      accumulateSparse2(-1,r2max,g2,e2,er2,m);
      accumulateSparseNew( 1,r2min,r2max,g2,e2,ef2,null);
      accumulateSparseNew(-1,r2min,r2max,g2,e2,er2,null);
      float scale = 1.0f/ng2;
      for (int i2=0; i2<ng2; ++i2) {
        for (int il=0; il<_nel; ++il) {
          es2[i2][il] = scale*(ef2[i2][il]+er2[i2][il]-e2[g2[i2]][il]);
        }
      }
    }});
    normalizeErrors(es);
    System.out.println("Finished 2nd Dimension Smoothing in "+s.time()+
        " seconds");
    return es;
  }
  
  public float[][][][] smoothErrorsSparseM(
      final float r1min, final float r1max, final int[][][] g1,
      final float r2min, final float r2max, final int[] g2,
      final float r3min, final float r3max, final int[] g3)
  {
    final int ng1 = g1[0][0].length;
    final int ng2 = g2.length;
    final int ng3 = g3.length;
    final float[][][][] es12 = new float[_n3][ng2][ng1][_nel]; // smoothed errors
    for (int i3=0; i3<_n3; i3++) {
      System.out.println("Starting i3 "+i3);
      es12[i3] = 
          smoothErrorsSparseM(_pp3[i3],_ps3[i3],r1min,r1max,g1[i3],r2min,r2max,g2);
    }
    normalizeErrors(es12);

    final float[][][][] es = new float[ng3][ng2][ng1][_nel]; // smoothed errors
    final Parallel.Unsafe<float[][][]> e3u = 
        new Parallel.Unsafe<float[][][]>();
    Parallel.loop(ng1,new Parallel.LoopInt() {
    public void compute(int i1) {
      for (int i2=0; i2<ng2; i2++) {
        float[][][] ee = e3u.get();
        if (ee==null) e3u.set(ee=new float[4][ng3][_nel]);
        float[][]   m = ee[0]; // moves
        float[][] es3 = ee[1]; // smooth sparse errors at index i1 
        float[][] ef3 = ee[2]; // forward errors
        float[][] er3 = ee[3]; // reverse errors
        float[][]  e3 = new float[_n3][_nel]; // smooth errors at index i1
        for (int i3=0; i3<ng3; i3++) {
          es3[i3] = es[i3][i2][i1];
          for (int il=0; il<_nel; ++il) {
            ef3[i3][il] = 0.0f;
            er3[i3][il] = 0.0f;
          }
        }
        for (int i3=0; i3<_n3; i3++) {
          e3[i3] = es12[i3][i2][i1];
        }
//        accumulateSparse2( 1,r3max,g3,e3,ef3,m);
//        accumulateSparse2(-1,r3max,g3,e3,er3,m);
        accumulateSparseNew( 1,r3min,r3max,g3,e3,ef3,null);
        accumulateSparseNew(-1,r3min,r3max,g3,e3,er3,null);
        for (int i3=0; i3<ng3; i3++) {
          for (int il=0; il<_nel; il++) {
            es3[i3][il] = ef3[i3][il]+er3[i3][il]-e3[g3[i3]][il];
          }
        }  
      }
    }});
    normalizeErrors(es);
    return es;
  }
  
  public float[][][] accumulateForward(float[][] e) {
    float[][] d = like(e);
    float[][] m = like(e);
    accumulate(1,e,d,m);
    normalizeErrors(d);
    return new float[][][]{d,m};
  }
  
  public float[][][] accumulateReverse(float[][] e) {
    float[][] d = like(e);
    float[][] m = like(e);
    accumulate(-1,e,d,m);
    return new float[][][]{d,m};
  }
  
  public float[][][] accumulateForward(
      float[][] e, int[] g, float rmin, float rmax)
  {
    float[][] d = like(e);
    float[][] m = like(e);
    accumulate(1,e,d,m,g,rmin,rmax);
//    normalizeErrors(d);
    return new float[][][]{d,m};
  }
  
  public static float[][][] accumulateForwardSparse(
      float[][] e, float rmin, float rmax, int[] g)
  {
    int nl = e[0].length;
    int ng = g.length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
//    Map<Integer,float[][]> mxbMap = Subsample.getMXB(g,rmin,rmax,false);
//    System.out.println("slope map size="+mxbMap.size());
    Stopwatch s = new Stopwatch();
    s.start();
//    accumulateSparse(1,rmin,rmax,g,mxbMap,e,d,m);
    accumulateSparseNew(1,rmin,rmax,g,e,d,m);
    System.out.println("accumulateSparse time: "+s.time()+" seconds");
    normalizeErrors(d);
    return new float[][][]{d,m};
  }
  
  public static float[][][] accumulateReverseSparse(
      float[][] e, float rmin, float rmax, int[] g)
  {
    int nl = e[0].length;
    int ng = g.length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
//    Map<Integer,float[][]> mxbMap = Subsample.getMXB(g,rmin,rmax,true);
//    System.out.println("slope map size="+mxbMap.size());
    Stopwatch s = new Stopwatch();
    s.start();
//    accumulateSparse(-1,rmin,rmax,g,mxbMap,e,d,m);
    accumulateSparseNew(-1,rmin,rmax,g,e,d,m);
    System.out.println("accumulateSparse time: "+s.time()+" seconds");
    normalizeErrors(d);
    return new float[][][]{d,m};
  }
  
  public float[] backtrackReverse(float[][] m) {
    float[] u = new float[m.length];
    backtrack(-1,_shifts1,m,u);
    return u;
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

//  public float[] backtrackForward(float[][] d, float[][] e) {
//    float[] u = new float[d.length];
//    backtrack(1,1,_shifts1,d,e,u);
//    return u;
//  }

  public float[] backtrackReverse2(float[][] s, float[][] m) {
    int n1 = s.length;
    int n1m = n1-1;
    float[] u = new float[n1];
    int nl = s[0].length;
    int ml = 0;
    float min = m[n1m][ml];
    for (int il=1; il<nl; il++) {
      if (m[n1m][il]<min) {
        min = m[n1m][il];
        ml = il;
      }
    }
    u[n1m] = s[n1m][ml];
    for (int i1=n1m-1; i1>=0; i1--) {
      int il = ml+(int)m[i1][ml];
      u[i1] = s[i1][il];
      ml = il;
    }
    return u;
  }
  
//  public float[][] backtrackReverse(float[][][] d, float[][][] e) {
//    float[] u1 = new float[d.length];
//    float[] uS = new float[d.length];
//    backtrackReverse(d,e,u1,uS);
//    return new float[][]{u1,uS};
//  }

  /**
   * Uses cubic interpolation to interpolate u, sampled on grid g,
   * to a uniformly sampled array of length n.
   * @param n
   * @param g
   * @param u
   * @return the interpolated shifts.
   */
  public static float[] interpolateSparseShifts(
      int n, int[] g, float[] u, Method m) 
  {
    int ng = g.length;
    float[] gf = new float[ng];
    float[] ui = new float[n ];
    for (int ig=0; ig<ng; ig++)
      gf[ig] = (float)g[ig];
    CubicInterpolator ci = new CubicInterpolator(m,ng,gf,u);
    for (int i=0; i<n; i++)
      ui[i] = ci.interpolate(i);
    return ui;
  }

  /**
   * Interpolates sparse shifts on a uniform grid (g1 indices are
   * constant for all g2). Computes the derivatives of the shifts
   * in the first dimension, interpolates the shifts and derivatives
   * in the n2 dimension with spline interpolation, and finally, 
   * constructs a cubic interpolator from those values to interpolate
   * the n1 values.
   * @param n1
   * @param n2
   * @param g1
   * @param g2
   * @param u
   * @return
   */
  public static float[][] interpolateSparseShifts(
      int n1, int n2, int[] g1, int[] g2, float[][] u)
  {
    int ng1 = g1.length;
    int ng2 = g2.length;
    Check.argument(ng1==u[0].length,"ng1==u[0].length: "+ng1+"=="+u[0].length);
    Check.argument(ng2==u.length,"ng2==u.length: "+ng2+"=="+u.length);
    float[] g1f = new float[ng1];
    float[] g2f = new float[ng2];
    for (int ig=0; ig<ng1; ig++)
      g1f[ig] = (float)g1[ig];
    for (int ig=0; ig<ng2; ig++)
      g2f[ig] = (float)g2[ig];

    // compute derivative in first dimension
    float[][] d = new float[ng2][ng1];
    for (int i2=0; i2<ng2; i2++) {
      CubicInterpolator ci = new CubicInterpolator(Method.LINEAR,g1f,u[i2]);
      ci.interpolate1(g1f,d[i2]);
    }

    // interpolate shifts and derivatives in the second dimension
    float[] u2 = new float[ng2];
    float[] d2 = new float[ng2];
    float[][] ui2 = new float[n2][ng1];
    float[][] di2 = new float[n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i2=0; i2<ng2; i2++) {
        u2[i2] = u[i2][i1];
        d2[i2] = d[i2][i1];
      }
      CubicInterpolator ciu = new CubicInterpolator(Method.SPLINE,g2f,u2);
      CubicInterpolator cid = new CubicInterpolator(Method.SPLINE,g2f,d2);
      for (int i2=0; i2<n2; i2++) {
        ui2[i2][i1] = ciu.interpolate(i2);
        di2[i2][i1] = cid.interpolate(i2);
      }
    }
    
    // interpolate in the first dimension, using derivatives
    float[][] ui = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      CubicInterpolator ci = new CubicInterpolator(g1f,ui2[i2],di2[i2]);
      for (int i1=0; i1<n1; i1++) {
        ui[i2][i1] = ci.interpolate(i1);
      }
    }
    return ui;
  }
  
  /**
   * Interpolates sparse shifts on a uniform grid (g1 indices are
   * constant for all g2). Computes the derivatives of the shifts
   * in the first dimension, interpolates the shifts and derivatives
   * in the n2 dimension with spline interpolation, and finally, 
   * constructs a cubic interpolator from those values to interpolate
   * the n1 values.
   * @param n1
   * @param n2
   * @param g1
   * @param g2
   * @param u
   * @return
   */
  public static float[][][] interpolateSparseShifts(
      int n1, int n2, int n3, int[] g1, int[] g2, int[] g3, float[][][] u)
  {
    int ng1 = g1.length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    Check.argument(
        ng1==u[0][0].length,"ng1==u[0][0].length: "+ng1+"=="+u[0][0].length);
    Check.argument(ng2==u[0].length,"ng2==u[0].length: "+ng2+"=="+u[0].length);
    Check.argument(ng3==u.length,"ng3==u.length: "+ng3+"=="+u.length);
    float[] g1f = new float[ng1];
    for (int ig1=0; ig1<ng1; ig1++)
      g1f[ig1] = (float)g1[ig1];
    float[] g2f = new float[ng2];
    for (int ig2=0; ig2<ng2; ig2++)
      g2f[ig2] = (float)g2[ig2];
    float[] g3f = new float[ng3];
    for (int ig3=0; ig3<ng3; ig3++)
      g3f[ig3] = (float)g3[ig3];
    
    // compute derivative in first dimension
    float[][][] d = new float[ng3][ng2][ng1];
    for (int i3=0; i3<ng3; i3++) {
      for (int i2=0; i2<ng2; i2++) {
        CubicInterpolator ci = 
            new CubicInterpolator(Method.LINEAR,g1f,u[0][i2]);
        ci.interpolate1(g1f,d[i3][i2]);
      }
    }
    
    // interpolate shifts and derivatives in the second dimension
    float[][] u23 = new float[ng3][ng2];
    float[][] d23 = new float[ng3][ng2];
    float[][][] ui23 = new float[n3][n2][ng1];
    float[][][] di23 = new float[n3][n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i3=0; i3<ng3; i3++) {
        for (int i2=0; i2<ng2; i2++) {
          u23[i3][i2] = u[i3][i2][i1];
          d23[i3][i2] = d[i3][i2][i1];
        }
      }
      BicubicInterpolator2 ciu = new BicubicInterpolator2(
          BicubicInterpolator2.Method.SPLINE,
          BicubicInterpolator2.Method.SPLINE,
          g2f,g3f,u23);
      BicubicInterpolator2 cid = new BicubicInterpolator2(
          BicubicInterpolator2.Method.SPLINE,
          BicubicInterpolator2.Method.SPLINE,
          g2f,g3f,d23);
      for (int i3=0; i3<n3; i3++) {
//        CubicInterpolator ciu = new CubicInterpolator(Method.SPLINE,g2f,u23[i3]);
//        CubicInterpolator cid = new CubicInterpolator(Method.SPLINE,g2f,d23[i3]);
        
        for (int i2=0; i2<n2; i2++) {
          ui23[i3][i2][i1] = ciu.interpolate(i2,i3);
          di23[i3][i2][i1] = cid.interpolate(i2,i3);
//          ui23[i3][i2][i1] = ciu.interpolate(i2);
//          di23[i3][i2][i1] = cid.interpolate(i2);
        }  
      }
    }
    
    // interpolate in the first dimension, using derivatives
    float[][][] ui = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        CubicInterpolator ci = 
            new CubicInterpolator(g1f,ui23[i3][i2],di23[i3][i2]);
        for (int i1=0; i1<n1; i1++) {
          ui[i3][i2][i1] = ci.interpolate(i1);
        }
      }
    }
    return ui;
  }
  
  /**
   * Interpolates sparse shifts on an irregular grid (g1 indices are
   * not constant for all g2). Computes the derivatives of the shifts
   * in the first dimension, interpolates the shifts and derivatives
   * in the n2 dimension with spline interpolation, and finally, 
   * constructs a cubic interpolator from those values to interpolate
   * the n1 values.
   * @param n1
   * @param n2
   * @param g1Flat
   * @param g2
   * @param u
   * @param x1
   * @return
   */
  public static float[][] interpolateSparseShifts(
      int n1, int n2, int[] g1Flat, int[] g2, float[][] u, float[][] x1)
  {
    int ng1 = g1Flat.length;
    int ng2 = g2.length;
    Check.argument(ng1==u[0].length,"ng1==u[0].length: "+ng1+"=="+u[0].length);
    Check.argument(ng2==u.length,"ng2==u.length: "+ng2+"=="+u.length);
    float[][] g1f = new float[n2][ng1];
    for (int i2=0; i2<n2; i2++) {
      for (int ig=0; ig<ng1; ig++)
        g1f[i2][ig] = x1[i2][g1Flat[ig]];
    }
    float[] g2f = new float[ng2];
    for (int ig=0; ig<ng2; ig++)
      g2f[ig] = (float)g2[ig];
    
    // compute derivative in first dimension
    float[][] d = new float[ng2][ng1];
    for (int i2=0; i2<ng2; i2++) {
      CubicInterpolator ci = 
          new CubicInterpolator(Method.LINEAR,g1f[g2[i2]],u[i2]);
      ci.interpolate1(g1f[g2[i2]],d[i2]);
    }
    
    // interpolate shifts and derivatives in the second dimension
    float[] u2 = new float[ng2];
    float[] d2 = new float[ng2];
//    float[] g12 = new float[ng2];
    float[][] ui2 = new float[n2][ng1];
    float[][] di2 = new float[n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i2=0; i2<ng2; i2++) {
        u2[i2] = u[i2][i1];
        d2[i2] = d[i2][i1];
//        g12[i2] = g1f[g2[i2]][i1];
      }
//      if (i1==11 || i1==12) {
//        dump(u2);
//        dump(d2);
//        dump(g12);
//      }
      CubicInterpolator ciu = new CubicInterpolator(Method.LINEAR,g2f,u2);
      CubicInterpolator cid = new CubicInterpolator(Method.LINEAR,g2f,d2);
      for (int i2=0; i2<n2; i2++) {
        ui2[i2][i1] = ciu.interpolate(i2);
        di2[i2][i1] = cid.interpolate(i2);
      }
    }
    
//    for (int i2=0; i2<ng2; i2++) {
//      for (int i1=0; i1<ng1; i1++) {
//        assert u[i2][i1]==ui2[g2[i2]][i1];
//        assert d[i2][i1]==di2[g2[i2]][i1];
//      }
//    }
    
    // interpolate in the first dimension, using derivatives
    float[][] ui = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      CubicInterpolator ci = new CubicInterpolator(g1f[i2],ui2[i2],di2[i2]);
      for (int i1=0; i1<n1; i1++) {
        ui[i2][i1] = ci.interpolate(i1);
      }
    }
    return ui;
  }
  
  public static float[][][] interpolateSparseShifts(
      int n1, int n2, int n3, int[] g1Flat, int[] g2, int[] g3, 
      float[][][] u, float[][][] x1)
  {
    int ng1 = g1Flat.length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    Check.argument(
        ng1==u[0][0].length,"ng1==u[0][0].length: "+ng1+"=="+u[0][0].length);
    Check.argument(ng2==u[0].length,"ng2==u[0].length: "+ng2+"=="+u[0].length);
    Check.argument(ng3==u.length,"ng3==u.length: "+ng3+"=="+u.length);
    float[][][] g1f = new float[n3][n2][ng1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        for (int ig=0; ig<ng1; ig++)
          g1f[i3][i2][ig] = x1[i3][i2][g1Flat[ig]];
      }
    }
    float[] g2f = new float[ng2];
    for (int ig=0; ig<ng2; ig++)
      g2f[ig] = (float)g2[ig];
    float[] g3f = new float[ng3];
    for (int ig=0; ig<ng3; ig++)
      g3f[ig] = (float)g3[ig];
    
    // compute derivative in first dimension
    float[][][] d = new float[ng3][ng2][ng1];
    for (int i3=0; i3<ng3; i3++) {
      for (int i2=0; i2<ng2; i2++) {
        CubicInterpolator ci = 
            new CubicInterpolator(Method.LINEAR,g1f[g3[i3]][g2[i2]],u[i3][i2]);
        ci.interpolate1(g1f[g3[i3]][g2[i2]],d[i3][i2]);
      }
    }
    
    // interpolate shifts and derivatives in the second dimension
    float[][] u23 = new float[ng3][ng2];
    float[][] d23 = new float[ng3][ng2];
    float[][][] ui23 = new float[n3][n2][ng1];
    float[][][] di23 = new float[n3][n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i3=0; i3<ng3; i3++) {
        for (int i2=0; i2<ng2; i2++) {
          u23[i3][i2] = u[i3][i2][i1];
          d23[i3][i2] = d[i3][i2][i1];
        }
      }
      BicubicInterpolator2 ciu = new BicubicInterpolator2(
          BicubicInterpolator2.Method.SPLINE,
          BicubicInterpolator2.Method.SPLINE,
          g2f,g3f,u23);
      BicubicInterpolator2 cid = new BicubicInterpolator2(
          BicubicInterpolator2.Method.SPLINE,
          BicubicInterpolator2.Method.SPLINE,
          g2f,g3f,d23);
      for (int i3=0; i3<n3; i3++) {
        for (int i2=0; i2<n2; i2++) {
          ui23[i3][i2][i1] = ciu.interpolate(i2,i3);
          di23[i3][i2][i1] = cid.interpolate(i2,i3);
        }  
      }
    }
    
    // interpolate in the first dimension, using derivatives
    float[][][] ui = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        CubicInterpolator ci = 
            new CubicInterpolator(g1f[i3][i2],ui23[i3][i2],di23[i3][i2]);
        for (int i1=0; i1<n1; i1++) {
          ui[i3][i2][i1] = ci.interpolate(i1);
        }
      }
    }
    return ui;
  }
  
  public static float[][] interpolateSparseShifts(
      int n1, int n2, int[] g1, int[] g2, float[][] u, 
      BicubicInterpolator2.Method m)
  {
    int ng1 = g1.length;
    int ng2 = g2.length;
    Check.argument(ng1==u[0].length,"ng1==u[0].length: "+ng1+"=="+u[0].length);
    Check.argument(ng2==u.length,"ng2==u.length: "+ng2+"=="+u.length);
    float[] g1f = new float[ng1];
    float[] g2f = new float[ng2];
    for (int ig=0; ig<ng1; ig++) {
      g1f[ig] = (float)g1[ig];
    }
    for (int ig=0; ig<ng2; ig++) {
      g2f[ig] = (float)g2[ig];
    }
    
    Sampling s1 = new Sampling(n1,1.0,0.0);
    Sampling s2 = new Sampling(n2,1.0,0.0);
    BilinearInterpolator2 bli = new BilinearInterpolator2(g1f,g2f,u);
    return bli.interpolate(s1,s2);
//    BicubicInterpolator2 bci = new BicubicInterpolator2(m,m,g1f,g2f,u);
//    return bci.interpolate(s1,s2);
  }
  
  public float[][] interpolateSparseShiftsFlat(
      int n1, int n2, int[] g1Flat, int[] g2, float[][] u)
  {
    int ng1 = g1Flat.length;
    int ng2 = g2.length;
    Check.argument(ng1==u[0].length,"ng1==u[0].length: "+ng1+"=="+u[0].length);
    Check.argument(ng2==u.length,"ng2==u.length: "+ng2+"=="+u.length);
    float[] g1f = new float[ng1];
    for (int ig=0; ig<ng1; ig++)
      g1f[ig] = (float)g1Flat[ig];
    float[] g2f = new float[ng2];
    for (int ig=0; ig<ng2; ig++)
      g2f[ig] = (float)g2[ig];

    Viewer vu = new Viewer(u);
    vu.setTitle("Shifts Before Interpolation");
    vu.setVLabel("Time (sparse samples)");
    vu.setHLabel("Crossline (sparse samples)");
    vu.setColorModel1(ColorMap.JET);
    vu.addColorBar("shift");
    vu.setSize(600,900);
    vu.show();
    Sampling s1 = new Sampling(n1,1.0,0.0);
    Sampling s2 = new Sampling(n2,1.0,0.0);
    BilinearInterpolator2 bli = new BilinearInterpolator2(g1f,g2f,u);
    float[][] uiFlat = bli.interpolate(s1,s2);
    Viewer vuif = new Viewer(uiFlat);
    vuif.setTitle("Interpolated Shifts Flattened");
    vuif.setVLabel("Time (samples)");
    vuif.setHLabel("Crossline (samples)");
    vuif.setColorModel1(ColorMap.JET);
    vuif.addColorBar("shift");
    vuif.setSize(600,900);
    vuif.show();
    return uiFlat;
  }
  
  public float[][][] interpolateSparseShiftsFlat(
      int n1, int n2, int n3, int[] g1Flat, int[] g2, int[] g3, float[][][] u)
  {
    int ng1 = g1Flat.length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    Check.argument(
        ng1==u[0][0].length,"ng1==u[0][0].length: "+ng1+"=="+u[0][0].length);
    Check.argument(ng2==u[0].length,"ng2==u[0].length: "+ng2+"=="+u[0].length);
    Check.argument(ng3==u.length,"ng3==u.length: "+ng2+"=="+u.length);
    float[] g1f = new float[ng1];
    for (int ig=0; ig<ng1; ig++)
      g1f[ig] = (float)g1Flat[ig];
    float[] g2f = new float[ng2];
    for (int ig=0; ig<ng2; ig++)
      g2f[ig] = (float)g2[ig];
    float[] g3f = new float[ng3];
    for (int ig=0; ig<ng3; ig++)
      g3f[ig] = (float)g3[ig];
    
    Sampling s1 = new Sampling(n1,1.0,0.0);
    Sampling s2 = new Sampling(n2,1.0,0.0);
    Sampling s3 = new Sampling(n3,1.0,0.0);
    TrilinearInterpolator3 tli = new TrilinearInterpolator3(g1f,g2f,g3f,u);
    float[][][] uiFlat = tli.interpolate(s1,s2,s3);
    Viewer vuif = new Viewer(uiFlat);
    vuif.setTitle("Interpolated Shifts Flattened");
    vuif.setVLabel("Time (samples)");
    vuif.setHLabel("Crossline (samples)");
    vuif.setColorModel1(ColorMap.JET);
    vuif.addColorBar("shift");
    vuif.setClips1(100f,400f);
    vuif.setSize(600,900);
    vuif.show();
    return uiFlat;

//    float[][] ui = new float[n2][n1];
//    float[] r = new float[n1];
//    ramp(0.0f,1.0f,r);
//    for (int i2=0; i2<n2; i2++) {
//      float[] u1 = new float[n1];
//      for (int i1=0; i1<n1; i1++)
//        u1[i1] = fm.u1[i2][i1];
//      CubicInterpolator ciu = new CubicInterpolator(r,uiFlat[i2]);
//      ui[i2] = ciu.interpolate(u1);
//    }
//    Viewer vui = new Viewer(ui);
//    vui.setTitle("Interpolated Shifts Unflattened");
//    vui.setVLabel("Sparse Time (samples)");
//    vui.setHLabel("Sparse Crossline (samples)");
//    vui.setColorModel1(ColorMap.JET);
//    vui.addColorBar("shift");
//    vui.setClips1(100f,400f);
//    vui.setSize(600,900);
//    vui.show();
//    return ui;
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
  
  public static float[] extrapolateShifts(float[] u, int n1) {
    int nu = u.length;
    float v = u[nu-1];
    float[] eu = new float[n1];
    for (int iu=0; iu<nu; ++iu) {
      eu[iu] = u[iu];
    }
    for (int iu=nu; iu<n1; ++iu) {
      eu[iu] = v;
    }
    return eu;
  }
  
  public static float[][] extrapolateShifts(float[][] u, int n1) {
    int nu = u[0].length;
    int n2 = u.length;
    float[][] eu = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float v = u[i2][nu-1];
      for (int iu=0; iu<nu; ++iu) {
        eu[i2][iu] = u[i2][iu];
      }
      for (int iu=nu; iu<n1; ++iu) {
        eu[i2][iu] = v;
      }
    }
    return eu;
  }
  
  public static float[][][] extrapolateShifts(float[][][] u, int n1) {
    int nu = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    float[][][] eu = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        float v = u[i3][i2][nu-1];
        for (int iu=0; iu<nu; iu++) {
          eu[i3][i2][iu] = u[i3][i2][iu];
        }
        for (int iu=nu; iu<n1; ++iu) {
          eu[i3][i2][iu] = v;
        }
      }  
    }
    return eu;
  }
  
  public static void normalize(float[] f) {
    float min = min(f);
    float max = max(f);
    int n1 = f.length;
    float shift = min;
    float scale = (max>min)?1.0f/(max-min):1.0f;
    for (int i1=0; i1<n1; ++i1) {
      f[i1] = (f[i1]-shift)*scale;
    }
  }
  
  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
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
  private CubicInterpolator.Method _m1 =  
      CubicInterpolator.Method.MONOTONIC; // 1D shift interpolation method
  private BicubicInterpolator2.Method _m2 = 
      BicubicInterpolator2.Method.MONOTONIC; // 2D shift interpolation method
  private TricubicInterpolator3.Method _m3 =
      TricubicInterpolator3.Method.MONOTONIC; // 3D shift interpolation method
  private int _fr; // fractional shift factor (1.0/_fr is the shift interval)
  private int _frMax; // fractional shift max. Controls maximum strain.
  private int _k2Min; // minimum constraint on second derivative of u
  private int _k2Max; // maximum constraint on second derivative of u
  private Sampling _shifts1; // sampling of pp to ps1 shift values
  private Sampling _shiftsS; // sampling of ps1 to ps2 shift values
  private float _scale; // computed scaling for pp traces
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private SincInterp _si; // for warping with non-integer shifts

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
    Check.argument(nps>=npp,"nps>=npp");
    int nppMax = (int)(nps*_scale);
    int nl = nps-nppMax;
    if (maxShift!=0) {
      nl = maxShift;
      nppMax = nps-nl;
    }
    nppMax = (nppMax>npp)?npp:nppMax;
//    int ppN1Max = npp-nl+1;
    nl = (nl-1)*_fr+1;
    System.out.println("nppMax="+nppMax+", nl1="+nl);
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
  
  private void computeErrors(float[] f, float[] g, float[][] e) {
    int n1max = e.length;
    int nl1 = e[0].length;
    for (int i1=0; i1<n1max; ++i1) {
      for (int il1=0,j1=i1*_fr; il1<nl1; ++il1,++j1) {
        e[i1][il1] = error(f[i1],g[j1]);
      }
    }
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
  
  private void computeErrorsFrac(
      float[] pp, float[] ps1, float[] ps2, float[][][] e)
  {
    int n1 = pp.length;
    int n1max = e.length;
    int nl1 = e[0].length;
    int nlS = e[0][0].length;
    int nx = (n1-1)*_fr+1;
    float[] ps2i = new float[nx];
    _si.interpolate(n1,1.0,0.0,ps2,nx,_shiftsS.getDelta(),0.0,ps2i);

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

  private static void smoothSparseErrors(
      float[][][] es, float r1min, float r1max, 
      int[][] g1, int[] g2, float r2max) 
  {
    int ng1 = g1[0].length;
    int ng2 = g2.length;
    int nl = es[0][0].length;
    float[][] e2 = new float[ng2][];
    for (int i1=0; i1<ng1; i1++) {
      for (int i2=0; i2<ng2; i2++)
        e2[i2] = es[i2][i1];
      float[][] ef = like(e2);
      float[][] er = like(e2);
      float[][] m = like(e2);
      accumulate2( 1,e2,ef,m,g2,r2max);
      accumulate2(-1,e2,er,m,g2,r2max);
      float scale = 1.0f/ng2;
      for (int i2=0; i2<ng2; i2++) {
        for (int il=0; il<nl; il++) {
          e2[i2][il] = scale*(ef[i2][il]+er[i2][il]-e2[i2][il]);   
        }
      }
      for (int i2=0; i2<ng2; i2++)
        es[i2][i1] = e2[i2];
    }
    for (int i2=0; i2<ng2; i2++) {
      float[][] ef = like(es[i2]);
      float[][] er = like(es[i2]);
      float[][] m = like(es[i2]);
      accumulate( 1,es[i2],ef,m,g1[i2],r1min,r1max);
      accumulate(-1,es[i2],er,m,g1[i2],r1min,r1max);
      float scale = 1.0f/ng1;
      for (int i1=0; i1<ng1; i1++) {
        for (int il=0; il<nl; il++) {
          es[i2][i1][il] = scale*(ef[i1][il]+er[i1][il]-es[i2][i1][il]);   
        }
      }
    }
  }
  
  private static void smoothSparseErrors(
      float[][][][] es, float r1min, float r1max, 
      int[][][] g1, int[] g2, int[] g3, float r2max,float r3max) 
  {
    int ng1 = g1[0][0].length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    int nl = es[0][0][0].length;
    float[][] e3 = new float[ng3][];
    for (int i1=0; i1<ng1; i1++) {
      for (int i2=0; i2<ng2; i2++) {
        for (int i3=0; i3<ng3; i3++) {
          e3[i3] = es[i3][i2][i1];
        }
        float[][] ef = like(e3);
        float[][] er = like(e3);
        float[][] m = like(e3);
        accumulate2( 1,e3,ef,m,g3,r3max);
        accumulate2(-1,e3,er,m,g3,r3max);
        float scale = 1.0f/ng3;
        for (int i3=0; i3<ng3; i3++) {
          for (int il=0; il<nl; il++) {
            e3[i3][il] = scale*(ef[i3][il]+er[i3][il]-e3[i3][il]);
          }
        }
        for (int i3=0; i3<ng3; i3++)
          es[i3][i2][i1] = e3[i3];
      }
    }
    for (int i3=0; i3<ng3; i3++) {
      smoothSparseErrors(es[i3],r1min,r1max,g1[i3],g2,r2max);  
    }
  }
  
  /**
   * Non-linear accumulation of alignment errors constraining the path
   * to increase monotonically with increasing time, or decrease 
   * monotonically with decreasing time depending on the accumulation 
   * direction. 
   * <p>
   * This method uses a simple stencil that only checks
   * the previous sample at the current lag index and one lag index
   * down or up (depeding on dir).
   * @param dir accumulation direction, positive or negative.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private static void accumulate(int dir, float[][] e, float[][] d) {
    int nl = e[0].length;
    int nlm1 = nl-1;
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    int ii = ib;
    for (int il=0; il<nl; ++il)
      d[ii][il] = e[ii][il];
    ii += is;
    for (; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is)); // next trace index from ii 
      for (int il=0; il<nl; ++il) {
        int ilc = il+ic;
        ilc = (ilc==-1)?0:(ilc==nl)?nlm1:ilc; // index of lag constraint
        float dc = d[ji][ilc];
        float di = d[ji][il ];
        d[ii][il] = min(dc,di)+e[ii][il];
      }
    }
  }
  
  /**
   * Non-linear accumulation of alignment errors constraining the path
   * to increase monotonically with increasing time, or decrease 
   * monotonically with decreasing time depending on the accumulation 
   * direction. 
   * <p>
   * This method uses a stencil that checks the previous sample at 
   * a range of lag indices. This range is from 0, the current lag index,
   * to ic*_frMax where ic is -1 if accumulating forward or +1 if
   * accumulating reverse.
   * @param dir accumulation direction, positive or negative.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   * @param m moves array[ni][nl] of recorded moves used for backtracking.
   */
  private void accumulate(int dir, float[][] e, float[][] d, float[][] m) {
    int nl = e[0].length;
    int nlm1 = nl-1;
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    int ii=ib;
    for (int il=0; il<nl; ++il)
      d[ii][il] = e[ii][il];
    ii+=is;
    int kmin = ic*_frMax;
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      for (int il=0; il<nl; ++il) {
        int k = kmin;
        while ((il+k)<0) k++;
        while ((il+k)>nlm1) k--;
        float dm = d[ji][il];
        int mi = 0;
        for (; k!=0; k+=-ic) {
          float dc = d[ji][il+k];
          if (dc<dm) {
            dm = dc;
            mi = k;
          }
        }
        m[ji][il] = mi;
        d[ii][il] = dm+e[ii][il];
      }
    }
  }
  
  /**
   * Non-linear accumulation of alignment errors constraining the path
   * to increase monotonically with increasing time, or decrease 
   * monotonically with decreasing time depending on the accumulation 
   * direction. 
   * <p>
   * This method uses a stencil that checks the previous sample at 
   * a range of lag indices. This range is from 0, the current lag index,
   * to ic*_frMax where ic is -1 if accumulating forward or +1 if
   * accumulating reverse.
   * @param dir accumulation direction, positive or negative.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   * @param m moves array[ni][nl] of recorded moves used for backtracking.
   */
  private static void accumulate(
      int dir, float[][] e, float[][] d, float[][] m, int kmax) {
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
    for (int il=0; il<nl; ++il) {
      d[ii][il] = e[ii][il];
    }
    ii+=is;
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      for (int il=0; il<nl; ++il) {
        int ke = kmax;
        while ((il+ke*ic)<0) ke+=ic;
        while ((il+ke*ic)>nlm1) ke-=ic;
        float dm = dmax;
        int mi = 0;
        for (int k=0; k<=ke; k++) {
          float dc = d[ji][il+k*ic];
          if (dc<dm) {
            dm = dc;
            mi = k*ic;
          }
        }
        m[ji][il] = mi;
        d[ii][il] = dm+e[ii][il];
      }
    }
  }
  
  private static void accumulate(
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
  
  private static void accumulate2(
      int dir, float[][] e, float[][] d, float[][] m, int[] g, float rmax) 
  {
    int nl = e[0].length;
    int nlm1 = nl-1;
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    float dmax = ni;
    int ii=ib;
    for (int il=0; il<nl; ++il) {
      d[ii][il] = e[ii][il];
    }
    ii+=is;
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      int iemax = g[ii];
      int iemin = g[ji];
      int dg = abs(iemax-iemin); // sparse grid delta
      // Loop over all lags.
      for (int il=0; il<nl; ++il) {
        int kmax = (int)floor(rmax*dg);
        int kmin = -kmax;
        while ((il+kmin)<0) kmin++;
        while ((il+kmax)>nlm1) kmax--;
        float dm = dmax; // Initialize minimum accumulation error.
        int mi = 0; // Initalize move index
        for (int k=kmin; k<=kmax; k++) {
          float dc = d[ji][il+k];
          if (dc<dm) {
            dm = dc;
            mi = k;
          }
        }
        m[ji][il] = mi;
        d[ii][il] = dm+e[ii][il];
      }
    }
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
  
  //TODO remove
  private static void accumulateSparse(
      int dir, float rmin, float rmax, int[] g, Map<Integer,float[][]> mxbMap,
      float[][] e, float[][] d, float[][] m)
  {
    int nl   = e[0].length;
    int nlm1 = nl-1;
    int ng   = g.length;
    int ngm1 = ng-1;
    int ib = (dir>0)?0:ngm1; // beginning index
    int ie = (dir>0)?ng:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    int isp = ib; // sparse grid index
    int ier = g[isp]; // error index
    float dmax = Float.MAX_VALUE; // default accumulation max value
    
    // Initialize accumulation values
    for (int il=0; il<nl; ++il)
      d[isp][il] = e[ier][il];
    isp += is;

    // Loop over all sparse grid points.
    for (; isp!=ie; isp+=is) {
      int ispm = isp-is; // previous sparse grid index
      ier = g[isp]; // new error index
      int iermin = g[ispm]; // min error index for interpolation.
      int dg = abs(ier-iermin); // sparse grid delta
      float[][] dgmxb = mxbMap.get(dg);
      int kmin = (int)ceil( rmin*dg);
      int kmax = (int)floor(rmax*dg);
//      System.out.println("ei="+ier+", dg="+dg+", kmin="+kmin+", kmax="+kmax);
      
      // Loop over all lags.
      for (int il=0; il<nl; ++il) {
        float dm = dmax; // Initialize minimum accumulation error.
        int   mi = 0;    // Initalize move index
        
        // Loop over all possible slopes, interpolating the alignment
        // errors between the sparse grid points.
        for (int k=kmin; k<=kmax; k++) {
          int rk = k*ic;
          int ik = il+rk;
          if (ik<0 || ik>nlm1)
            continue;
//          float dc = interpSlope(rk,is,dg,il,iermin,ier,d[ispm][ik],e);
          float dc = interpSlope(rk,is,il,iermin,ier,dgmxb[k],d[ispm][ik],e);
//          float dc = interpSlopeNN(rk,is,il,iermin,ier,dgmxb[k],d[ispm][ik],e);
          if (dc<dm) {
            dm = dc;
            mi = rk;
          }
        }
        m[ispm][il] = mi;
        d[isp ][il] = dm+e[ier][il];
      }
    }
  }
  
  private static void accumulateSparseNew(
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
  
  private static float interpSlope(
      int k, int is, int dg, int il, int iemin, int iemax, 
      float dp, float[][] e) 
  {
    float dc = dp;
    if (k==0) {
      for (int x=iemin+is; x!=iemax; x+=is) {
        dc += e[x][il]; 
      }
    } else {
      double slope = (double)-k*is/dg;
      for (int x=iemin+is; x!=iemax; x+=is) {
        double y = il+(slope*(x-iemin)+k);
        int y1 = (int)y;
        int y2 = y1+1;
//        dc += (y2-y)*e[x][y1]+(y-y1)*e[x][y2];
        dc += e[x][y1]+(y-y1)*(e[x][y2]-e[x][y1]);
      }
    }
    return dc;
  }
  
  private static float interpSlope(
      int k, int is, int il, int iemin, int iemax, 
      float[] mxb, float dp, float[][] e) 
  {
    float dc = dp;
    if (k==0) {
      for (int x=iemin+is; x!=iemax; x+=is) {
        dc += e[x][il]; 
      }
    } else {
      int n = mxb.length;
      int nm = n-1;
      for (int x=1; x<nm; x++) {
        int ex = iemin+is*x;
        double y = il+mxb[x];
        int y1 = (int)y;
        int y2 = y1+1;
//        dc += (y2-y)*e[ex][y1]+(y-y1)*e[ex][y2];
        dc += e[ex][y1]+(y-y1)*(e[ex][y2]-e[ex][y1]);
      }
    }
    return dc;
  }
  
  private static float interpSlopeNN(
      int k, int is, int dg, int il, int iemin, int iemax, 
      float dp, float[][] e) 
  {
    float dc = dp;
    if (k==0) {
      for (int x=iemin+is; x!=iemax; x+=is) {
        dc += e[x][il]; 
      }
    } else {
      double slope = (double)-k*is/dg;
      for (int x=iemin+is; x!=iemax; x+=is) {
        int y = (int)(il+(slope*(x-iemin)+k)+0.5);
        dc += e[x][y];
      }
    }
    return dc;
  }
  
  private static float interpSlopeNN(
      int k, int is, int il, int iemin, int iemax, 
      float[] mxb, float dp, float[][] e)
  {
    float dc = dp;
    if (k==0) {
      for (int x=iemin+is; x!=iemax; x+=is) {
        dc += e[x][il]; 
      }
    } else {
      int n = mxb.length;
      int nm = n-1;
      for (int x=1; x<nm; x++) {
        int ex = iemin+is*x;
        int y = (int)(il+mxb[x]+0.5);
        dc += e[ex][y];
      }
    }
    return dc;
  }
  
  //TODO remove
  private static void accumulateSparse2(
      int dir, float rmax, int[] g, float[][] e, float[][] d, float[][] m)
  {
    int n1   = e.length;
    int nl   = e[0].length;
    int nlm1 = nl-1;
    int ng   = g.length;
    int ngm1 = ng-1;
    int ib = (dir>0)?0:ngm1; // beginning index
    int ie = (dir>0)?ng:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int isp = ib; // sparse grid index
    int ier = g[isp]; // error index
    float dmax = n1; // default accumulation max value
    
    // Initialize accumulation values
    for (int il=0; il<nl; ++il) {
      d[isp][il] = e[ier][il];
    }
    isp += is;      

    // Loop over all sparse grid points.
    for (; isp!=ie; isp+=is) {
      int ispm = isp-is; // previous sparse grid index
      ier = g[isp]; // new error index
      int iermin = g[ispm]; // min error index for interpolation.
      int dg = abs(ier-g[ispm]); // sparse grid delta
      
      // Loop over all lags.
      for (int il=0; il<nl; ++il) {
        int kmax = (int)floor(rmax*dg);
        int kmin = -kmax;
        while ((il+kmin)<0) kmin++;
        while ((il+kmax)>nlm1) kmax--;
        float dm = dmax; // Initialize minimum accumulation error.
        int mi = 0; // Initalize move index
        
        // Loop over all possible slopes, interpolating the alignment
        // errors between the sparse grid points.
        for (int k=kmin; k<=kmax; k++) {
          float dc = interpSlope(k,is,dg,il,iermin,ier,d[ispm][il+k],e);
//          float dc = interpSlopeNN(k,is,dg,il,iermin,ier,d[ispm][il+k],e);
          if (dc<dm) {
            dm = dc;
            mi = k;
          }
        }
        m[ispm][il] = mi;
        d[isp ][il] = dm+e[ier][il];
      }
    }
  }
  
  /**
   * Finds shifts by backtracking in accumulated alignment errors.
   * Backtracking must be performed in the direction opposite to
   * that for which accumulation was performed.
   * @param dir backtrack direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param shifts a Sampling of the shift values that starts at
   *  shiftMin, has length nl and has an interval of 1.0/fr
   * @param d input array[ni][nl] of accumulated errors.
   * @param e input array[ni][nl] of alignment errors.
   * @param u output array[ni] of computed shifts.
   */
  private static void backtrack(
      int dir, int b, Sampling shifts, float[][] d, float[][] e, float[] u) 
    {
      float ob = 1.0f/b;
      int nl = d[0].length;
      int ni = d.length;
      int nlm1 = nl-1;
      int nim1 = ni-1;
      int ib = (dir>0)?0:nim1; // begining index
      int ie = (dir>0)?nim1:0; // end index
      int is = (dir>0)?1:-1;   // stride
      int ic = (dir>0)?1:-1;   // constraint dir, forward=+lag, reverse=-lag
      int ii = ib;
      // Set initial lag for the case that all errors at ii are equal.
      int il = (dir>0)?0:nlm1;
//      int il = 0;
      float dl = d[ii][il]; // Current accumulated error value.
          
//    Find minimum lag value(dl) and index(il) at trace index ii.
      for (int jl=0; jl<nl; ++jl) {
        if (d[ii][jl]<dl) {
          dl = d[ii][jl];
          il = jl;
        }
      }
      u[ii] = (float)shifts.getValue(il);
        
      // Iterate through traces, finding the minimum path that satisfies
      // the constraints.
      while (ii!=ie) {
        int ji = max(0,min(nim1,ii+is));   // next trace index from ii
        int jb = max(0,min(nim1,ii+is*b)); // ji scaled by b, for strain limits 
        int ilc = il+ic;
        ilc = (ilc==-1)?0:(ilc==nl)?nlm1:ilc; // index of lag constraint
        float dc = d[jb][ilc];
        float di = d[ji][il ];
        for (int kb=ji; kb!=jb; kb+=is) {
          dc += e[kb][ilc];
        }
        dl = min(dc,di);
        if (dl!=di)
          il = ilc;
        ii += is;
        u[ii] = (float)shifts.getValue(il);
          
        // If the minimum path is changing lags, set the shift value.
        // This value can be non-integer if ob<1 and jb!=ji 
        if (il==ilc) {
          float du = (u[ii]-u[ii-is])*ob;
          u[ii] = u[ii-is]+du;
          for (int kb=ji; kb!=jb; kb+=is) {
            ii += is;
            u[ii] = u[ii-is]+du;
          }
        }
      }
    }
  
//  private void backtrack(
//      int dir, Sampling shifts, float[][] d, float[] u) 
//    {
//      int nl = d[0].length;
//      int ni = d.length;
//      int nlm1 = nl-1;
//      int nim1 = ni-1;
//      int ib = (dir>0)?0:nim1; // begining index
//      int ie = (dir>0)?nim1:0; // end index
//      int is = (dir>0)?1:-1;   // stride
//      int ic = (dir>0)?1:-1;   // constraint dir, forward=+lag, reverse=-lag
//      int ii = ib;
//      // Set initial lag for the case that all errors at ii are equal.
//      int il = (dir>0)?0:nlm1; 
//      float dl = d[ii][il]; // Current accumulated error value.
//          
//      // Find minimum lag value(dl) and index(il) at trace index ii.
//      for (int jl=0; jl<nl; ++jl) {
//        if (d[ii][jl]<dl) {
//          dl = d[ii][jl];
//          il = jl;
//        }
//      }
//      u[ii] = (float)shifts.getValue(il);
//        
//      // Iterate through traces, finding the minimum path that satisfies
//      // the constraints.
//      while (ii!=ie) {
//        int ji = max(0,min(nim1,ii+is));   // next trace index from ii
//        float dm = d[ji][il];
//        int ilc = il+ic*_frMax;
//        ilc = (ilc<0)?0:(ilc>nlm1)?nlm1:ilc; // index of lag constraint
//        int ilf = 0;
//        while (ilf<_frMax) {
//          float dc = d[ji][ilc];
//          if (dc<dm) {
//            dm = dc;
//            il = ilc;
//          }
//          ilc+=-ic;
//          if (ilc<0 || ilc>nlm1)
//            break;
//          ilf++;
//        }
//        ii += is;
//        u[ii] = (float)shifts.getValue(il);
//      }
//    }
  
  private void backtrack(
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
  
  private void backtrack(
      int dir, Sampling shifts, float[][] m, float[] u) 
    {
      int nl = m[0].length;
      int ni = m.length;
      int nlm1 = nl-1;
      int nim1 = ni-1;
      int ib = (dir>0)?0:nim1; // begining index
      int ie = (dir>0)?nim1:0; // end index
      int is = (dir>0)?1:-1;   // stride
      int ii = ib;
      // Set initial lag for the case that all errors at ii are equal.
      int il = (dir>0)?0:nlm1; 
      float dl = m[ii][il]; // Current accumulated error value.
          
      // Find minimum lag value(dl) and index(il) at trace index ii.
      for (int jl=0; jl<nl; ++jl) {
        if (m[ii][jl]<dl) {
          dl = m[ii][jl];
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
    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
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
    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
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
    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
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

  private static float[] getEnvelope(float[] x) {
    HilbertTransformFilter htf = new HilbertTransformFilter();
    int n = x.length;
    float[] y = new float[n];
    float[] e = new float[n];
    htf.apply(n,x,y);
    for (int i=0; i<n; i++)
      e[i] = sqrt(x[i]*x[i]+y[i]*y[i]);
    return e;
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

  private static float[] like(float[] a) {
    return new float[a.length];
  }
  private static float[][] like(float[][] a) {
    return new float[a.length][a[0].length];
  }
  private static float[][][] like(float[][][] a) {
    return new float[a.length][a[0].length][a[0][0].length];
  }

  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }
  
}
