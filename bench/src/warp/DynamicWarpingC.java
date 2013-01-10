package warp;

import edu.mines.jtk.dsp.LinearInterpolator;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.util.*;
import edu.mines.jtk.util.CubicInterpolator.Method;
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
    _si1 = new SincInterpolator();
  }

  public DynamicWarpingC(
      float vpvsAvg, int fr, int frMax, int k2Min, int k2Max)
  {
    Check.argument(vpvsAvg>=1,"vpvsAvg>=1");
  Check.argument(fr>0,"fr>0");
  Check.argument(k2Min<k2Max, "k2Min<k2Max");
  _fr = fr;
  _frMax = frMax;
  _k2Min = k2Min;
  _k2Max = k2Max;
  _scale = 2.0f/(vpvsAvg+1.0f);
  _si1 = new SincInterpolator();
  _siS = new SincInterpolator();
  }
  
  public void setShifts(int n1pp, int n1ps) {
    computeErrorLengths(n1pp,n1ps,0);
  }
  
  /**
   * Sets the exponent used to compute alignment errors |f-g|^e.
   * The default exponent is 2.
   * @param e the exponent.
   */
  public void setErrorExponent(double e) {
    _epow = (float)e;
  }

  /**
   * Sets the number of nonlinear smoothings of alignment errors.
   * In dynamic warping, alignment errors are smoothed the specified 
   * number of times, along all dimensions (in order 1, 2, ...), 
   * before estimating shifts by accumulating and backtracking along 
   * only the 1st dimension. 
   * <p> 
   * The default number of smoothings is zero, which is best for 1D
   * sequences. For 2D and 3D images, two smoothings are recommended.
   * @param esmooth number of nonlinear smoothings.
   */
  public void setErrorSmoothing(int esmooth) {
    _esmooth = esmooth;
  }

//  public float[] findShifts(float[] f, float[] g) {
//    float[][] e = computeErrors(f,g);
//    for (int is=0; is<_esmooth; ++is)
//      smoothErrors(e,e);
//    float[][] d = accumulateForward(e);
//    float[] u = backtrackReverse(d,e);
//    smoothShifts(u,u);
//    return u;
//  }
//  
//  public float[] findShifts1(float[][] f, float[][] g) {
//    float[][] e = computeErrors1(f,g);
//    for (int is=0; is<_esmooth; ++is)
//      smoothErrors(e,e);
//    float[][] d = accumulateForward(e);
//    float[] u = backtrackReverse(d,e);
//    smoothShifts(u,u);
//    return u;
//  }
  
  public float[] findShifts(
      float[] f, float[] g, float rmin, float rmax, float dr) {
    float[][] e = computeErrors(f,g);
    float[] u0 = getU0(e);
    float[][][] sm = accumulateForward(e,u0[0],rmin,rmax,dr);
    return backtrackReverse2(sm[0],sm[1]);
  }
  
  public float[] findShifts(
      float[][] f, float[][] g, float rmin, float rmax, float dr) {
    float[][] e = computeErrors1(f,g);
    float[] u0 = getU0(e);
    float[][][] sm = accumulateForward(e,u0[0],rmin,rmax,dr);
    return backtrackReverse2(sm[0],sm[1]);
  }
  
  public float[] findShifts(
      float[][][] f, float[][][] g, float rmin, float rmax, float dr) {
    float[][] e = computeErrors1(f,g);
    float[] u0 = getU0(e);
    float[][][] sm = accumulateForward(e,u0[0],rmin,rmax,dr);
    return backtrackReverse2(sm[0],sm[1]);
  }
  
  public float[] findShifts(
      float[][] e, float rmin, float rmax, float dr) {
    float[] u0 = getU0(e);
    float[][][] sm = accumulateForward(e,u0[0],rmin,rmax,dr);
    return backtrackReverse2(sm[0],sm[1]);
  }
  
  public float[][] findShifts(
      float[][] e, float rmin, float rmax, float dr, float u0Perc) {
    Check.argument(u0Perc >= 0 && u0Perc <= 100,"u0Perc >= 0 && u0Perc <= 100");
    float[] u0 = getU0(e);
    int nu = (int)(u0.length*u0Perc/100.0f);
    nu = (nu==0)?1:nu;
    float[][] u = new float[nu][];
    System.out.println("Number of initial shifts(u0) to test: "+nu);
    for (int iu=0; iu<nu; iu++) {
      float[][][] sm = accumulateForward(e,u0[iu],rmin,rmax,dr);
      u[iu] = backtrackReverse2(sm[0],sm[1]);
    }
    return u;
  }
  
  public float[] findShifts(float[][] e, int kmin, int kmax) {
    float[][][] dm = accumulateForward(e,kmin,kmax);
    return backtrackReverse(dm[0],dm[1]);
  }
  
  public float[][] findShifts(float[][][] e, int kmin, int kmax) {
    int n2 = e.length;
    final float[][][] fe = e;
    final int fkmin= kmin;
    final int fkmax = kmax;
    final float[][] u = new float[n2][];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(fe[i2],fkmin,fkmax);
      u[i2] = backtrackReverse(dm[0],dm[1]);
    }});
    
    return u;
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
  
  public float[] applyShifts(float[] u, float[] g) {
    int n1 = g.length;
    float[] h = new float[n1]; 
    applyShifts(u,g,h);
    return h;
  }
  
  public void applyShifts(float[] u, float[] g, float[] h) {
    int n1 = g.length;
    int nu = u.length;
    int num = nu-1;
    _si1.setUniformSampling(n1,1.0,0.0);
    _si1.setUniformSamples(g);
    for (int iu=0; iu<nu; ++iu) {
      h[iu] = _si1.interpolate(iu+u[iu]);
    }
    for (int i1=nu; i1<n1; ++i1) {
      h[i1] = _si1.interpolate(i1+u[num]);
    }
  }
  
  public float[][] applyShifts(float[][] u, float[][] g) {
    int n1 = g[0].length;
    int n2 = g.length;
    float[][] h = new float[n2][n1]; 
    applyShifts(u,g,h);
    return h;
  }
  
  public void applyShifts(float[][] u, float[][] g, float[][] h) {
    Check.argument(g.length==u.length,"g.length==u.length");
    final int n1 = g[0].length;
    final int n2 = g.length;
    final int n1u = u[0].length;
    final int n1um = n1u-1;
    final float[][] uf = u;
    final float[][] gf = g;
    final float[][] hf = h;
    final Parallel.Unsafe<SincInterpolator> siu =
      new Parallel.Unsafe<SincInterpolator>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      SincInterpolator si = siu.get();
      if (si==null) {
        si = new SincInterpolator();
        si.setUniformSampling(n1,1.0,0.0);
        siu.set(si);
      }
      si.setUniformSamples(gf[i2]);
      for (int i1=0; i1<n1u; ++i1) {
        hf[i2][i1] = si.interpolate(i1+uf[i2][i1]);
      }
      for (int i1=n1u; i1<n1; ++i1) {
        hf[i2][i1] = si.interpolate(i1+uf[i2][n1um]);
      }
    }});
  }
  
  public float[][] applyShifts(
      float[] u1, float[] uS, float[] ps1, float[] ps2)
  {
    int n1 = ps1.length;
    float[] ps1w = new float[n1]; 
    float[] ps2w = new float[n1];
    int nu = u1.length;
    int num = nu-1;
    _si1.setUniformSampling(n1,1.0,0.0);
    _si1.setUniformSamples(ps1);
    _siS.setUniformSampling(n1,1.0,0.0);
    _siS.setUniformSamples(ps2);
    for (int iu=0; iu<nu; ++iu) {
      ps1w[iu] = _si1.interpolate(iu+u1[iu]);
      ps2w[iu] = _siS.interpolate(iu+u1[iu]+uS[iu]);
    }
    for (int i1=nu; i1<n1; ++i1) {
      ps1w[i1] = _si1.interpolate(i1+u1[num]);
      ps2w[i1] = _siS.interpolate(i1+u1[num]+uS[num]);
    }
    return new float[][]{ps1w, ps2w};
  }
  
  public static float[] vpvs(float[] u) {
    int nu = u.length;
    float[] vpvs = new float[nu];
    float um1 = u[0] - (u[1]-u[0]);
//    System.out.println(um1);
    vpvs[0] = 1.0f + 2*(u[0]-um1);
    for (int iu=1; iu<nu; ++iu) {
      vpvs[iu] = 1.0f + 2*(u[iu]-u[iu-1]);
    }
    return vpvs;
  }
  
  public static float[][] vpvs(float[][] u) {
    int nu = u[0].length;
    int n2 = u.length;
    float[][] vpvs = new float[n2][nu];
    for (int i2=0; i2<n2; ++i2) {
      float um1 = u[i2][0] - (u[i2][1]-u[i2][0]);
      vpvs[i2][0] = 1.0f + 2*(u[i2][0]-um1);
      for (int iu=1; iu<nu; ++iu) {
        vpvs[i2][iu] = 1.0f+2*(u[i2][iu]-u[i2][iu-1]);
      }
    }
    return vpvs;
  }
  
  public static float[][][] vpvs(float[][][] u) {
    int n1 = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    float[][][] vpvs = new float[n3][n2][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        vpvs[i3][i2][0] = 1.0f;
        for (int i1=1; i1<n1; ++i1) {
          vpvs[i3][i2][i1] = 1.0f+2*(u[i3][i2][i1]-u[i3][i2][i1-1]);
        }
      }  
    }
    return vpvs;
  }

  ///////////////////////////////////////////////////////////////////////////
  // for research and atypical applications

  public float[][] computeErrors(float[] pp, float[] ps) {
    int npp = pp.length;
    int nps = ps.length;
    int[] el = computeErrorLengths(npp,nps,0);
    float[][] e = new float[el[0]][el[1]];
    if (_fr>1) {
      int nx = (npp-1)*_fr+1;
      SincInterpolator si = new SincInterpolator();
      si.setUniformSampling(npp,1.0f,0.0f);
      si.setUniformSamples(ps);
      float[] psi = new float[nx];
      si.interpolate(nx,_shifts1.getDelta(),0.0,psi);
      computeErrors(pp,psi,e);
    } else {
      computeErrors(pp,ps,e);
    }
    normalizeErrors(e);
    return e;
  }
  
  public float[][][] computeErrors(float[][] pp, float[][] ps) {
    int npp1 = pp[0].length;
    int npp2 = pp.length;
    int nps1 = ps[0].length;
    int nps2 = ps.length;
    Check.argument(npp2==nps2,"npp2==nps2");
    final float[][] fpp = pp;
    final float[][] fps = ps;
    int[] el = computeErrorLengths(npp1,nps1,0);
    final float[][][] e = new float[npp2][el[0]][el[1]];
    Parallel.loop(npp2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrors(fpp[i2],fps[i2],e[i2]);  
    }});
    normalizeErrors(e);
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
        SincInterpolator si = new SincInterpolator();
        si.setUniformSampling(n1,1.0f,0.0f);
        si.setUniformSamples(fps1[i2]);
        float[] ps1i = new float[nx];
        si.interpolate(nx,_shifts1.getDelta(),0.0,ps1i);
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
        SincInterpolator si = new SincInterpolator();
        si.setUniformSampling(n1pp,1.0f,0.0f);
        si.setUniformSamples(fps1[i3][i2]);
        float[] ps1i = new float[nx];
        si.interpolate(nx,_shifts1.getDelta(),0.0,ps1i);
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
      SincInterpolator si = new SincInterpolator();
      si.setUniformSampling(n1,1.0f,0.0f);
      si.setUniformSamples(ps2);
      float[] ps2i = new float[nx];
      si.interpolate(nx,_shifts1.getDelta(),0.0,ps2i);
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
      SincInterpolator si = new SincInterpolator();
      si.setUniformSampling(n1,1.0f,0.0f);
      si.setUniformSamples(ps2);
      float[] ps2i = new float[nx];
      si.interpolate(nx,_shiftsS.getDelta(),0.0,ps2i);
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
          SincInterpolator si = new SincInterpolator();
          si.setUniformSampling(n1,1.0f,0.0f);
          si.setUniformSamples(fps2[i2]);
          float[] ps2i = new float[nx];
          si.interpolate(nx,_shifts1.getDelta(),0.0,ps2i);
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
        SincInterpolator si = new SincInterpolator();
        si.setUniformSampling(n1,1.0f,0.0f);
        si.setUniformSamples(fps2[i2]);
        float[] ps2i = new float[nx];
        si.interpolate(nx,_shifts1.getDelta(),0.0,ps2i);
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

  public float[] getU0(float[][] e) {
    float[][] d = like(e);
    float[][] m = like(e);
    accumulate(-1,e,d,m);
    int nl = d[0].length;
    float [] u0 = new float[nl];
    int[] i = rampint(0,1,nl);
    quickIndexSort(d[0],i);
    for (int il=0; il<nl; il++) {
      u0[il] = (float)_shifts1.getValue(i[il]);
//      System.out.println("index i="+i[il]+", d="+d[0][i[il]]+", u0="+u0[il]);
    }
    return u0;
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
  
  public static float[][] smoothErrorsSparse(
      float[][] e, float rmin, float rmax, float dr, boolean extrapolate)
  {
    int n1 = e.length;
    int dg = (int)ceil(1.0f/dr);
    int nl = e[0].length;
    int[] g = extrapolate?
        Decompose.decomposeUniform(n1,dg):Decompose.decompose(n1,dg);
    int ng = g.length;        
    float[][] ef = new float[ng][nl];
    float[][] er = new float[ng][nl];
    float[][] es = new float[ng][nl];
    float[][] m  = new float[ng][nl];
    if (extrapolate) {
      float[][] e2 = Decompose.extrapolateErrors(dg,e);
      accumulateSparse( 1,rmin,rmax,g,e2,ef,m);
      accumulateSparse(-1,rmin,rmax,g,e2,er,m);
      for (int i1=0; i1<ng; i1++) {
        for (int il=0; il<nl; il++) {
          es[i1][il] = ef[i1][il]+er[i1][il]-e2[g[i1]][il];
        }
      }
    } else {
      accumulateSparse( 1,rmin,rmax,g,e,ef,m);  
      accumulateSparse(-1,rmin,rmax,g,e,er,m);
      for (int i1=0; i1<ng; i1++) {
        for (int il=0; il<nl; il++) {
          es[i1][il] = ef[i1][il]+er[i1][il]-e[g[i1]][il];
        }
      }
    }
    normalizeErrors(es);
    return es;
  }
  
  public static float[][][] smoothErrorsSparse(
      float[][][] e, 
      final float r1min, final float r1max, float dr1,
      final float r2min, final float r2max, float dr2)
  {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    int dg1 = (int)ceil(1.0f/dr1);
    int dg2 = (int)ceil(1.0f/dr2);
    
    final int[] g1 = Decompose.decomposeUniform(n1,dg1);    
    dump(g1);
    final int ng1 = g1.length;
    final float[][][] es1 = new float[n2][ng1][nl];
    final float[][][] ee1 = Decompose.extrapolateErrors1(dg1,e);
    Parallel.loop(n2, new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] m  = new float[ng1][nl];
      float[][] ef1 = new float[ng1][nl];
      float[][] er1 = new float[ng1][nl];
      accumulateSparse( 1,r1min,r1max,g1,ee1[i2],ef1,m);
      accumulateSparse(-1,r1min,r1max,g1,ee1[i2],er1,m);
      for (int i1=0; i1<ng1; i1++) {
        for (int il=0; il<nl; il++) {
          es1[i2][i1][il] = ef1[i1][il]+er1[i1][il]-ee1[i2][g1[i1]][il];
        }
      }        
    }});
    normalizeErrors(es1);
    
    final int[] g2 = Decompose.decomposeUniform(n2,dg2);
    dump(g2);
    final int ng2 = g2.length;
    final float[][][] es = new float[ng2][ng1][nl];
    final float[][][] ee2 = Decompose.extrapolateErrors2(dg2,es1);
    final int n2e = ee2.length;
    final Parallel.Unsafe<float[][][]> eeu = 
      new Parallel.Unsafe<float[][][]>();
    Parallel.loop(ng1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][][] ee = eeu.get();
      if (ee==null) eeu.set(ee=new float[4][ng2][nl]);
      float[][] ee1 = new float[n2e][nl];
      float[][]   m = ee[0];
      float[][] es2 = ee[1];
      float[][] ef2 = ee[2];
      float[][] er2 = ee[3];
      for (int i2=0; i2<ng2; ++i2) {
        es2[i2] = es[i2][i1];
      }
      for (int i2=0; i2<n2e; ++i2) {
        ee1[i2] = ee2[i2][i1];
      }
      accumulateSparse2( 1,r2min,r2max,g2,ee1,ef2,m);
      accumulateSparse2(-1,r2min,r2max,g2,ee1,er2,m);
      for (int i2=0; i2<ng2; ++i2) {
        for (int il=0; il<nl; ++il) {
          es2[i2][il] = ef2[i2][il]+er2[i2][il]-ee1[i2][il];
        }
      }
    }});
    normalizeErrors(es);
    return es;
  }
  
  public static float[][][][] smoothErrorsSparse(
      float[][][][] e, 
      final float r1min, final float r1max, float dr1,
      final float r2min, final float r2max, float dr2,
      final float r3min, final float r3max, float dr3)
  {
    final int n3 = e.length;
    float[][][][] e3 = new float[n3][][][];
    for (int i3=0; i3<n3; i3++) {
      e3[i3] = smoothErrorsSparse(e[i3],r1min,r1max,dr1,r2min,r2max,dr2);
    }
    final int nl  = e3[0][0][0].length;
    final int ng1 = e3[0][0].length;
    final int ng2 = e3[0].length;
    int dg3 = (int)ceil(1.0f/dr3);
    final int[] g3 = Decompose.decomposeUniform(n3,dg3);
    dump(g3);
    final int ng3 = g3.length;
    final float[][][][] es = new float[ng3][ng2][ng1][nl];
    final float[][][][] ee3 = Decompose.extrapolateErrors2(dg3,e3);
    final int n3e = ee3.length;
    final Parallel.Unsafe<float[][][]> eeu = 
      new Parallel.Unsafe<float[][][]>();
    Parallel.loop(ng1,new Parallel.LoopInt() {
    public void compute(int i1) {
      for (int i2=0; i2<ng2; ++i2) {
        float[][][] ee = eeu.get();
        if (ee==null) eeu.set(ee=new float[4][ng3][nl]);
        float[][] ee1 = new float[n3e][nl];
        float[][]   m = ee[0];
        float[][] es3 = ee[1];
        float[][] ef3 = ee[2];
        float[][] er3 = ee[3];
        for (int i3=0; i3<ng3; ++i3) {
          es3[i3] = es[i3][i2][i1];  
        }
        for (int i3=0; i3<n3e; ++i3) {
          ee1[i3] = ee3[i3][i2][i1];  
        }
        accumulateSparse2( 1,r2min,r2max,g3,ee1,ef3,m);
        accumulateSparse2(-1,r2min,r2max,g3,ee1,er3,m);
        for (int i3=0; i3<ng3; ++i3) {
          for (int il=0; il<nl; ++il) {
            es3[i3][il] = ef3[i3][il]+er3[i3][il]-ee1[i3][il];
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
  
  public float[][][] accumulateForward(float[][] e, int kmin, int kmax) {
    float[][] d = like(e);
    float[][] m = like(e);
    accumulate(1,e,d,m,kmin,kmax);
    normalizeErrors(d);
    return new float[][][]{d,m};
  }
  
  public float[][][] accumulateForwardK2(float[][] e) {
    float[][] d = like(e);
    float[][] m = like(e);
    accumulateK2(1,e,d,m);
    return new float[][][]{d,m};
  }
  
  public float[][][] accumulateForward(
      float[][] e, float u0, float rmin, float rmax, float dr) {
    Sampling sr = new Sampling((int)((rmax-rmin)/dr+1),dr,rmin);
    int n1 = e.length;
    int nr = sr.getCount();
    float[][] moves = new float[n1][nr];
    float[][] s = new float[n1][nr];
    accumulateR(1,e,s,moves,u0,sr);
    return new float[][][]{s,moves};
  }
  
  public float[][][] accumulateReverse(
      float[][] e, float u0, float rmin, float rmax, float dr) {
    Sampling sr = new Sampling((int)((rmax-rmin)/dr+1),dr,rmin);
    int n1 = e.length;
    int nr = sr.getCount();
    float[][] moves = new float[n1][nr];
    float[][] s = new float[n1][nr];
    accumulateR(-1,e,s,moves,u0,sr);
    return new float[][][]{s,moves};
  }
  
  public static int[] getSparseGrid(int n, float dr, boolean extrapolate) {
    int delta = (int)ceil(1.0f/dr);
    int[] g = extrapolate?
        Decompose.decomposeUniform(n,delta):Decompose.decompose(n,delta);
    return g;
  }
  
  public float[][][] accumulateForwardSparse(
      float[][] e, float rmin, float rmax, float dr, boolean extrapolate)
  {
    int n1 = e.length;
    int dg = (int)ceil(1.0f/dr);
//  int kmin = (int)ceil(rmin/dr);
//  int kmax = (int)floor(rmax/dr);
    int[] g = extrapolate?
        Decompose.decomposeUniform(n1,dg):Decompose.decompose(n1,dg);
    int ng = g.length;
    int nl = e[0].length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    if (extrapolate) {
      float[][] e2 = Decompose.extrapolateErrors(dg,e);
      accumulateSparse(1,rmin,rmax,g,e2,d,m);
    } else {
      accumulateSparse(1,rmin,rmax,g,e,d,m);  
    }
    normalizeErrors(d);
    return new float[][][]{d,m};
  }
  
  public float[][][] accumulateForwardSparse2(
      float[][] e, float rmin, float rmax, float dr, boolean extrapolate)
  {
    int n1 = e.length;
    int dg = (int)ceil(1.0f/dr);
//    int kmin = (int)ceil(rmin/dr);
//    int kmax = (int)floor(rmax/dr);
    int[] g = extrapolate?
        Decompose.decomposeUniform(n1,dg):Decompose.decompose(n1,dg);
    int ng = g.length;
    int nl = e[0].length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    if (extrapolate) {
      float[][] e2 = Decompose.extrapolateErrors(dg,e);
      accumulateSparse2(1,rmin,rmax,g,e2,d,m);
    } else {
      accumulateSparse2(1,rmin,rmax,g,e,d,m);  
    }
    normalizeErrors(d);
    return new float[][][]{d,m};
  }
  
  public float[][][] accumulateReverseSparse(
      float[][] e, float rmin, float rmax, float dr, boolean extrapolate)
  {
    int n1 = e.length;
    int dg = (int)ceil(1.0f/dr);
//    int kmin = (int)ceil(rmin/dr);
//    int kmax = (int)floor(rmax/dr);
    int[] g = extrapolate?
        Decompose.decomposeUniform(n1,dg):Decompose.decompose(n1,dg);
    int ng = g.length;
    int nl = e[0].length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    if (extrapolate) {
      float[][] e2 = Decompose.extrapolateErrors(dg,e);
      accumulateSparse(-1,rmin,rmax,g,e2,d,m);
    } else {
      accumulateSparse(-1,rmin,rmax,g,e,d,m);  
    }
    normalizeErrors(d);
    return new float[][][]{d,m};
  }
  
  public float[][][] accumulateReverseSparse2(
      float[][] e, float rmin, float rmax, float dr, boolean extrapolate)
  {
    int n1 = e.length;
    int dg = (int)ceil(1.0f/dr);
//    int kmin = (int)ceil(rmin/dr);
//    int kmax = (int)floor(rmax/dr);
    int[] g = extrapolate?
        Decompose.decomposeUniform(n1,dg):Decompose.decompose(n1,dg);
    int ng = g.length;
    int nl = e[0].length;
    System.out.println("dg="+dg+", ng="+ng);
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    if (extrapolate) {
      float[][] e2 = Decompose.extrapolateErrors(dg,e);
      accumulateSparse2(-1,rmin,rmax,g,e2,d,m);
    } else {
      accumulateSparse2(-1,rmin,rmax,g,e,d,m);  
    }
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

  public static float[] interpolateSparseShifts(int n, int[] g, float[] u) {
    int ng = g.length;
    float[] gf = new float[ng];
    for (int ig=0; ig<ng; ig++) {
      gf[ig] = (float)g[ig];
    }
    CubicInterpolator ci = new CubicInterpolator(Method.MONOTONIC,ng,gf,u);
//    CubicInterpolator ci = new CubicInterpolator(Method.LINEAR,ng,gf,u);
    float[] ui = new float[n];
    for (int i=0; i<n; i++) {
      ui[i] = ci.interpolate(i);
    }
    return ui;
  }

  public static float[][] interpolateSparseShifts(
      int n1, int n2, int[] g1, int[] g2, float[][] u)
  {
    final int n1f = n1;
    final int n2f = n2;
    final int ng1 = g1.length;
    final int ng2 = g2.length;
    Check.argument(ng1==u[0].length,"ng1==u[0].length");
    Check.argument(ng2==u.length,"ng2==u.length");
    final float[][] uf = u;
    final float[] gf1 = new float[ng1];
    final float[] gf2 = new float[ng2];
    for (int ig=0; ig<ng1; ig++) {
      gf1[ig] = (float)g1[ig];
    }
    for (int ig=0; ig<ng2; ig++) {
      gf2[ig] = (float)g2[ig];
    }
    
    final float[][] ui1 = new float[ng2][n1];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      CubicInterpolator ci1 = 
          new CubicInterpolator(Method.MONOTONIC,ng1,gf1,uf[i2]);
      for (int i1=0; i1<n1f; i1++) {
        ui1[i2][i1] = ci1.interpolate(i1);
      }
    }});
    
    final float[][] ui2 = new float[n2][n1];
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[] u1 = new float[ng2];
      for (int i2=0; i2<ng2; i2++) {
        u1[i2] = ui1[i2][i1];
      }
      CubicInterpolator ci2 = 
          new CubicInterpolator(Method.MONOTONIC,ng2,gf2,u1);
      for (int i2=0; i2<n2f; i2++) {
        ui2[i2][i1] = ci2.interpolate(i2);
      } 
    }});
    return ui2;
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
  
  public static float[][][] transposeLag01(float[][][] e) {
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
  
  public static float[][][][] transposeLag01(float[][][][] e) {
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

  private float _vpvsMin;
  private int _fr; // fractional shift factor (1.0/_fr is the shift interval)
  private int _frMax; // fractional shift max. Controls maximum strain.
  private int _k2Min; // minimum constraint on second derivative of u
  private int _k2Max; // maximum constraint on second derivative of u
  private Sampling _shifts1; // sampling of pp to ps1 shift values
  private Sampling _shiftsS; // sampling of ps1 to ps2 shift values
  private float _scale; // computed scaling for pp traces
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private int _esmooth = 0; // number of nonlinear smoothings of errors
  private double _usmooth1 = 0.0; // extent of smoothing shifts in 1st dim
  private double _usmooth2 = 0.0; // extent of smoothing shifts in 2nd dim
  private double _usmooth3 = 0.0; // extent of smoothing shifts in 3rd dim
  private SincInterpolator _si1; // for warping with non-integer shifts
  private SincInterpolator _siS; // for warping with non-integer shifts

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
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0f,0.0f);
    si.setUniformSamples(ps2);
    float[] ps2i = new float[nx];
    si.interpolate(nx,_shiftsS.getDelta(),0.0,ps2i);

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
    for (int il=0; il<nl; ++il) {
      d[ii][il] = e[ii][il];
    }
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
      int dir, float[][] e, float[][] d, float[][] m, int kmin, int kmax) {
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
        int ks = kmin;
        int ke = kmax;
//        while ((il+ke*ic)<0) ke+=ic;
//        while ((il+ke*ic)>nlm1) ke-=ic;
        if (dir>0) {
          if ((il+ks*ic)<0) {
            m[ji][il] = 0;
            d[ii][il] = dmax+e[ii][il];
            continue;
          }
          while ((il+ke*ic)<0) ke+=ic;
        } else {
          if ((il+ks*ic)>nlm1) {
            m[ji][il] = 0;
            d[ii][il] = dmax+e[ii][il];
            continue;
          }
          while ((il+ke*ic)>nlm1) ke-=ic;
        }
        
        float dm = dmax;
        int mi = 0;
        for (int k=ks; k<=ke; k++) {
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
  
  @Deprecated
  /**
   * Non-linear accumulation of alignment errors constraining the path
   * to increase monotonically with increasing time, or decrease 
   * monotonically with decreasing time depending on the accumulation 
   * direction. 
   * <p>
   * This method uses a stencil that checks the previous sample at 
   * a range of lag indices. This range is from 0, the current lag index,
   * to ic*_frMax where ic is -1 if accumulating forward or +1 if
   * accumulating reverse. However, in order for lag indices to be checked
   * they must satisfy the _k2Min and _k2Max constraints, that control how
   * quickly shift are allowed to change from sample to sample.
   * @param dir accumulation direction, positive or negative.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   * @param m moves array[ni][nl] of recorded moves used for backtracking.
   */
  private void accumulateK2(int dir, float[][] e, float[][] d, float[][] m) {
    int nl = e[0].length;
    int nlm1 = nl-1;
    int ni = e.length;
    int nim1 = ni-1;
    float dmax = ni;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    int fcount = 0;
    int ii=ib;
    for (int il=0; il<nl; ++il) {
      d[ii][il] = e[ii][il];
      m[ii][il] = 0;
    }
    ii+=is;
    int kmin = ic*_frMax;
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      for (int il=0; il<nl; ++il) {
        int k = kmin;
        while ((il+k)<0) k++;
        while ((il+k)>nlm1) k--;
        float dm = dmax;
        int mi = 0;
        for (; k!=-ic; k+=-ic) {
          float ds = k-m[ji][il+k];
          if (ds>=_k2Min && ds<=_k2Max) {
            float dc = d[ji][il+k];
            if (dc<dm) {
              dm = dc;
              mi = k;
            }
          }
        }
        if (dm==dmax) {
          fcount++;
        }
        m[ii][il] = mi;
        d[ii][il] = dm+e[ii][il];
      }
    }
    System.out.println("Number of accumulation failures = "+fcount);
  }
  
  @Deprecated
  /**
   * Accumulation using shift increments.
   * @param dir accumulation direction, positive or negative.
   * @param e input array[ni][nl] of alignment errors.
   * @param s 
   * @param moves moves array[ni][nl] of recorded moves used for backtracking.
   * @param u0
   * @param sr
   */
  private void accumulateR(int dir, float[][] e, float[][] s,
      float[][] moves,float u0, Sampling sr)
  {
    int n1 = e.length;
    int n1m1 = n1-1;
    int nl = e[0].length;
    int umax = nl/_fr;
    int nr = sr.getCount();
    int nrm1 = nr-1;
    float[][] d = new float[2][nr];
    SincInterpolator si = new SincInterpolator();
    si.setUniform(nl,0.5,0.0,e[0]);
    int ib = (dir>0)?0:n1m1; // beginning index
    int ie = (dir>0)?n1:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int i1 = ib;
    for (int il=0; il<nr; il++) {
      s[i1][il] = u0;
      d[i1][il] = si.interpolate(u0);
    }
    i1 += is;
    for (; i1!=ie; i1+=is) {
      si.setUniform(nl,0.5,0.0,e[i1]);
      float[] d0 = d[0]; d[0] = d[1]; d[1] = d0;
      for (int ir=0; ir<nr; ir++) {
        float r = (float)sr.getValue(ir);
        int j1 = max(0,min(n1m1,i1-is));
        int m = ir;
        int iqm = max(0,min(nrm1,ir-1));
        int iqp = min(nrm1,ir+1);
        float dm = d[1][iqm]+si.interpolate(s[j1][iqm]+r);
        float di = d[1][ir ]+si.interpolate(s[j1][ir ]+r);
        float dp = d[1][iqp]+si.interpolate(s[j1][iqp]+r);
        if (dm<di && dm<dp)
          m = iqm;
        if (dp<di && dp<dm)
          m = iqp;
        float ui = s[j1][m] + r;
        ui = (ui>umax)?umax:ui;
        s[i1][ir] = ui;
        d[0][ir] = d[1][m] + si.interpolate(s[i1][ir]);
        moves[j1][ir] = m-ir;
      }
    }
    moves[n1m1] = d[0];
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
    int n1   = e.length;
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
    float dmax = n1; // default accumulation max value
    
    // Setup interpolation for alignment errors
    LinearInterpolator li = new LinearInterpolator();
    li.setUniform(nl,1.0,0.0,n1,1.0,0.0,e);
    
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
        int kmin = (int)ceil(rmin*dg);
        int kmax = (int)floor(rmax*dg);
        if (dir>0) {
          if ((il+kmin*ic)<0) {
            m[ispm][il] = 0;
            d[isp ][il] = dmax+e[ier][il];
            continue;
          }
//          while ((il+kmin*ic)<0) kmin+=ic;
          while ((il+kmax*ic)<0) kmax+=ic;
        } else {
          if ((il+kmin*ic)>nlm1) {
            m[ispm][il] = 0;
            d[isp ][il] = dmax+e[ier][il];
            continue;
          }
//          while ((il+kmin*ic)>nlm1) kmin-=ic;
          while ((il+kmax*ic)>nlm1) kmax-=ic;
        }
        float dm = dmax; // Initialize minimum accumulation error.
        int mi = 0; // Initalize move index
        
        // Loop over all possible slopes, interpolating the alignment
        // errors between the sparse grid points.
        for (int k=kmin; k<=kmax; k++) {
          int rk = k*ic;
          float slope = (float)k/dg;
          float dc = d[ispm][il+rk];
          for (int j=iermin+is; j!=ier; j+=is) {
            double x1 = il+(slope*(j-iermin)+rk); 
            double x2 = j;
            dc += li.interpolate(x1,x2);
          }
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
  
  private static void accumulateSparse2(
      int dir, float rmin, float rmax, int[] g,
      float[][] e, float[][] d, float[][] m)
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
    
    // Setup interpolation for alignment errors
    LinearInterpolator li = new LinearInterpolator();
    li.setUniform(nl,1.0,0.0,n1,1.0,0.0,e);
    
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
        int kminU = (int)ceil(rmin*dg);
        int kmaxU = (int)floor(rmax*dg);
        int kminD = -kminU;
        int kmaxD = -kmaxU;
        float dm = dmax; // Initialize minimum accumulation error.
        int mi = 0; // Initalize move index
        
        // Loop over all possible slopes, interpolating the alignment
        // errors between the sparse grid points.
        for (int k=kmaxD; k<=kminD; k++) {
          if ((il+k)<0)
            continue;
          float slope = (float)k/dg;
          float dc = d[ispm][il+k];
          for (int j=iermin+is; j!=ier; j+=is) {
            double x1 = il+(slope*(j-iermin)+k);
            double x2 = j;
            dc += li.interpolate(x1,x2);
          }
          if (dc<dm) {
            dm = dc;
            mi = k;
          }
        }
        for (int k=kminU; k<=kmaxU; k++) {
          if ((il+k)>nlm1)
            continue;
          float slope = (float)k/dg;
          float dc = d[ispm][il+k];
          for (int j=iermin+is; j!=ier; j+=is) {
            double x1 = il+(slope*(j-iermin)+k);
            double x2 = j;
            dc += li.interpolate(x1,x2);
          }
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
