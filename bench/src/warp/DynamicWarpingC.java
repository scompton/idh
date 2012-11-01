package warp;

import edu.mines.jtk.dsp.RecursiveExponentialFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * 
 */
public class DynamicWarpingC {

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
  
  public void setShifts(int n1) {
    computeErrorLengths(n1,0);
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
    vpvs[0] = 1.0f;
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
      vpvs[i2][0] = 1.0f;
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

  public float[][] computeErrors(float[] pp, float[] ps1) {
    int n1 = pp.length;
    int[] el = computeErrorLengths(n1,0);
    float[][] e = new float[el[0]][el[1]];
    if (_fr>1) {
      int nx = (n1-1)*_fr+1;
      SincInterpolator si = new SincInterpolator();
      si.setUniformSampling(n1,1.0f,0.0f);
      si.setUniformSamples(ps1);
      float[] ps1i = new float[nx];
      si.interpolate(nx,_shifts1.getDelta(),0.0,ps1i);
      computeErrors(pp,ps1i,e);
    } else {
      computeErrors(pp,ps1,e);
    }
    normalizeErrors(e);
    return e;
  }
  
  public float[][] computeErrors1(float[][] pp, float[][] ps1) {
    final int n1 = pp[0].length;
    int[] el = computeErrorLengths(n1,0);
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
  
  public float[][] computeErrors1(float[][][] pp, float[][][] ps1) {
    final int n1 = pp[0][0].length;
    final int n2 = pp[0].length;
    final int n3 = pp.length;
    int[] el = computeErrorLengths(n1,0);
    final int n1M = el[0]; // maximum pp time index
    final int nl1 = el[1]; // number of pp-ps1 lags
    final float[][][] fpp = pp;
    final float[][][] fps1 = ps1;
    float[][] e = Parallel.reduce(n2*n3,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i23) {
      int i2 = i23%n2;
      int i3 = i23/n2;
      float[][] e = new float[n1M][nl1];
      if (_fr>1) {
        int nx = (n1-1)*_fr+1;
        SincInterpolator si = new SincInterpolator();
        si.setUniformSampling(n1,1.0f,0.0f);
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
    float[][] d1 = like(e);
    float[][] d2 = like(e);
    accumulate( 1, e,d1);
    accumulate(-1,d1,d2);
    normalizeErrors(d2);
    return d2;
  }
  
  public float[][][] accumulateForward(float[][] e) {
    float[][] d = like(e);
    float[][] m = like(e);
    accumulate(1,e,d,m);
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
  
//  public float[][][] accumulateForward(
//      float[][] e, float u0, float rmin, float rmax, float dr) {
//    int n1 = e.length;
//    int n1m1 = n1-1;
//    int nl = e[0].length;
//    int umax = nl/_fr;
//    Sampling sr = new Sampling((int)((rmax-rmin)/dr+1),dr,rmin);
//    int nr = sr.getCount();
//    int nrm1 = nr-1;
//    float[][] d = new float[2][nr];
//    float[][] s = new float[n1][nr];
//    float[][] moves = new float[n1][nr];
//    SincInterpolator si = new SincInterpolator();
//    si.setUniform(nl,0.5,0.0,e[0]);
//    for (int il=0; il<nr; il++) {
//      s[0][il] = u0;
//      d[0][il] = si.interpolate(u0);
//    }
//    for (int i1=1; i1<n1; i1++) {
//      si.setUniform(nl,0.5,0.0,e[i1]);
//      float[] d0 = d[0]; d[0] = d[1]; d[1] = d0;
//      for (int ir=0; ir<nr; ir++) {
//        float r = (float)sr.getValue(ir);
//        int j1 = max(0,min(n1m1,i1-1));
//        int m = ir;
//        int iqm = max(0,min(nrm1,ir-1));
//        int iqp = min(nrm1,ir+1);
//        float dm = d[1][iqm]+si.interpolate(s[j1][iqm]+r);
//        float di = d[1][ir ]+si.interpolate(s[j1][ir ]+r);
//        float dp = d[1][iqp]+si.interpolate(s[j1][iqp]+r);
//        if (dm<di && dm<dp)
//          m = iqm;
//        if (dp<di && dp<dm)
//          m = iqp;
//        float ui = s[j1][m] + r;
//        ui = (ui>umax)?umax:ui;
//        s[i1][ir] = ui;
//        d[0][ir] = d[1][m] + si.interpolate(s[i1][ir]);
//        moves[j1][ir] = m-ir;
//      }
//    }
//    moves[n1m1] = d[0];
//    return new float[][][]{s,moves};
//  }
  
  public float[] backtrackReverse(float[][] m) {
    float[] u = new float[m.length];
    backtrack(-1,_shifts1,m,u);
    return u;
  }
  
  public float[] backtrackReverse(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrackReverse(d,m,u);
    return u;
  }
  
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

  public void backtrackReverse(float[][] d, float[][] m, float[] u) {
    backtrack(-1,_shifts1,d,m,u);
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
   * @param n1
   * @param maxShift
   * @return
   */
  private int[] computeErrorLengths(int n1, int maxShift) {
    int nl1;
    if (maxShift==0)
      nl1 = n1-(int)(n1*_scale);
    else
      nl1 = maxShift;
    int ppN1Max = n1-nl1+1;
    nl1 = (nl1-1)*_fr+1;
//    System.out.println("ppN1Max="+ppN1Max+", nl1="+nl1);
    _shifts1 = new Sampling(nl1,1.0/_fr,0.0);
    return new int[]{ppN1Max,nl1};
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
//    System.out.println("ppN1Max="+ppN1Max+", nl1="+nl1+", nlS="+nlS);
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
      m[ii][il] = 0;
    }
    ii+=is;
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      for (int il=0; il<nl; ++il) {
        int k=ic*_frMax;
        while ((il+k)<0) k++;
        while ((il+k)>nlm1) k--;
        float dm = d[ji][il];
        int mi = 0;
        for (; k!=-ic; k+=-ic) {
          float dc = d[ji][il+k];
          if (dc<dm) {
            dm = dc;
            mi = k;
          }
        }
        m[ii][il] = mi;
        d[ii][il] = dm+e[ii][il];
      }
    }
  }
  
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
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      for (int il=0; il<nl; ++il) {
        int k=ic*_frMax;
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
  
  private void accumulateR(
      int dir, float[][] e, float[][] s, float[][] moves,
      float u0, Sampling sr) {
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
      float dl = d[ii][il]; // Current accumulated error value.
          
      // Find minimum lag value(dl) and index(il) at trace index ii.
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
        il += (int)m[ii][il];
        ii += is;
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
