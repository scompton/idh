package warp;

import edu.mines.jtk.dsp.RecursiveExponentialFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * 
 */
public class DynamicWarpingG {

  private int _ng;
  private int _glMax;

  public DynamicWarpingG(
      float gAvg, float gMin, float gMax, float gDelta, float dgMax) {
    Check.argument(gAvg>1,"gAvg>1");
    Check.argument(gMin>=1,"gMin>=1");
    Check.argument(gDelta>0,"gDelta>0");
    Check.argument(dgMax>=gDelta,"dgMax>=gDelta");
    _scale = 2.0f/(gAvg+1.0f);
    _ng = 1 + (int)((gMax-gMin)/gDelta);
    _fr = (int)(1.0f/(gDelta/2.0f));
    _glMax = (int)(1.0f/dgMax);
    _si1 = new SincInterpolator();
  }
 
  /**
   * Sets extent of smoothing filters used to smooth integer shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factor. Default 
   * factor is zero, for no smoothing.
   * @param usmooth extent of smoothing filter in all dimensions.
   */
  public void setShiftSmoothing(double usmooth) {
    setShiftSmoothing(usmooth,usmooth);
  }

  /**
   * Sets extents of smoothing filters used to smooth integer shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   */
  public void setShiftSmoothing(double usmooth1, double usmooth2) {
    setShiftSmoothing(usmooth1,usmooth2,usmooth2);
  }

  /**
   * Sets extents of smoothing filters used to smooth integer shifts.
   * Half-widths of smoothing filters are inversely proportional to
   * strain limits, and are scaled by the specified factors. Default 
   * factors are zero, for no smoothing.
   * @param usmooth1 extent of smoothing filter in 1st dimension.
   * @param usmooth2 extent of smoothing filter in 2nd dimension.
   * @param usmooth3 extent of smoothing filter in 3rd dimension.
   */
  public void setShiftSmoothing(
    double usmooth1, double usmooth2, double usmooth3) 
  {
    _usmooth1 = usmooth1;
    _usmooth2 = usmooth2;
    _usmooth3 = usmooth3;
    updateSmoothingFilters();
  }

  public float[] findShifts(float[] f, float[] g) {
    float[][] e = computeErrors(f,g);
    float[][] d = accumulateForward(e);
    float[] u = backtrackReverse(d,e);
    smoothShifts(u,u);
    return u;
  }
  
  public float[] findShifts1(float[][] f, float[][] g) {
    float[][] e = computeErrors1(f,g);
    float[][] d = accumulateForward(e);
    float[] u = backtrackReverse(d,e);
    smoothShifts(u,u);
    return u;
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
  
  public static float[] vpvs(float[] u, double sigma) {
    int nu = u.length;
    float[] vpvs = new float[nu];
    vpvs[0] = 1.0f;
    for (int iu=1; iu<nu; ++iu) {
      vpvs[iu] = 1.0f + 2*(u[iu]-u[iu-1]);
    }
    if (sigma > 0.0) {
      RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
      ref.apply(vpvs,vpvs);
    }
    return vpvs;
  }
  
  public static float[][] vpvs(float[][] u, double sigma) {
    int nu = u[0].length;
    int n2 = u.length;
    float[][] vpvs = new float[n2][nu];
    for (int i2=0; i2<n2; ++i2) {
      vpvs[i2][0] = 1.0f;
      for (int iu=1; iu<nu; ++iu) {
        vpvs[i2][iu] = 1.0f+2*(u[i2][iu]-u[i2][iu-1]);
      }
    }
    if (sigma > 0.0) {
      RecursiveExponentialFilter ref = new RecursiveExponentialFilter(sigma);
      ref.apply(vpvs,vpvs);
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

  /**
   * Returns smoothed shifts.
   * @param u array of shifts to be smoothed.
   * @return array of smoothed shifts
   */
  public float[] smoothShifts(float[] u) {
    float[] us = like(u);
    smoothShifts(u,us);
    return us;
  }

  /**
   * Returns smoothed shifts.
   * @param u array of shifts to be smoothed.
   * @return array of smoothed shifts
   */
  public float[][] smoothShifts(float[][] u) {
    float[][] us = like(u);
    smoothShifts(u,us);
    return us;
  }

  /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothShifts(float[] u, float[] us) {
    if (_ref1!=null) {
      _ref1.apply(u,us); 
    } else if (u!=us) {
      copy(u,us);
    }
  }

  /**
   * Smooths the specified shifts. Smoothing can be performed 
   * in place; input and output arrays can be the same array.
   * @param u input array of shifts to be smoothed.
   * @param us output array of smoothed shifts.
   */
  public void smoothShifts(float[][] u, float[][] us) {
    if (_ref1!=null) {
      _ref1.apply1(u,us);
    } else {
      copy(u,us);
    }
    if (_ref2!=null)
      _ref2.apply2(us,us);
  }

  public float[][] accumulateForward(float[][] e) {
    float[][] d = like(e);
    accumulateForward(e,d);
    return d;
  }
  
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,e,d);
  }

  public float[] backtrackReverse(float[][] d, float[][] e) {
    float[] g = new float[d.length];
    backtrackReverse(d,e,g);
    return g;
  }
 
  public void backtrackReverse(float[][] d, float[][] e, float[] g) {
    backtrack(-1,_glMax,_shifts1,d,e,g);
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

//  private int _nl1; // number of pp to ps1 lags
//  private int _nlS; // number of ps1 to ps2 lags 
  private int _fr; // fractional shift factor (1.0/_fr is the shift interval)
  private Sampling _shifts1; // sampling of pp to ps1 shift values
  private Sampling _shiftsS; // sampling of ps1 to ps2 shift values
  private float _scale; // computed scaling for pp traces
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private double _usmooth1 = 0.0; // extent of smoothing shifts in 1st dim
  private double _usmooth2 = 0.0; // extent of smoothing shifts in 2nd dim
  private double _usmooth3 = 0.0; // extent of smoothing shifts in 3rd dim
  private RecursiveExponentialFilter _ref1; // for smoothing shifts
  private RecursiveExponentialFilter _ref2; // for smoothing shifts
  private RecursiveExponentialFilter _ref3; // for smoothing shifts
  private SincInterpolator _si1; // for warping with non-integer shifts
  private SincInterpolator _siS; // for warping with non-integer shifts

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

  private void updateSmoothingFilters() {
    _ref1 = (_usmooth1<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth1);
    _ref2 = (_usmooth2<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth2);
    _ref3 = (_usmooth3<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth3);
  }
  
  private int[] computeErrorLengths(int n1, int maxShift) {
    int nl1;
    if (maxShift==0)
      nl1 = n1-(int)(n1*_scale);
    else
      nl1 = maxShift;
    int ppN1Max = n1-nl1+1;
    nl1 = (nl1-1)*_fr+1;
    System.out.println("ppN1Max="+ppN1Max+", nl1="+nl1);
    _shifts1 = new Sampling(nl1,1.0/_fr,0.0);
    return new int[]{ppN1Max,nl1};
  }
  
  private int[] computeErrorLengths(int n1, float fsp) {
    int nl1 = n1-(int)(n1*_scale);
    int nlS = (int)(nl1*fsp);
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
  
  /**
   * Non-linear accumulation of alignment errors.
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
    for (int il=0; il<nl; ++il)
      d[ib][il] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));   // next trace index from ii 
      for (int il=0; il<nl; ++il) {
        int ilc = il+ic;
        ilc = (ilc==-1)?0:(ilc==nl)?nlm1:ilc; // index of lag constraint
        float dc = d[ji][ilc];
        float di = d[ji][il ];
        d[ii][il] = min(dc,di)+e[ii][il];
      }
    }
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
   * @param g output array[ni] of computed shifts.
   */
  private static void backtrack(
      int dir, int gmax, Sampling shifts, float[][] d, float[][] e, float[] g) 
    {
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
      g[ii] = 1.0f;
      
      // Iterate through traces, finding the minimum path that satisfies
      // the constraints.
      while (ii!=ie) {
        int ji = max(0,min(nim1,ii+is));   // next trace index from ii
        int ilc = il+ic;
        ilc = (ilc==-1)?0:(ilc==nl)?nlm1:ilc; // index of lag constraint
        float dc = d[ji][ilc];
        float di = d[ji][il ];
        int gc = 1;
        int ig = ilc+ic;
        while (gc<gmax) {
          if (ig==-1 || ig==nl)
            break;
          if (d[ji][ig]<dc) {
            dc = d[ji][ig];
            ilc = ig;
          }
          gc++;
          ig += ic;
        }
        dl = min(dc,di);
        int du;
        if (dl!=di) {
          du = il-ilc;
          il = ilc;
        }else {
          du = 0;
        }
        ii += is;
        g[ii] = 1.0f+2.0f*du;
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
