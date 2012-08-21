package warp;

import edu.mines.jtk.dsp.RecursiveExponentialFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterpolator;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic warping of sequences and images.
 * <p>
 * For sequences f and g, dynamic warping finds a sequence of 
 * integer shifts u such that f[i1] ~ g[i1+u[i1]], subject to a 
 * bound b1 on strain, the rate at which the shifts u[i1] vary 
 * with sample index i1.
 * <p>
 * An increasing u[i1] = u[i1-1] + 1 implies that, between indices
 * i1-1 and i1, g[i1] is a stretched version of f[i1] ~ g[i1+u[i1]].
 * For, in this case, values in f for indices i1 and i1-1 are one 
 * sample apart, but corresponding values in g are two samples 
 * apart, which implies stretching by 100%. Likewise, a decreasing 
 * u[i1] = u[i1-1] - 1 implies squeezing by 100%.
 * <p>
 * In practice, 100% strain (stretching or squeezing) may be extreme.
 * Therefore, the upper bound on strain may be smaller than one. For 
 * example, if the bound b1 = 0.5, then the local average strain is 
 * bounded by 0.5. This constraint is complicated by the fact that 
 * the shifts u[i1] are integers. The actual constraint for the bound 
 * b1 = 0.5 is |u[i1-1]-u[i1-2]| + |u[i1]-u[i1-1]| &le; 1.
 * <p>
 * For 2D images f and g, dynamic warping finds a 2D array of integer 
 * shifts u[i2][i1] such that f[i2][i1] ~ g[i2][i1+u[i2][i1]], 
 * subject to bounds b1 and b2 on strains, the rates at which shifts 
 * u[i2][i1] vary with samples indices i1 and i2, respectively.
 * <p>
 * For 3D images f and g, dynamic warping finds a 3D array of integer 
 * shifts u[i3][i2][i1] in a similar way. However, finding shifts for 
 * 3D images may require an excessive amount of memory. Dynamic image 
 * warping requires a temporary array of nlag*nsample floats, where 
 * the number of lags nlag = 1+shiftMax-shiftMin and nsample is the 
 * number of image samples. For 3D images, the product nlag*nsample 
 * is likely to be too large for the temporary array to fit in random-
 * access memory (RAM). In this case, shifts u are obtained by blending 
 * together shifts computed from overlapping subsets of the 3D image.
 * <p>
 * Estimated shifts u are integers, but can be smoothed to obtain
 * non-integer shifts. The extent of smoothing along each dimension 
 * is inversely proportional to the strain limit for that dimension, 
 * and these extents can be scaled by specified factors for more or 
 * less smoothing. The default scale factors are zero, for no 
 * smoothing.
 * <p>
 * This class provides numerous methods, but typical applications
 * require only several of these, usually only the methods that find
 * and apply shifts. The many other methods are provided only for 
 * atypical applications and research.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.07.05
 */
public class DynamicWarpingC {

  
  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   */
  public DynamicWarpingC(int shiftMin, int shiftMax) {
    this(shiftMin,shiftMax,1);
  }
  
  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * Shifts vary at the interval given by 1.0/fr.
   * @param shiftMin lower bound on shift u.
   * @param shiftMax upper bound on shift u.
   * @param fr factor that determines the fractional interval
   *  of shifts.
   */
  public DynamicWarpingC(int shiftMin, int shiftMax, int fr) {
    Check.argument(shiftMin>=0,"shiftMin>=0");
    Check.argument(shiftMax-shiftMin>1,"shiftMax-shiftMin>1");
    Check.argument(fr>0,"fr>0");
    _l1min = shiftMin;
    _l1max = shiftMax;
    _fr = fr;
    _nl1 = (_l1max-_l1min)*_fr+1;
    _shifts1 = new Sampling(_nl1,1.0/_fr,_l1min);
    _si = new SincInterpolator();
  }
  
  public DynamicWarpingC(
      int shift1Min, int shift1Max, int shiftSMin, int shiftSMax, int fr)
  {
    Check.argument(shift1Min>=0,"shift1Min>=0");
    Check.argument(shift1Max-shift1Min>1,"shift1Max-shift1Min>1");
    Check.argument(shiftSMin>=0,"shiftSMin>=0");
    Check.argument(shiftSMax-shiftSMin>1,"shiftSMax-shiftSMin>1");
    Check.argument(fr>0,"fr>0");
    _l1min = shift1Min;
    _l1max = shift1Max;
    _lSmin = shiftSMin;
    _lSmax = shiftSMax;
    _fr = fr;
    _nl1 = 1+_l1max-_l1min;
    _nlS = (_lSmax-_lSmin)*_fr+1;
    _shifts1 = new Sampling(_nl1,1.0,_l1min);
    _shiftsS = new Sampling(_nlS,1.0/_fr,_lSmin);
    _si = new SincInterpolator();
  }

  /**
   * Sets bound on strain for all dimensions. Must be in (0,1].
   * The actual bound on strain is 1.0/ceil(1.0/strainMax), which
   * is less than the specified strainMax when 1.0/strainMax is not
   * an integer. The default bound on strain is 1.0 (100%).
   * @param strainMax the bound, a value less than or equal to one.
   */
  public void setStrainMax(double strainMax) {
    Check.argument(strainMax<=1.0,"strainMax<=1.0");
    Check.argument(strainMax>0.0,"strainMax>0.0");
    setStrainMax(strainMax,strainMax);
  }

  /**
   * Sets bound on strains in 1st and 2nd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   */
  public void setStrainMax(double strainMax1, double strainMax2) {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    Check.argument(strainMax2>0.0,"strainMax2>0.0");
    setStrainMax(strainMax1,strainMax2,strainMax2);
  }

  /**
   * Sets bound on strains in 1st, 2nd and 3rd dimensions.
   * @param strainMax1 bound on strain in the 1st dimension.
   * @param strainMax2 bound on strain in the 2nd dimension.
   * @param strainMax3 bound on strain in the 3rd dimension.
   */
  public void setStrainMax(
    double strainMax1, double strainMax2, double strainMax3)
  {
    Check.argument(strainMax1<=1.0,"strainMax1<=1.0");
    Check.argument(strainMax2<=1.0,"strainMax2<=1.0");
    Check.argument(strainMax3<=1.0,"strainMax3<=1.0");
    Check.argument(strainMax1>0.0,"strainMax1>0.0");
    Check.argument(strainMax2>0.0,"strainMax2>0.0");
    Check.argument(strainMax3>0.0,"strainMax3>0.0");
    _bstrain1 = (int)ceil(1.0/strainMax1);
    _bstrain2 = (int)ceil(1.0/strainMax2);
    _bstrain3 = (int)ceil(1.0/strainMax3);
    updateSmoothingFilters();
  }

  public void setStrainSMax(double strainSMax1) {
    Check.argument(strainSMax1<=1.0,"strainSMax1<=1.0");
    Check.argument(strainSMax1>0.0,"strainSMax1>0.0");
    setStrainSMax(strainSMax1, strainSMax1);
  }
  
  public void setStrainSMax(double strainSMax1, double strainSMax2) {
    Check.argument(strainSMax1<=1.0,"strainSMax1<=1.0");
    Check.argument(strainSMax2<=1.0,"strainSMax2<=1.0");
    Check.argument(strainSMax1>0.0,"strainSMax1>0.0");
    Check.argument(strainSMax2>0.0,"strainSMax2>0.0");
    setStrainSMax(strainSMax1, strainSMax2, strainSMax2);
  }
  
  public void setStrainSMax(
      double strainSMax1, double strainSMax2, double strainSMax3)
  {
    Check.argument(strainSMax1<=1.0,"strainSMax1<=1.0");
    Check.argument(strainSMax2<=1.0,"strainSMax2<=1.0");
    Check.argument(strainSMax3<=1.0,"strainSMax3<=1.0");
    Check.argument(strainSMax1>0.0,"strainSMax1>0.0");
    Check.argument(strainSMax2>0.0,"strainSMax2>0.0");
    Check.argument(strainSMax3>0.0,"strainSMax3>0.0");
    _bstrainS1 = (int)ceil(1.0/strainSMax1);
    _bstrainS2 = (int)ceil(1.0/strainSMax2);
    _bstrainS3 = (int)ceil(1.0/strainSMax3);
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

  /**
   * Sets the size and overlap of windows used for 3D image warping.
   * Window size determines the amount of memory required to store
   * a temporary array[l3][l2][n1][nl] of floats used to compute
   * shifts. Here, nl is the number of lags and n1, l2 and l3 are 
   * the numbers of samples in the 1st, 2nd and 3rd dimensions of 
   * 3D image subsets. Let n1, n2 and n3 denote numbers of samples 
   * for all three dimensions of a complete 3D image. Typically, 
   * l2&lt;n2 and l3&lt;n3, because insufficient memory is available 
   * for a temporary array of nl*n1*n2*n3 floats. 
   * <p>
   * Image subsets overlap in the 2nd and 3rd dimensions by specified 
   * fractions f2 and f3, which must be less than one. Because window 
   * sizes are integers, the actual overlap be greater than (but never 
   * less than) these fractions.
   * <p>
   * Default window sizes are 50 samples; default overlap fractions 
   * are 0.5, which corresponds to 50% overlap in both dimensions.
   * @param l2 length of window in 2nd dimension.
   * @param l3 length of window in 3rd dimension.
   * @param f2 fraction of window overlap in 2nd dimension.
   * @param f3 fraction of window overlap in 3rd dimension.
   */
  public void setWindowSizeAndOverlap(int l2, int l3, double f2, double f3) {
    _owl2 = l2;
    _owl3 = l3;
    _owf2 = f2;
    _owf3 = f3;
  }

  /**
   * Computes and returns shifts for specified sequences.
   * @param f array for the sequence f.
   * @param g array for the sequence g.
   * @return array of shifts u.
   */
  public float[] findShifts(float[] f, float[] g) {
    float[] u = like(f);
    findShifts(f,g,u);
    return u;
  }

  /**
   * Computes and returns shifts for specified images.
   * @param f array for the image f.
   * @param g array for the image g.
   * @return array of shifts u.
   */
  public float[][] findShifts(float[][] f, float[][] g) {
    float[][] u = like(f);
    findShifts(f,g,u);
    return u;
  }

  /**
   * Computes and returns shifts for specified images.
   * @param f array for the image f.
   * @param g array for the image g.
   * @return array of shifts u.
   */
  public float[][][] findShifts(float[][][] f, float[][][] g) {
    float[][][] u = like(f);
    findShifts(f,g,u);
    return u;
  }

  /**
   * Computes and returns 1D shifts u for specified 2D images f and g.
   * This method is useful in the case that shifts vary only slightly 
   * (or perhaps not at all) in the 2nd image dimension.
   * @param f array[n2][n1] for the image f.
   * @param g array[n2][n1] for the image g.
   * @return array[n1] of shifts u.
   */
  public float[] findShifts1(float[][] f, float[][] g) {
    float[] u = like(f[0]);
    findShifts1(f,g,u);
    return u;
  }

  /**
   * Computes and returns 1D shifts u for specified 3D images f and g.
   * This method is useful in the case that shifts vary only slightly 
   * (or perhaps not at all) in the 2nd and 3rd image dimensions.
   * @param f array[n3][n2][n1] for the image f.
   * @param g array[n3][n2][n1] for the image g.
   * @return array[n1] of shifts u.
   */
  public float[] findShifts1(float[][][] f, float[][][] g) {
    float[] u = like(f[0][0]);
    findShifts1(f,g,u);
    return u;
  }

  /**
   * Computes shifts for specified sequences.
   * @param f input array for the sequence f.
   * @param g input array for the sequence g.
   * @param u output array of shifts u.
   */
  public void findShifts(float[] f, float[] g, float[] u) {
    float[][] e = computeErrors(f,g);
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(e,e);
    float[][] d = accumulateForward(e);
    backtrackReverse(d,e,u);
    smoothShifts(u,u);
  }

  /**
   * Computes shifts for specified images.
   * @param f input array for the image f.
   * @param g input array for the image g.
   * @param u output array of shifts u.
   */
  public void findShifts(float[][] f, float[][] g, float[][] u) {
    final float[][][] e = computeErrors(f,g);
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][] uf = u;
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(e,e);
    final Parallel.Unsafe<float[][]> du = new Parallel.Unsafe<float[][]>();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] d = du.get();
      if (d==null) du.set(d=new float[n1][nl]);
      accumulateForward(e[i2],d);
      backtrackReverse(d,e[i2],uf[i2]);
    }});
    smoothShifts(u,u);
  }

  /**
   * Computes shifts for specified images.
   * @param f input array for the image f.
   * @param g input array for the image g.
   * @param u output array of shifts u.
   */
  public void findShifts(float[][][] f, float[][][] g, float[][][] u) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    OverlappingWindows2 ow = 
      new OverlappingWindows2(n2,n3,_owl2,_owl3,_owf2,_owf3);
    int m2 = ow.getM1();
    int m3 = ow.getM2();
    int l2 = ow.getL1();
    int l3 = ow.getL2();
    float[][][] fw = new float[l3][l2][];
    float[][][] gw = new float[l3][l2][];
    float[][][] uw = new float[l3][l2][n1];
    float[][][][] ew = new float[l3][l2][n1][_nl1];
    for (int k3=0; k3<m3; ++k3) {
      int i3 = ow.getI2(k3);
      for (int k2=0; k2<m2; ++k2) {
        int i2 = ow.getI1(k2);
        for (int j3=0; j3<l3; ++j3) {
          for (int j2=0; j2<l2; ++j2) {
            fw[j3][j2] = f[i3+j3][i2+j2];
            gw[j3][j2] = g[i3+j3][i2+j2];
          }
        }
        computeErrors(fw,gw,ew);
        normalizeErrors(ew);
        for (int is=0; is<_esmooth; ++is)
          smoothErrors(ew);
        computeShifts(ew,uw);
        for (int j3=0; j3<l3; ++j3) {
          for (int j2=0; j2<l2; ++j2) {
            float wij = ow.getWeight(i2,i3,j2,j3);
            float[] u32 = u[i3+j3][i2+j2];
            for (int i1=0; i1<n1; ++i1)
              u32[i1] += wij*uw[j3][j2][i1];
          }
        }
      }
    }
    smoothShifts(u);
  }

  /**
   * Computes 1D shifts u for specified 2D images f and g.
   * This method is useful in the case that shifts vary only slightly 
   * (or perhaps not at all) in the 2nd image dimension.
   * @param f input array[n2][n1] for the image f.
   * @param g input array[n2][n1] for the image g.
   * @param u output array[n1] of shifts u.
   */
  public void findShifts1(float[][] f, float[][] g, float[] u) {
    float[][] e = computeErrors1(f,g);
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(e,e);
    float[][] d = accumulateForward(e);
    backtrackReverse(d,e,u);
    smoothShifts(u,u);
  }

  /**
   * Computes 1D shifts u for specified 3D images f and g.
   * This method is useful in the case that shifts vary only slightly 
   * (or perhaps not at all) in the 2nd and 3rd image dimensions.
   * @param f input array[n3][n2][n1] for the image f.
   * @param g input array[n3][n2][n1] for the image g.
   * @param u output array[n1] of shifts u.
   */
  public void findShifts1(float[][][] f, float[][][] g, float[] u) {
    float[][] e = computeErrors1(f,g);
    for (int is=0; is<_esmooth; ++is)
      smoothErrors(e,e);
    float[][] d = accumulateForward(e);
    backtrackReverse(d,e,u);
    smoothShifts(u,u);
  }

  /**
   * Returns a sequence warped by applying specified shifts.
   * @param u array of shifts.
   * @param g array for the sequence to be warped.
   * @return array for the warped sequence.
   */
  public float[] applyShifts(float[] u, float[] g) {
    float[] h = like(g);
    applyShifts(u,g,h);
    return h;
  }

  /**
   * Returns an image warped by applying specified shifts.
   * @param u array of shifts.
   * @param g array for the image to be warped.
   * @return array for the warped image.
   */
  public float[][] applyShifts(float[][] u, float[][] g) {
    float[][] h = like(g);
    applyShifts(u,g,h);
    return h;
  }

  /**
   * Returns an image warped by applying specified shifts.
   * @param u array of shifts.
   * @param g array for the image to be warped.
   * @return array for the warped image.
   */
  public float[][][] applyShifts(float[][][] u, float[][][] g) {
    float[][][] h = like(g);
    applyShifts(u,g,h);
    return h;
  }

  /**
   * Computes a sequence warped by applying specified shifts.
   * @param u input array of shifts.
   * @param g input array for the sequence to be warped.
   * @param h output array for the warped sequence.
   */
  public void applyShifts(float[] u, float[] g, float[] h) {
    int n1 = u.length;
    _si.setUniformSampling(n1,1.0,0.0);
    _si.setUniformSamples(g);
    for (int i1=0; i1<n1; ++i1) {
      h[i1] = _si.interpolate(i1+u[i1]);
    }
  }

  /**
   * Computes an image warped by applying specified shifts.
   * @param u input array of shifts.
   * @param g input array for the image to be warped.
   * @param h output array for the warped image.
   */
  public void applyShifts(float[][] u, float[][] g, float[][] h) {
    final int n1 = u[0].length;
    final int n2 = u.length;
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
      for (int i1=0; i1<n1; ++i1) {
        hf[i2][i1] = si.interpolate(i1+uf[i2][i1]);
      }
    }});
  }

  /**
   * Computes an image warped by applying specified shifts.
   * @param u input array of shifts.
   * @param g input array for the image to be warped.
   * @param h output array for the warped image.
   */
  public void applyShifts(float[][][] u, float[][][] g, float[][][] h) {
    int n3 = u.length;
    final float[][][] uf = u;
    final float[][][] gf = g;
    final float[][][] hf = h;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      applyShifts(uf[i3],gf[i3],hf[i3]);
    }});
  }

  ///////////////////////////////////////////////////////////////////////////
  // for research and atypical applications

  /**
   * Returns normalized alignment errors for all samples and lags.
   * The number of lags nl = 1+shiftMax-shiftMin. Lag indices 
   * il = 0, 1, 2, ..., nl-1 correspond to integer shifts in 
   * [shiftMin,shiftMax]. Alignment errors are a monotonically
   * increasing function of |f[i1]-g[i1+il+shiftMin]|.
   * @param f array[n1] for the sequence f[i1].
   * @param g array[n1] for the sequence g[i1].
   * @return array[n1][nl] of alignment errors.
   */
  public float[][] computeErrors(float[] f, float[] g) {
    int n1 = f.length;
    float[][] e = new float[n1][_nl1];
    computeErrors(f,g,e);
    normalizeErrors(e);
    return e;
  }
  
  /**
   * Returns normalized alignment errors for all samples and fractional 
   * lags. The number of lags nl = (shiftMax-shiftMin)*_fr+1. Lag indices
   * il = 0, 1, 2, ..., nl-1 correspond to fractional shifts in 
   * [shiftMin,shiftMax], at intervals dx = 1.0/fr. Alignment errors are
   * a monotonically increasing function of |f[i1]-g[i1+il*dx+shiftMin]|
   * where values at g[i1+il*dx+shiftMin] are interpolated using a
   * SincInterpolator.
   * @param f array[n1] for the sequence f[i1].
   * @param g array[n1] for the sequence g[i1].
   * @return array[n1][nl] of alignment errors.
   */
  public float[][] computeErrorsFrac(float[] f, float[] g) {
    int n1 = f.length;
    float[][] e = new float[n1][_nl1];
    computeErrorsFrac(f,g,e);
    normalizeErrors(e);
    return e;
  }
  
  public float[][][] computeErrors(float[] pp, float[] ps1, float[] ps2) {
    int n1 = pp.length;
    float[][][] e = new float[n1][_nl1][_nlS];
    computeErrors(pp, ps1, ps2, e);
    normalizeErrors(e);
    return e;
  }

  /**
   * Returns normalized alignment errors for all samples and lags.
   * The number of lags nl = 1+shiftMax-shiftMin. Lag indices 
   * il = 0, 1, 2, ..., nl-1 correspond to integer shifts in 
   * [shiftMin,shiftMax]. Alignment errors are a monotonically
   * increasing function of |f[i2][i1]-g[i2][i1+il+shiftMin]|.
   * @param f array[n2][n1] for the image f[i2][i1].
   * @param g array[n2][n1] for the image g[i2][i1].
   * @return array[n2][n1][nl] of alignment errors.
   */
  public float[][][] computeErrors(float[][] f, float[][] g) {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] ff = f;
    final float[][] gf = g;
    final float[][][] ef = new float[n2][n1][_nl1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrors(ff[i2],gf[i2],ef[i2]);
    }});
    normalizeErrors(ef);
    return ef;
  }
  
  /**
   * Returns normalized alignment errors for all samples and fractional 
   * lags. The number of lags nl = (shiftMax-shiftMin)*_fr+1. Lag indices
   * il = 0, 1, 2, ..., nl-1 correspond to fractional shifts in 
   * [shiftMin,shiftMax], at intervals dx = 1.0/fr. Alignment errors are
   * a monotonically increasing function of 
   * |f[i2][i1]-g[i2][i1+il*dx+shiftMin]| where values at 
   * g[i2][i1+il*dx+shiftMin] are interpolated using a SincInterpolator.
   * @param f array[n2][n1] for the sequence f[i2][i1].
   * @param g array[n2][n1] for the sequence g[i2][i1].
   * @return array[n2][n1][nl] of alignment errors.
   */
  public float[][][] computeErrorsFrac(float[][] f, float[][] g) {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] ff = f;
    final float[][] gf = g;
    final float[][][] ef = new float[n2][n1][_nl1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      computeErrorsFrac(ff[i2],gf[i2],ef[i2]);
    }});
    normalizeErrors(ef);
    return ef;
  }

  /**
   * Returns normalized 1D alignment errors for 2D images.
   * The number of lags nl = 1+shiftMax-shiftMin. Lag indices 
   * il = 0, 1, 2, ..., nl-1 correspond to integer shifts in 
   * [shiftMin,shiftMax].
   * @param f array[n2][n1] for the image f[i2][i1].
   * @param g array[n2][n1] for the image g[i2][i1].
   * @return array[n1][nl] of alignment errors.
   */
  public float[][] computeErrors1(float[][] f, float[][] g) {
    final float[][] ff = f;
    final float[][] gf = g;
    final int nl = 1+_l1max-_l1min;
    final int n1 = f[0].length;
    final int n2 = f.length;
    float[][] e = Parallel.reduce(n2,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i2) {
      float[][] e = new float[n1][nl];
      computeErrors(ff[i2],gf[i2],e);
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    normalizeErrors(e);
    return e;
  }
  
  public float[][][] computeErrors1(
      float[][] pp, float[][] ps1, float[][] ps2)
  {
    final int n1 = pp[0].length;
    final int n2 = pp.length;
    final float[][] fpp = pp;
    final float[][] fps1 = ps1;
    final float[][] fps2 = ps2;
    float[][][] e = Parallel.reduce(n2, new Parallel.ReduceInt<float[][][]>() {
    public float[][][] compute(int i2) {
      float[][][] e = new float[n1][_nl1][_nlS];
      computeErrors(fpp[i2],fps1[i2],fps2[i2],e);
      return e;
    }
    public float[][][] combine(float[][][] ea, float[][][] eb) {
      return add(ea,eb);
    }});
    normalizeErrors(e);
    return e;
  }

  /**
   * Returns normalized 1D alignment errors for 3D images.
   * The number of lags nl = 1+shiftMax-shiftMin. Lag indices 
   * il = 0, 1, 2, ..., nl-1 correspond to integer shifts in 
   * [shiftMin,shiftMax].
   * @param f array[n3][n2][n1] for the image f[i3][i2][i1].
   * @param g array[n3][n2][n1] for the image g[i3][i2][i1].
   * @return array[n1][nl] of alignment errors.
   */
  public float[][] computeErrors1(float[][][] f, float[][][] g) {
    final float[][][] ff = f;
    final float[][][] gf = g;
    final int nl = 1+_l1max-_l1min;
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    float[][] e = Parallel.reduce(n2*n3,new Parallel.ReduceInt<float[][]>() {
    public float[][] compute(int i23) {
      int i2 = i23%n2;
      int i3 = i23/n2;
      float[][] e = new float[n1][nl];
      computeErrors(ff[i3][i2],gf[i3][i2],e);
      return e;
    }
    public float[][] combine(float[][] ea, float[][] eb) {
      return add(ea,eb);
    }});
    normalizeErrors(e);
    return e;
  }

  /**
   * Returns smoothed (and normalized) alignment errors.
   * @param e array[n1][nl] of alignment errors.
   * @return array[n1][nl] of smoothed errors.
   */
  public float[][] smoothErrors(float[][] e) {
    float[][] es = like(e);
    smoothErrors(e,es);
    return es;
  }

  /**
   * Returns smoothed (and normalized) alignment errors.
   * @param e array[n2][n1][nl] of alignment errors.
   * @return array[n2][n1][nl] of smoothed errors.
   */
  public float[][][] smoothErrors(float[][][] e) {
    float[][][] es = like(e);
    smoothErrors(e,es);
    return es;
  }

  /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n1][nl] of alignment errors.
   * @param es output array[n1][nl] of smoothed errors.
   */
  public void smoothErrors(float[][] e, float[][] es) {
    smoothErrors1(_bstrain1,e,es);
    normalizeErrors(es);
  }

  /**
   * Smooths (and normalizes) alignment errors.
   * Input and output arrays can be the same array.
   * @param e input array[n2][n1][nl] of alignment errors.
   * @param es output array[n2][n1][nl] of smoothed errors.
   */
  public void smoothErrors(float[][][] e, float[][][] es) {
    smoothErrors1(_bstrain1,e,es);
    normalizeErrors(es);
    smoothErrors2(_bstrain2,es,es);
    normalizeErrors(es);
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

  /**
   * Returns errors accumulated in forward direction.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][] accumulateForward(float[][] e) {
    float[][] d = like(e);
    accumulateForward(e,d);
    return d;
  }
  
  public float[][][] accumulateForward(float[][][] e) {
    float[][][] d = like(e);
    accumulateForward(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][] accumulateReverse(float[][] e) {
    float[][] d = like(e);
    accumulateReverse(e,d);
    return d;
  }
  
  public float[][][] accumulateReverse(float[][][] e) {
    float[][][] d = like(e);
    accumulateReverse(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in forward direction in 1st dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateForward1(float[][][] e) {
    float[][][] d = like(e);
    accumulateForward1(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction in 1st dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateReverse1(float[][][] e) {
    float[][][] d = like(e);
    accumulateReverse1(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in forward direction in 2nd dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateForward2(float[][][] e) {
    float[][][] d = like(e);
    accumulateForward2(e,d);
    return d;
  }

  /**
   * Returns errors accumulated in reverse direction in 2nd dimension.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   */
  public float[][][] accumulateReverse2(float[][][] e) {
    float[][][] d = like(e);
    accumulateReverse2(e,d);
    return d;
  }

  /**
   * Accumulates alignment errors in forward direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward(float[][] e, float[][] d) {
    accumulate( 1,_bstrain1,e,d);
  }
  
  public void accumulateForward(float[][][] e, float[][][] d) {
//    accumulate( 1,_bstrain1,_bstrainS1,e,d);
    accumulate( 1,_bstrainS1,e,d);
  }

  /**
   * Accumulates alignment errors in reverse direction.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse(float[][] e, float[][] d) {
    accumulate(-1,_bstrain1,e,d);
  }
  
  public void accumulateReverse(float[][][] e, float[][][] d) {
//    accumulate(-1,_bstrain1,_bstrainS1,e,d);
    accumulate(-1,_bstrainS1,e,d);
  }

  /**
   * Accumulates alignment errors in forward direction in 1st dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward1(float[][][] e, float[][][] d) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulateForward(e[i2],d[i2]);
  }

  /**
   * Accumulates alignment errors in reverse direction in 1st dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse1(float[][][] e, float[][][] d) {
    int n2 = e.length;
    for (int i2=0; i2<n2; ++i2)
      accumulateReverse(e[i2],d[i2]);
  }

  /**
   * Accumulates alignment errors in forward direction in 2nd dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateForward2(float[][][] e, float[][][] d) {
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] di1 = new float[n2][];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        ei1[i2] = e[i2][i1];
        di1[i2] = d[i2][i1];
      }
      accumulate( 1,_bstrain2,ei1,di1);
    }
  }

  /**
   * Accumulates alignment errors in reverse direction in 2nd dimension.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   */
  public void accumulateReverse2(float[][][] e, float[][][] d) {
    int n1 = e[0].length;
    int n2 = e.length;
    float[][]  ei1 = new float[n2][];
    float[][] di1 = new float[n2][];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        ei1[i2] = e[i2][i1];
        di1[i2] = d[i2][i1];
      }
      accumulate(-1,_bstrain2,ei1,di1);
    }
  }

  public float[] backtrackForward(float[][] d, float[][] e) {
    float[] u = new float[d.length];
    backtrackForward(d,e,u);
    return u;
  }
  
  public float[][] backtrackForward(float[][][] d, float[][][] e) {
    float[] u1 = new float[d.length];
    float[] uS = new float[d.length];
    backtrackForward(d,e,u1,uS);
    return new float[][]{u1,uS};
  }
  
  /**
   * Returns shifts found by backtracking in reverse.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[] backtrackReverse(float[][] d, float[][] e) {
    float[] u = new float[d.length];
    backtrackReverse(d,e,u);
    return u;
  }
  
  public float[][] backtrackReverse(float[][][] d, float[][][] e) {
    float[] u1 = new float[d.length];
    float[] uS = new float[d.length];
    backtrackReverse(d,e,u1,uS);
    return new float[][]{u1,uS};
  }

  /**
   * Returns shifts found by backtracking in reverse in 1st dimension.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[][] backtrackReverse1(float[][][] d, float[][][] e) {
    float[][] u = new float[d.length][d[0].length];
    backtrackReverse1(d,e,u);
    return u;
  }

  /**
   * Returns shifts found by backtracking in reverse in 2nd dimension.
   * @param d array of accumulated errors.
   * @param e array of alignment errors.
   */
  public float[][] backtrackReverse2(float[][][] d, float[][][] e) {
    float[][] u = new float[d.length][d[0].length];
    backtrackReverse2(d,e,u);
    return u;
  }

  public void backtrackForward(float[][] d, float[][] e, float[] u) {
    backtrack(1,_bstrain1,_shifts1,d,e,u);
  }
  
  public void backtrackForward(
      float[][][] d, float[][][] e, float[] u1, float[] uS)
  {
//    backtrack(1,_bstrain1,_bstrainS1,_shifts1,_shiftsS,d,e,u1,uS);
    backtrack(1,_bstrainS1,_shifts1,_shiftsS,d,e,u1,uS);
  }
  
  /**
   * Computes shifts by backtracking in reverse direction.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse(float[][] d, float[][] e, float[] u) {
    backtrack(-1,_bstrain1,_shifts1,d,e,u);
  }
  
  public void backtrackReverse(
      float[][][] d, float[][][] e, float[] u1, float[] uS)
  {
//    backtrack(-1,_bstrain1,_bstrainS1,_shifts1,_shiftsS,d,e,u1,uS);
    backtrack(-1,_bstrainS1,_shifts1,_shiftsS,d,e,u1,uS);
  }

  /**
   * Computes shifts by backtracking in reverse direction in 1st dimension.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse1(float[][][] d, float[][][] e, float[][] u) {
    int n2 = d.length;
    for (int i2=0; i2<n2; ++i2)
      backtrackReverse(d[i2],e[i2],u[i2]);
  }

  /**
   * Computes shifts by backtracking in reverse direction in 2nd dimension.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   */
  public void backtrackReverse2(float[][][] d, float[][][] e, float[][] u) {
    int n1 = d[0].length;
    int n2 = d.length;
    float[][] di1 = new float[n2][];
    float[][] ei1 = new float[n2][];
    float[] ui1 = new float[n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        di1[i2] = d[i2][i1];
        ei1[i2] = e[i2][i1];
      }
      backtrack(-1,_bstrain2,_shifts1,di1,ei1,ui1);
      for (int i2=0; i2<n2; ++i2)
        u[i2][i1] = ui1[i2];
    }
  }

//  public float[][][] compositeShifts(float[] u1, float[] uS) {
//    int n1 = u1.length;
//    float[][][] u = zerofloat(_nlS,_nl1,n1);
//    for (int i1=0; i1<n1; ++i1) {
//      for (int il1=0; il1<_nl1; ++il1) {
//        for (int ilS=0; ilS<_nlS; ++ilS) {
//          u[]
//        }
//      }
//    }
//  }
  
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
    System.out.println("min="+mm.emin+", max="+mm.emax);
    shiftAndScale(mm.emin,mm.emax,e);
  }

  /**
   * Returns the sum of errors for specified shifts, rounded to integers.
   * @param e array[n1][nl] of errors.
   * @param u array[n1] of shifts.
   * @return the sum of errors.
   */
  public float sumErrors(float[][] e, float[] u) {
    int n1 = e.length;
    int nl = e[0].length;
    float ul = 0.5f-_l1min;
    double sum = 0.0;
    for (int i1=0; i1<n1; ++i1) {
      int il = (int)(u[i1]+ul);
      il = max(0,min(nl-1,il));
      sum += e[i1][il];
    }
    return (float)sum;
  }

  /**
   * Returns the sum of errors for specified shifts, rounded to integers.
   * @param e array[n2][n1][nl] of errors.
   * @param u array[n2][n1] of shifts.
   * @return the sum of errors.
   */
  public float sumErrors(float[][][] e, float[][] u) {
    int n2 = e.length;
    double sum = 0.0;
    for (int i2=0; i2<n2; ++i2)
      sum += sumErrors(e[i2],u[i2]);
    return (float)sum;
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

  private int _nl1; // number of pp to ps1 lags
  private int _nlS; // number of ps1 to ps2 lags 
  private int _fr; // fractional shift factor (1.0/_fr is the shift interval)
  private int _l1min,_l1max; // min,max pp to ps1 lags (corresponds to u1)
  private int _lSmin,_lSmax; // min,max ps1 to ps2 lags (corresponds to us)
  private Sampling _shifts1; // sampling of pp to ps1 shift values
  private Sampling _shiftsS; // sampling of ps1 to ps2 shift values
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private int _esmooth = 0; // number of nonlinear smoothings of errors
  private double _usmooth1 = 0.0; // extent of smoothing shifts in 1st dim
  private double _usmooth2 = 0.0; // extent of smoothing shifts in 2nd dim
  private double _usmooth3 = 0.0; // extent of smoothing shifts in 3rd dim
  private int _bstrain1 = 1; // inverse of bound on strain 1 in 1st dimension
  private int _bstrain2 = 1; // inverse of bound on strain 1 in 2nd dimension
  private int _bstrain3 = 1; // inverse of bound on strain 1 in 3rd dimension
  private int _bstrainS1 = 1; // inverse of bound on strain S in 1st dimension
  private int _bstrainS2 = 1; // inverse of bound on strain S in 2nd dimension
  private int _bstrainS3 = 1; // inverse of bound on strain S in 3rd dimension
  private RecursiveExponentialFilter _ref1; // for smoothing shifts
  private RecursiveExponentialFilter _ref2; // for smoothing shifts
  private RecursiveExponentialFilter _ref3; // for smoothing shifts
  private SincInterpolator _si; // for warping with non-integer shifts
  private int _owl2 = 50; // window size in 2nd dimension for 3D images
  private int _owl3 = 50; // window size in 3rd dimension for 3D images
  private double _owf2 = 0.5; // fraction of window overlap in 2nd dimension
  private double _owf3 = 0.5; // fraction of window overlap in 3rd dimension

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

  private void updateSmoothingFilters() {
    _ref1 = (_usmooth1<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth1*_bstrain1);
    _ref2 = (_usmooth2<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth2*_bstrain2);
    _ref3 = (_usmooth3<=0.0) ? null :
      new RecursiveExponentialFilter(_usmooth3*_bstrain3);
  }

  /**
   * Computes alignment errors, not normalized.
   * @param f input array[ni] for sequence f.
   * @param g input array[ni] for sequence g.
   * @param e output array[ni][nl] of alignment errors.
   */
  private void computeErrors(float[] f, float[] g, float[][] e) {
    int n1 = f.length;
    int nl = _nl1;
    int n1m = n1-1;

    // Notes for indexing:
    // 0 <= il < nl, where il is index for lag
    // 0 <= i1 < n1, where i1 is index for sequence f
    // 0 <= j1 < n1, where j1 is index for sequence g
    // j1 = i1+il+lmin, where il+lmin = lag
    // 0 <= i1+il+lmin < n1, so that j1 is in bounds
    // max(0,-lmin-i1) <= il < min(nl,n1-lmin-i1)
    // max(0,-lmin-il) <= i1 < min(n1,n1-lmin-il)
    // j1 = 0    => i1 =     -lmin-il
    // j1 = n1-1 => i1 = n1-1-lmin-il

    // Compute errors where indices are in bounds for both f and g.
    for (int i1=0; i1<n1; ++i1) {
      int ilhi = min(nl,n1-_l1min-i1); // see notes above
      int k1 = n1m-_l1min-ilhi;
      for (int il=0,j1=i1+_l1min; il<nl; ++il,++j1) {
        if (il>=ilhi)
          e[i1][il] = e[k1][il];
        else {
          float ei = error(f[i1],g[j1]);
          e[i1][il] = ei;
        }
      }
    }
  }
  
  private void computeErrors(
      float[] pp, float[] ps1, float[] ps2, float[][][] e)
  {
    int n1 = pp.length;
    float e1m=0.0f,e2m=0.0f,e3m=0.0f; // last computed values for extrapolation

    // Notes for indexing:
    // 0 <= il1 < _nl1, where il1 is index for lag between pp and ps1
    // 0 <= ilS < _nlS, where ilS is index for lag between ps1 and ps2    
    // 0 <= i1 < n1, where i1 is index for sequence pp
    // 0 <= j1 < n1, where j1 is index for sequence ps1
    // 0 <= jS < n1, where jS is index for sequence ps2
    // j1 = i1+_l1min, where _l1min is the minimum shift u1 where u2 = u1 + uS
    // jS = j1+_lSmin, where _lSmin is the minimum shift uS where u2 = u1 + uS

    // Compute errors for pp, ps1, and ps2.
    for (int i1=0; i1<n1; ++i1) {
      for (int il1=0,j1=i1+_l1min; il1<_nl1; ++il1,++j1) {
        float e1;
        if (j1<n1) {
          e1 = error(pp[i1],ps1[j1]);
          e1m = e1;
        } else
          e1 = e1m;
        for (int ilS=0,jS=j1+_lSmin; ilS<_nlS; ++ilS,++jS) {
          float e2,e3;
          if (jS<n1) {
            e2 = error(pp [i1],ps2[jS]);
            e3 = error(ps1[j1],ps2[jS]);
            e2m = e2;
            e3m = e3;
          } else {
            e2 = e2m;
            e3 = e3m;
          }
          e[i1][il1][ilS] = e1+e2+e3;
        }
      }
    }
  }
  
  /**
   * Computes alignment errors for fractional shifts, not normalized
   * @param f input array[ni] for sequence f.
   * @param g input array[ni] for sequence g.
   * @param e output array[ni][_nl] of alignment errors.
   */
  private void computeErrorsFrac(float[] f, float[] g, float[][] e) { 
    int n1 = f.length;
    int nx = (n1-1)*_fr+1;
    int n1m = n1-1;
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,1.0f,0.0f);
    si.setUniformSamples(g);
    float[] gi = new float[nx];
    si.interpolate(nx,_shifts1.getDelta(),0.0,gi);

    // Notes for indexing:
    // _fr = fractional factor. Error is measured at interval 1.0/_fr
    // 0 <= ilx < _nl, where ilx is the fractional lag index, and _nl has
    //                 length (_lmax-_lmin)*fr+1
    // 0 <=  i1 <  n1, where i1 is index for sequence f
    // 0 <=  j1 <  nx, where j1 is the index for interpolated sequence gi, and
    //                 nx is the length of gi
    // j1 = _fr*(i1+_lmin), the current index in gi
    // ilxhi = min(_nl,nx+_fr*(-_lmin-i1)), where ilxhi is the maximum ilx that
    //                                      is in bounds of sequence g, 
    //                                      relative to the current index of
    //                                      sequence f

    // Compute errors where indices are in bounds for both f and g.
    for (int i1=0; i1<n1; ++i1) {
      int ilxhi = min(_nl1,nx+_fr*(-_l1min-i1)); // see notes above
      int k1 = n1m-_l1min-(ilxhi/_fr);
      for (int ilx=0,j1=_fr*(i1+_l1min); ilx<_nl1; ++ilx,++j1) {
        if (ilx>=ilxhi)
          e[i1][ilx] = e[k1][ilx];
        else {
          float ei = error(f[i1], gi[j1]);
          e[i1][ilx] = ei;
        }  
      }
    }
  }

  /**
   * Non-linear accumulation of alignment errors.
   * @param dir accumulation direction, positive or negative.
   * @param b sample offset used to constrain changes in lag.
   * @param e input array[ni][nl] of alignment errors.
   * @param d output array[ni][nl] of accumulated errors.
   */
  private static void accumulate(int dir, int b, float[][] e, float[][] d) {
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
      int jb = max(0,min(nim1,ii-is*b)); // ji scaled by b, for strain limits
      for (int il=0; il<nl; ++il) {
        int ilc = il+ic;
        ilc = (ilc==-1)?0:(ilc==nl)?nlm1:ilc; // index of lag constraint
        float dc = d[jb][ilc];
        float di = d[ji][il ];
        for (int kb=ji; kb!=jb; kb-=is) { // adjust for strain limits
          dc += e[kb][ilc];
        }
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
  
  private static void accumulate(
      int dir, int bS, float[][][] e, float[][][] d)
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
      int jbS = max(0,min(nim1,ii-is*bS));
      for (int il1=0; il1<nl1; ++il1) {
        int il1c = il1+ic;
        il1c = (il1c==-1)?0:(il1c==nl1)?nl1m1:il1c;
        for (int ilS=0; ilS<nlS; ++ilS) {
          int ilSc = ilS+ic;
          ilSc = (ilSc==-1)?0:(ilSc==nlS)?nlSm1:ilSc;
          float dc1 = d[jbS][il1c][ilS ];
          float dcS = d[jbS][il1 ][ilSc];
          float dc1S= d[jbS][il1c][ilSc];
          float di  = d[jbS][il1 ][ilS ];
          for (int kbS=ji; kbS!=jbS; kbS-=is) {
            dc1 += e[kbS][il1c][ilS ];
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
          if (dl==dc1)
            il1 = ilc1;
          else if (dl==dcS)
            ilS = ilcS;
          else if (dl==dc1S) {
            il1 = ilc1;
            ilS = ilcS;
          }
        }
        ii += is;
        u1[ii] = (float)shifts1.getValue(il1);
        uS[ii] = (float)shifts1.getValue(ilS);
        
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
          for (int kb1=ji; kb1!=jbS; kb1+=is) {
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
  
  private static void backtrack(
      int dir, int bS, Sampling shifts1, Sampling shiftsS,
      float[][][] d, float[][][] e, float[] u1, float[] uS) 
    {
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
        int jbS = max(0,min(nim1,ii+is*bS));
        int ilc1 = il1+ic;
        int ilcS = ilS+ic;
        ilc1 = (ilc1==-1)?0:(ilc1==nl1)?nl1m1:ilc1;
        ilcS = (ilcS==-1)?0:(ilcS==nlS)?nlSm1:ilcS;
        float dc1 = d[jbS][ilc1][ilS ];
        float dcS = d[jbS][il1 ][ilcS];
        float dc1S= d[jbS][ilc1][ilcS];
        float di  = d[jbS][il1 ][ilS ];
        for (int kbS=ji; kbS!=jbS; kbS+=is) {
          dc1 += e[kbS][ilc1][ilS ];
          dcS += e[kbS][il1 ][ilcS];
          dc1S+= e[kbS][ilc1][ilcS];
        }
        dl = min4(dc1,dcS,dc1S,di);
        if (dl!=di) {
          if (dl==dc1)
            il1 = ilc1;
          else if (dl==dcS)
            ilS = ilcS;
          else if (dl==dc1S) {
            il1 = ilc1;
            ilS = ilcS;
          }
        }
        ii += is;
        u1[ii] = (float)shifts1.getValue(il1);
        uS[ii] = (float)shifts1.getValue(ilS);
        
        // Adjust shifts for strain limits, if applicable.
        float du1 = (u1[ii]-u1[ii-is])*obS; // u1 must satisfy the strain
        float duS = (uS[ii]-uS[ii-is])*obS; // limits imposed by uS.
        if (il1==ilc1 && ilS==ilcS) {
          u1[ii] = u1[ii-is]+du1;
          uS[ii] = uS[ii-is]+duS;
          for (int kbS=ji; kbS!=jbS; kbS+=is) {
            ii += is;
            u1[ii] = u1[ii-is]+du1;
            uS[ii] = uS[ii-is]+duS;
          }
        } else if (il1==ilc1) {
          u1[ii] = u1[ii-is]+du1;
          for (int kb1=ji; kb1!=jbS; kb1+=is) {
            ii += is;
            u1[ii] = u1[ii-is]+du1;
            uS[ii] = uS[ii-is];
          }
        } else if (ilS==ilcS) {
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
   * Smooths alignment errors in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors1(int b, float[][] e, float[][] es) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] ef = new float[n1][nl];
    float[][] er = new float[n1][nl];
    accumulate( 1,b,e,ef);
    accumulate(-1,b,e,er);
    for (int i1=0; i1<n1; ++i1)
      for (int il=0; il<nl; ++il)
        es[i1][il] = ef[i1][il]+er[i1][il]-e[i1][il];
  }

  /**
   * Smooths alignment errors in 1st dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 1st dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors1(int b, float[][][] e, float[][][] es) {
    final int n2 = e.length;
    final int bf = b;
    final float[][][] ef = e;
    final float[][][] esf = es;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      smoothErrors1(bf,ef[i2],esf[i2]);
    }});
  }

  /**
   * Smooths alignment errors in 2nd dimension.
   * Does not normalize errors after smoothing.
   * @param b strain parameter in 2nd dimension.
   * @param e input array of alignment errors to be smooothed.
   * @param es output array of smoothed alignment errors.
   */
  private static void smoothErrors2(int b, float[][][] e, float[][][] es) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final int bf = b;
    final float[][][]  ef = e;
    final float[][][] esf = es;
    final Parallel.Unsafe<float[][][]> eeu = 
      new Parallel.Unsafe<float[][][]>();
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][][] ee = eeu.get();
      if (ee==null) eeu.set(ee=new float[4][n2][nl]);
      float[][]  e1 = ee[0];
      float[][] es1 = ee[1];
      float[][] ef1 = ee[2];
      float[][] er1 = ee[3];
      for (int i2=0; i2<n2; ++i2) {
         e1[i2] =  ef[i2][i1];
        es1[i2] = esf[i2][i1];
        for (int il=0; il<nl; ++il) {
          ef1[i2][il] = 0.0f;
          er1[i2][il] = 0.0f;
        }
      }
      accumulate( 1,bf,e1,ef1);
      accumulate(-1,bf,e1,er1);
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          es1[i2][il] = ef1[i2][il]+er1[i2][il]-e1[i2][il];
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

  ///////////////////////////////////////////////////////////////////////////
  // for 3D image warping

  private void computeErrors(float[][][] f, float[][][] g, float[][][][] e) {
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][] ff = f;
    final float[][][] gf = g;
    final float[][][][] ef = e;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        computeErrors(ff[i3][i2],gf[i3][i2],ef[i3][i2]);
      }
    }});
    normalizeErrors(e);
  }
  private static void normalizeErrors(float[][][][] e) {
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
  private void smoothErrors(float[][][][] e) {
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      smoothErrors1(_bstrain1,ef[i3],ef[i3]);
    }});
    normalizeErrors(e);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      smoothErrors2(_bstrain2,ef[i3],ef[i3]);
    }});
    normalizeErrors(e);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] ei2 = new float[n3][][];
      for (int i3=0; i3<n3; ++i3)
        ei2[i3] = ef[i3][i2];
      smoothErrors2(_bstrain3,ei2,ei2);
    }});
    normalizeErrors(e);
  }
  private void computeShifts(float[][][][] e, float[][][] u) {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    final float[][][] uf = u;
    final Parallel.Unsafe<float[][]> du = new Parallel.Unsafe<float[][]>();
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      float[][] d = du.get();
      if (d==null) du.set(d=new float[n1][nl]);
      for (int i2=0; i2<n2; ++i2) {
        accumulateForward(ef[i3][i2],d);
        backtrackReverse(d,ef[i3][i2],uf[i3][i2]);
      }
    }});
  }
  private void smoothShifts(float[][][] u) {
    if (_ref1!=null) _ref1.apply1(u,u);
    if (_ref2!=null) _ref2.apply2(u,u);
    if (_ref3!=null) _ref3.apply3(u,u);
  }
  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }
  private static class OverlappingWindows2 {
    public OverlappingWindows2(
      int n1, int n2, int l1, int l2, double f1, double f2) 
    {
      Check.argument(0.0<=f1 && f1<1.0,"0 <= f1 < 1");
      Check.argument(0.0<=f2 && f2<1.0,"0 <= f2 < 1");
      _n1 = n1;
      _n2 = n2;
      _l1 = min(l1,n1);
      _l2 = min(l2,n2);
      _m1 = 1+(int)ceil((_n1-_l1)/(_l1*(1.0-f1)));
      _m2 = 1+(int)ceil((_n2-_l2)/(_l2*(1.0-f2)));
      _s1 = (double)(_n1-_l1)/max(1,_m1-1);
      _s2 = (double)(_n2-_l2)/max(1,_m2-1);
      makeWeights();
      makeScalars();
    }
    public int getN1() { return _n1; }
    public int getN2() { return _n2; }
    public int getL1() { return _l1; }
    public int getL2() { return _l2; }
    public int getM1() { return _m1; }
    public int getM2() { return _m2; }
    public int getI1(int k1) { return (int)(k1*_s1+0.5); }
    public int getI2(int k2) { return (int)(k2*_s2+0.5); }
    public float getWeight(int i1, int i2, int j1, int j2) {
      return _w[j2][j1]*_s[i2+j2][i1+j1];
    }
    public float[][] getWeights() { return _w; }
    public float[][] getScalars() { return _s; }
    private int _n1,_n2; // numbers of samples
    private int _l1,_l2; // window lengths
    private int _m1,_m2; // numbers of windows
    private double _s1,_s2; // nominal window spacings
    private float[][] _w; // weights[l2][l1] for windowing
    private float[][] _s; // scalars[n2][n1] for normalization
    private void makeWeights() {
      _w = new float[_l2][_l1];
      for (int i2=0; i2<_l2; ++i2) {
        for (int i1=0; i1<_l1; ++i1) {
          double s1 = sin((i1+1.0)*PI/(_l1+1.0));
          double s2 = sin((i2+1.0)*PI/(_l2+1.0));
          _w[i2][i1] = (float)(s1*s1*s2*s2);
        }
      }
    }
    private void makeScalars() {
      _s = new float[_n2][_n1];
      for (int k2=0; k2<_m2; ++k2) {
        int i2 = getI2(k2);
        for (int k1=0; k1<_m1; ++k1) {
          int i1 = getI1(k1);
          for (int j2=0; j2<_l2; ++j2) {
            for (int j1=0; j1<_l1; ++j1) {
              _s[i2+j2][i1+j1] += _w[j2][j1];
            }
          }
        }
      }
      for (int i2=0; i2<_n2; ++i2) {
        for (int i1=0; i1<_n1; ++i1) {
          _s[i2][i1] = 1.0f/_s[i2][i1];
        }
      }
    }
  }
}
