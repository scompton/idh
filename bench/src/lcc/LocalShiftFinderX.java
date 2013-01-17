/****************************************************************************
Copyright (c) 2006, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
the Common Public License - v1.0, which accompanies this distribution, and is 
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package lcc;

import edu.mines.jtk.la.TridiagonalFMatrix;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Estimates displacement vector fields for two images. For example, given 
 * two 2-D images f(x1,x2) and g(x1,x2), a shift finder estimates two vector 
 * components of displacement u1(x1,x2) and u2(x1,x2) such that 
 * f(x1,x2) ~ g(x1+u1(x1,x2),x2+u2(x1,x2)).
 * <p>
 * Like the images f and g, the components of displacement are sampled
 * functions of coordinates x1 and x2. That is, displacements may vary 
 * from sample to sample. The components u1 and u2 represent displacements 
 * in the x1 and x2 coordinate directions, respectively.
 * <p>
 * This shift finder estimates each component of displacement using local
 * cross-correlations. For each image sample, the estimated shift is that 
 * which yields the maximum correlation coefficient. This coefficient is
 * found by quadratic interpolation of correlation functions sampled at
 * integer lags.
 * <p>
 * The peak (maximum) correlation coefficient may be used to measure 
 * quality of an estimated shift. However, because a correlation function 
 * may have more than one peak (local maxima), a better measure of quality 
 * may be the difference between the coefficients for the correlation peak 
 * and next highest peak. Both the peak coefficient and this difference may 
 * be computed with the shifts.
 * <p>
 * Methods are provided to find and compensate for each component of shift 
 * sequentially. As each component is found, that component can be removed 
 * from the image g before estimating another component. For example, again 
 * for 2-D images f(x1,x2) and g(x1,x2), we might first estimate u1(x1,x2). 
 * If we then compute an image h(x1,x2) = g(x1+u1(x1,x2),x2), we can use
 * f(x1,x2) and h(x1,x2) to estimate u2(x1,x2). By repeating this process
 * sequentially, we obtain estimates for both u1(x1,x2) and u2(x1,x2) such
 * that f(x1,x2) ~ g(x1+u1(x1,x2),x2+u2(x1,x2)).
 * <p>
 * Methods are also provided to whiten 2-D and 3-D images before estimating
 * displacement vectors. This (spectral) whitening improves estimates of
 * displacements parallel to image features that may be otherwise poorly
 * resolved. Whitening is performed with local prediction error filters
 * computed from local auto-correlations.
 *
 * @author Dave Hale, Colorado School of Mines
 * @version 2006.11.18
 */
public class LocalShiftFinderX {

  /**
   * Construct a shift estimator with specified parameters.
   * When applied to multi-dimensional arrays, the estimator has the
   * same correlation window half-width for all dimensions.
   * @param sigma the correlation window half-width; must not be less than 1.
   */
  public LocalShiftFinderX(double sigma) {
    this(sigma,sigma,sigma);
  }

  /**
   * Construct a shift estimator with specified parameters.
   * When applied to multi-dimensional arrays, the estimator has half-width 
   * sigma1 for the 1st dimension and half-width sigma2 for 2nd and higher 
   * dimensions.
   * @param sigma1 correlaton window half-width for 0st dimension; 
   *  must not be less than 1.
   * @param sigma2 correlation window half-width for 2nd and higher 
   *  dimensions; must not be less than 1.
   */
  public LocalShiftFinderX(double sigma1, double sigma2) {
    this(sigma1,sigma2,sigma2);
  }

  /**
   * Construct a shift estimator with specified parameters.
   * When applied to multi-dimensional arrays, the estimator has half-width 
   * sigma1 for the 1st dimension, half-width sigma2 for the 2nd dimension, 
   * and half-width sigma3 for 3rd and higher dimensions.
   * @param sigma1 correlation window half-width for 1st dimension; 
   *  must not be less than 1.
   * @param sigma2 correlation window half-width for 2nd dimension;
   *  must not be less than 1.
   * @param sigma3 correlation window half-width for 3rd and higher 
   *  dimensions; must not be less than 1.
   */
  public LocalShiftFinderX(double sigma1, double sigma2, double sigma3) {
    Check.argument(sigma1>=1.0,"sigma1>=1.0");
    Check.argument(sigma2>=1.0,"sigma2>=1.0");
    Check.argument(sigma3>=1.0,"sigma3>=1.0");
    _sigma1 = (float)sigma1;
    _sigma2 = (float)sigma2;
    _sigma3 = (float)sigma3;
    _lcfSimple = new LocalCorrelationFilter(
      LocalCorrelationFilter.Type.SIMPLE,
      LocalCorrelationFilter.Window.GAUSSIAN,
      sigma1,sigma2,sigma3);
    _si = new SincInterp();
    _si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
    _lsf = new LocalSmoothingFilter();
  }

  /**
   * Enables or disables smoothing of shifts. If enabled, shifts are 
   * smoothed to prevent discontinuities caused by cycle skipping.
   * <p>
   * The default is that smoothing is disabled, because smoothing is 
   * currently supported for only 1D shifts. Smoothing is currently
   * not supported for 2D and 3D arrays of shifts.
   * @param smoothShifts true, for smoothing; false, otherwise.
   */
  public void setSmoothShifts(boolean smoothShifts) {
    _smoothShifts = smoothShifts;
  }

  /**
   * Enables or disables interpolation of displacements when shifting.
   * The default is to interpolate displacements. This is the most
   * accurate method when sequentially applying non-constant shifts.
   * @param enable true, to enable interpolation; false, to disable.
   */
  public void setInterpolateDisplacements(boolean enable) {
    _interpolateDisplacements = enable;
  }

  /**
   * Finds shifts in the 1st (and only) dimension.
   * @param min1 the minimum shift.
   * @param max1 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find1(
    int min1, int max1, float[] f, float[] g, float[] u) 
  {
    if (_smoothShifts) {
      findShiftsSmooth(min1,max1,f,g,u);
    } else {
      findShifts(min1,max1,f,g,u);
    }
  }

  /**
   * Finds shifts in the 1st dimension.
   * @param min1 the minimum shift.
   * @param max1 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find1(
    int min1, int max1, float[][] f, float[][] g, float[][] u) 
  {
    findShifts(1,min1,max1,f,g,u);
  }

  /**
   * Finds shifts in the 2nd dimension.
   * @param min2 the minimum shift.
   * @param max2 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find2(
    int min2, int max2, float[][] f, float[][] g, float[][] u) 
  {
    findShifts(2,min2,max2,f,g,u);
  }

  /**
   * Finds shifts in the 1st dimension.
   * @param min1 the minimum shift.
   * @param max1 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find1(
    int min1, int max1, float[][][] f, float[][][] g, float[][][] u) 
  {
    findShifts(1,min1,max1,f,g,u);
  }

  /**
   * Finds shifts in the 2nd dimension.
   * @param min2 the minimum shift.
   * @param max2 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find2(
    int min2, int max2, float[][][] f, float[][][] g, float[][][] u) 
  {
    findShifts(2,min2,max2,f,g,u);
  }

  /**
   * Finds shifts in the 3rd dimension.
   * @param min3 the minimum shift.
   * @param max3 the maximum shift.
   * @param f the input array f.
   * @param g the input array g.
   * @param u output array of shifts.
   */
  public void find3(
    int min3, int max3, float[][][] f, float[][][] g, float[][][] u) 
  {
    findShifts(3,min3,max3,f,g,u);
  }

  /**
   * Applies specified shift in the 1st (and only) dimension.
   * @param du input array of changes to displacements in 1st dimension.
   * @param u1 input/output array of displacements in 1st dimension.
   * @param h input/output array of image samples.
   */
  public void shift1(float[] du, float[] u1, float[] h) {
    int n1 = h.length;
    float[] xu1 = new float[n1];
    float[] u1a = u1;
    float[] u1b = new float[n1];
    float[] ha = h;
    float[] hb = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      xu1[i1] = (float)(i1)+du[i1];
    _si.interpolate(n1,1.0,0.0,ha,n1,xu1,hb);
      copy(hb,h);
    if (_interpolateDisplacements) {
      _si.interpolate(n1,1.0,0.0,u1a,n1,xu1,u1b);
      copy(u1b,u1);
    }
  }

  /**
   * Applies specified shift in the 1st dimension.
   * @param du input array of changes to displacements in 1st dimension.
   * @param u1 input/output array of displacements in 1st dimension.
   * @param u2 input/output array of displacements in 2nd dimension.
   * @param h input/output array of image samples.
   */
  public void shift1(float[][] du, float[][] u1, float[][] u2, float[][] h) {
    int n1 = h[0].length;
    int n2 = h.length;
    float[] xu1 = new float[n1];
    float[] u1b = new float[n1];
    float[] u2b = new float[n1];
    float[] hb = new float[n1];
    for (int i2=0; i2<n2; ++i2) {
      float[] ha = h[i2];
      float[] u1a = u1[i2];
      float[] u2a = u2[i2];
      float[] du1 = du[i2];
      for (int i1=0; i1<n1; ++i1) {
        xu1[i1] = (float)(i1)+du1[i1];
      }
      _si.interpolate(n1,1.0,0.0,ha,n1,xu1,hb);
      if (_interpolateDisplacements) {
        _si.interpolate(n1,1.0,0.0,u1a,n1,xu1,u1b);
        _si.interpolate(n1,1.0,0.0,u2a,n1,xu1,u2b);
      } else {
        copy(u1a,u1b);
        copy(u2a,u2b);
      }
      for (int i1=0; i1<n1; ++i1) {
        h[i2][i1] = hb[i1];
        u1[i2][i1] = u1b[i1]+du1[i1];
        u2[i2][i1] = u2b[i1];
      }
    }
  }

  /**
   * Applies specified shift in the 2nd dimension.
   * @param du input array of changes to displacements in 2nd dimension.
   * @param u1 input/output array of displacements in 1st dimension.
   * @param u2 input/output array of displacements in 2nd dimension.
   * @param h input/output array of image samples.
   */
  public void shift2(float[][] du, float[][] u1, float[][] u2, float[][] h) {
    int n1 = h[0].length;
    int n2 = h.length;
    float[] du2 = new float[n2];
    float[] xu2 = new float[n2];
    float[] u1a = new float[n2];
    float[] u1b = new float[n2];
    float[] u2a = new float[n2];
    float[] u2b = new float[n2];
    float[] ha = new float[n2];
    float[] hb = new float[n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        ha[i2] = h[i2][i1];
        u1a[i2] = u1[i2][i1];
        u2a[i2] = u2[i2][i1];
        du2[i2] = du[i2][i1];
        xu2[i2] = (float)(i2)+du2[i2];
      }
      _si.interpolate(n2,1.0,0.0,ha,n2,xu2,hb);
      if (_interpolateDisplacements) {
        _si.interpolate(n2,1.0,0.0,u1a,n2,xu2,u1b);
        _si.interpolate(n2,1.0,0.0,u2a,n2,xu2,u2b);
      } else {
        copy(u1a,u1b);
        copy(u2a,u2b);
      }
      for (int i2=0; i2<n2; ++i2) {
        h[i2][i1] = hb[i2];
        u1[i2][i1] = u1b[i2];
        u2[i2][i1] = u2b[i2]+du2[i2];
      }
    }
  }

  /**
   * Applies specified shift in the 1st dimension.
   * @param du input array of changes to displacements in 1st dimension.
   * @param u1 input/output array of displacements in 1st dimension.
   * @param u2 input/output array of displacements in 2nd dimension.
   * @param u3 input/output array of displacements in 3rd dimension.
   * @param h input/output array of image samples.
   */
  public void shift1(
    float[][][] du, float[][][] u1, float[][][] u2, float[][][] u3,
    float[][][] h) 
  {
    int n1 = h[0][0].length;
    int n2 = h[0].length;
    int n3 = h.length;
    float[] xu1 = new float[n1];
    float[] u1b = new float[n1];
    float[] u2b = new float[n1];
    float[] u3b = new float[n1];
    float[] hb = new float[n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        float[] ha = h[i3][i2];
        float[] u1a = u1[i3][i2];
        float[] u2a = u2[i3][i2];
        float[] u3a = u3[i3][i2];
        float[] du1 = du[i3][i2];
        for (int i1=0; i1<n1; ++i1) {
          xu1[i1] = (float)(i1)+du1[i1];
        }
        _si.interpolate(n1,1.0,0.0,ha,n1,xu1,hb);
        if (_interpolateDisplacements) {
          _si.interpolate(n1,1.0,0.0,u1a,n1,xu1,u1b);
          _si.interpolate(n1,1.0,0.0,u2a,n1,xu1,u2b);
          _si.interpolate(n1,1.0,0.0,u3a,n1,xu1,u3b);
        } else {
          copy(u1a,u1b);
          copy(u2a,u2b);
          copy(u3a,u3b);
        }
        for (int i1=0; i1<n1; ++i1) {
          h[i3][i2][i1] = hb[i1];
          u1[i3][i2][i1] = u1b[i1]+du1[i1];
          u2[i3][i2][i1] = u2b[i1];
          u3[i3][i2][i1] = u3b[i1];
        }
      }
    }
  }

  /**
   * Applies specified shift in the 2nd dimension.
   * @param du input array of changes to displacements in 2nd dimension.
   * @param u1 input/output array of displacements in 1st dimension.
   * @param u2 input/output array of displacements in 2nd dimension.
   * @param u3 input/output array of displacements in 3rd dimension.
   * @param h input/output array of image samples.
   */
  public void shift2(
    float[][][] du, float[][][] u1, float[][][] u2, float[][][] u3,
    float[][][] h) 
  {
    int n1 = h[0][0].length;
    int n2 = h[0].length;
    int n3 = h.length;
    float[] du2 = new float[n2];
    float[] xu2 = new float[n2];
    float[] u1a = new float[n2];
    float[] u1b = new float[n2];
    float[] u2a = new float[n2];
    float[] u2b = new float[n2];
    float[] u3a = new float[n2];
    float[] u3b = new float[n2];
    float[] ha = new float[n2];
    float[] hb = new float[n2];
    for (int i3=0; i3<n3; ++i3) {
      for (int i1=0; i1<n1; ++i1) {
        for (int i2=0; i2<n2; ++i2) {
          ha[i2] = h[i3][i2][i1];
          u1a[i2] = u1[i3][i2][i1];
          u2a[i2] = u2[i3][i2][i1];
          u3a[i2] = u3[i3][i2][i1];
          du2[i2] = du[i3][i2][i1];
          xu2[i2] = (float)(i2)+du2[i2];
        }
        _si.interpolate(n2,1.0,0.0,ha,n2,xu2,hb);
        if (_interpolateDisplacements) {
          _si.interpolate(n2,1.0,0.0,u1a,n2,xu2,u1b);
          _si.interpolate(n2,1.0,0.0,u2a,n2,xu2,u2b);
          _si.interpolate(n2,1.0,0.0,u3a,n2,xu2,u3b);
        } else {
          copy(u1a,u1b);
          copy(u2a,u2b);
          copy(u3a,u3b);
        }
        for (int i2=0; i2<n2; ++i2) {
          h[i3][i2][i1] = hb[i2];
          u1[i3][i2][i1] = u1b[i2];
          u2[i3][i2][i1] = u2b[i2]+du2[i2];
          u3[i3][i2][i1] = u3b[i2];
        }
      }
    }
  }

  /**
   * Applies specified shift in the 3rd dimension.
   * @param du input array of changes to displacements in 3rd dimension.
   * @param u1 input/output array of displacements in 1st dimension.
   * @param u2 input/output array of displacements in 2nd dimension.
   * @param u3 input/output array of displacements in 3rd dimension.
   * @param h input/output array of image samples.
   */
  public void shift3(
    float[][][] du, float[][][] u1, float[][][] u2, float[][][] u3,
    float[][][] h) 
  {
    int n1 = h[0][0].length;
    int n2 = h[0].length;
    int n3 = h.length;
    float[] du3 = new float[n3];
    float[] xu3 = new float[n3];
    float[] u1a = new float[n3];
    float[] u1b = new float[n3];
    float[] u2a = new float[n3];
    float[] u2b = new float[n3];
    float[] u3a = new float[n3];
    float[] u3b = new float[n3];
    float[] ha = new float[n3];
    float[] hb = new float[n3];
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int i3=0; i3<n3; ++i3) {
          ha[i3] = h[i3][i2][i1];
          u1a[i3] = u1[i3][i2][i1];
          u2a[i3] = u2[i3][i2][i1];
          u3a[i3] = u3[i3][i2][i1];
          du3[i3] = du[i3][i2][i1];
          xu3[i3] = (float)(i3)+du3[i3];
        }
        _si.interpolate(n3,1.0,0.0,ha,n3,xu3,hb);
        if (_interpolateDisplacements) {
          _si.interpolate(n3,1.0,0.0,u1a,n3,xu3,u1b);
          _si.interpolate(n3,1.0,0.0,u2a,n3,xu3,u2b);
          _si.interpolate(n3,1.0,0.0,u3a,n3,xu3,u3b);
        } else {
          copy(u1a,u1b);
          copy(u2a,u2b);
          copy(u3a,u3b);
        }
        for (int i3=0; i3<n3; ++i3) {
          h[i3][i2][i1] = hb[i3];
          u1[i3][i2][i1] = u1b[i3];
          u2[i3][i2][i1] = u2b[i3];
          u3[i3][i2][i1] = u3b[i3]+du3[i3];
        }
      }
    }
  }

  /**
   * Applies local prediction-error (spectal whitening) filters.
   * The input and output arrays f and g can be the same array.
   * @param f the input array.
   * @param g the output array.
   */
  public void whiten(float[][] f, float[][] g) {
    whiten(1.0,f,g);
  }

  /**
   * Applies local prediction-error (spectal whitening) filters.
   * The input and output arrays f and g can be the same array.
   * @param sigma half-width of Gaussian smoothing applied after whitening;
   *  less than one for no smoothing.
   * @param f the input array.
   * @param g the output array.
   */
  public void whiten(double sigma, float[][] f, float[][] g) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[][] r00 = new float[n2][n1];
    float[][] rpm = new float[n2][n1];
    float[][] rm0 = new float[n2][n1];
    float[][] r0m = new float[n2][n1];
    _lcfSimple.setInputs(f,f);
    _lcfSimple.correlate( 0, 0,r00);
    _lcfSimple.correlate( 1,-1,rpm);
    _lcfSimple.correlate(-1, 0,rm0);
    _lcfSimple.correlate( 0,-1,r0m);
    float[][] s = rm0;
    float[][] t = r0m;
    for (int i2=0; i2<n2; ++i2)
      s[i2][0] = 0.0f;
    for (int i1=0; i1<n1; ++i1)
      s[0][i1] = 0.0f;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        double b1 = rm0[i2][i1];
        double b2 = r0m[i2][i1];
        double a11 = r00[i2][i1-1];
        double a21 = rpm[i2][i1-1];
        double a22 = r00[i2-1][i1];
        double l11 = sqrt(a11);
        double l21 = a21/l11;
        double d22 = a22-l21*l21;
        double x1 = 0.0;
        double x2 = 0.0;
        if (d22>0.0) {
          double l22 = sqrt(d22);
          double v1 = b1/l11;
          double v2 = (b2-l21*v1)/l22;
          x2 = v2/l22;
          x1 = (v1-l21*x2)/l11;
        }
        float a1 = (float)x1;
        float a2 = (float)x2;
        s[i2][i1] = f[i2][i1]
                    - a1*f[i2][i1-1]
                    - a2*f[i2-1][i1];
      }
    }
    if (sigma>=1.0) {
      RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
      rgf.apply0X(s,t);
      rgf.applyX0(t,g);
    } else {
      copy(s,g);
    }
  }

  /**
   * Applies local prediction-error (spectal whitening) filters.
   * The input and output arrays f and g can be the same array.
   * Smooths the output with a Gaussian filter with half-width sigma = 1.0.
   * @param f the input array.
   * @param g the output array.
   */
  public void whiten(float[][][] f, float[][][] g) {
    whiten(1.0,f,g);
  }

  /**
   * Applies local prediction-error (spectal whitening) filters.
   * The input and output arrays f and g can be the same array.
   * @param sigma half-width of Gaussian smoothing applied after whitening;
   *  less than one for no smoothing.
   * @param f the input array.
   * @param g the output array.
   */
  public void whiten(double sigma, float[][][] f, float[][][] g) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[][][] r000 = new float[n3][n2][n1];
    float[][][] rpm0 = new float[n3][n2][n1];
    float[][][] rp0m = new float[n3][n2][n1];
    float[][][] r0pm = new float[n3][n2][n1];
    float[][][] rm00 = new float[n3][n2][n1];
    float[][][] r0m0 = new float[n3][n2][n1];
    float[][][] r00m = new float[n3][n2][n1];
    float[][][] s = rm00;
    float[][][] t = r0m0;
    _lcfSimple.setInputs(f,f);
    _lcfSimple.correlate( 0, 0, 0,r000);
    _lcfSimple.correlate( 1,-1, 0,rpm0);
    _lcfSimple.correlate( 1, 0,-1,rp0m);
    _lcfSimple.correlate( 0, 1,-1,r0pm);
    _lcfSimple.correlate(-1, 0, 0,rm00);
    _lcfSimple.correlate( 0,-1, 0,r0m0);
    _lcfSimple.correlate( 0, 0,-1,r00m);
    for (int i3=0; i3<n3; ++i3)
      for (int i2=0; i2<n2; ++i2)
        s[i3][i2][0] = 0.0f;
    for (int i3=0; i3<n3; ++i3)
      for (int i1=0; i1<n1; ++i1)
        s[i3][0][i1] = 0.0f;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1)
        s[0][i2][i1] = 0.0f;
    for (int i3=1; i3<n3; ++i3) {
      for (int i2=1; i2<n2; ++i2) {
        for (int i1=1; i1<n1; ++i1) {
          double b1 = rm00[i3][i2][i1];
          double b2 = r0m0[i3][i2][i1];
          double b3 = r00m[i3][i2][i1];
          double a11 = r000[i3][i2][i1-1];
          double a21 = rpm0[i3][i2][i1-1];
          double a22 = r000[i3][i2-1][i1];
          double a31 = rp0m[i3][i2][i1-1];
          double a32 = r0pm[i3][i2-1][i1];
          double a33 = r000[i3-1][i2][i1];
          double x1 = 0.0;
          double x2 = 0.0;
          double x3 = 0.0;
          double l11 = sqrt(a11);
          double l21 = a21/l11;
          double l31 = a31/l11;
          double d22 = a22-l21*l21;
          if (d22>0.0) {
            double l22 = sqrt(d22);
            double l32 = (a32-l31*l21)/l22;
            double d33 = a33-l31*l31-l32*l32;
            if (d33>0.0) {
              double l33 = sqrt(d33);
              double v1 = b1/l11;
              double v2 = (b2-l21*v1)/l22;
              double v3 = (b3-l31*v1-l32*v2)/l33;
              x3 = v3/l33;
              x2 = (v2-l32*x3)/l22;
              x1 = (v1-l21*x2-l31*x3)/l11;
            }
          }
          float a1 = (float)x1;
          float a2 = (float)x2;
          float a3 = (float)x3;
          s[i3][i2][i1] = f[i3][i2][i1]
                          - a1*f[i3][i2][i1-1]
                          - a2*f[i3][i2-1][i1]
                          - a3*f[i3-1][i2][i1];
        }
      }
    }
    if (sigma>=1.0) {
      RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(sigma);
      rgf.apply0XX(s,t);
      rgf.applyX0X(t,s);
      rgf.applyXX0(s,g);
    } else {
      copy(s,g);
    }
  }


  ///////////////////////////////////////////////////////////////////////////
  // private

  private float _sigma1,_sigma2,_sigma3;
  private LocalCorrelationFilter _lcfSimple;
  private SincInterp _si;
  private boolean _smoothShifts = false;
  private boolean _interpolateDisplacements = true;
  private LocalSmoothingFilter _lsf;

  private void findShiftsSmooth(
    int min, int max, float[] f, float[] g, float[] u) 
  {
    int n1 = f.length;

    // Parameters for highest and 2nd highest correlation peaks.
    float[] u1 = new float[n1];
    float[] c1 = new float[n1];
    float[] c2 = new float[n1];

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][] cabc = new float[3][n1];

    // Correlate for min lag.
    LocalCorrelationFilter lcf = _lcfSimple;
    lcf.setInputs(f,g);
    int lag1 = min;
    lcf.correlate(lag1,cabc[1]);
    lcf.normalize(lag1,cabc[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[] ca = (lag>min)?cabc[(i  )%3]:cabc[(i+2)%3];
      float[] cb =           cabc[(i+1)%3];
      float[] cc = (lag<max)?cabc[(i+2)%3]:cabc[(i  )%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = lag+1;
        lcf.correlate(lag1,cc);
        lcf.normalize(lag1,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, interpolate three correlation values with a quadratic
      // and then update parameters for the highest and 2nd highest peaks.
      for (int i1=0; i1<n1; ++i1) {
        float ai = ca[i1];
        float bi = cb[i1];
        float ci = cc[i1];
        if (bi>=ai && bi>=ci) {
          double q0 = bi;
          double q1 = 0.5*(ci-ai);
          double q2 = 0.5*(ci+ai)-bi;
          double up = (q2<0.0)?-0.5*q1/q2:0.0;
          double qp = q0+up*(q1+up*q2);
          if (qp>0.0) {
            float cp = (float)qp;
            if (cp>c1[i1]) {
              c2[i1] = c1[i1];
              c1[i1] = cp;
              u1[i1] = (float)(lag+up);
            } else if (cp>c2[i1]) {
              c2[i1] = cp;
            }
          }
        }
      }
    }

    // Use parameters for correlation peaks to apply smoothing filter.
    float c = 0.01f*_sigma1*_sigma1;
    float tiny = 0.01f;
    float cpow = 2.0f;
    float[] s = cabc[0];
    float cc = pow(c1[0]-c2[0],cpow);
    s[0] = (1.0f-cc)/(tiny+cc);
    for (int i1=1; i1<n1-1; ++i1) {
      cc = pow(c1[i1]-c2[i1],cpow);
      s[i1] = (1.0f-cc)/(tiny+cc);
    }
    cc = pow(c1[n1-1]-c2[n1-1],cpow);
    s[n1-1] = (1.0f-cc)/(tiny+cc);
    _lsf.apply(c,s,u1,u);
  }

  private void findShiftsSmoothVersion2(
    int min, int max, float[] f, float[] g, float[] u) 
  {
    int n1 = f.length;

    // Parameters for highest and 2nd highest correlation peaks.
    float[] u1 = new float[n1];
    float[] c1 = new float[n1];
    float[] c2 = new float[n1];

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][] cabc = new float[3][n1];

    // Correlate for min lag.
    LocalCorrelationFilter lcf = _lcfSimple;
    lcf.setInputs(f,g);
    int lag1 = min;
    lcf.correlate(lag1,cabc[1]);
    lcf.normalize(lag1,cabc[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[] ca = (lag>min)?cabc[(i  )%3]:cabc[(i+2)%3];
      float[] cb =           cabc[(i+1)%3];
      float[] cc = (lag<max)?cabc[(i+2)%3]:cabc[(i  )%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = lag+1;
        lcf.correlate(lag1,cc);
        lcf.normalize(lag1,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, interpolate three correlation values with a quadratic
      // and then update parameters for the highest and 2nd highest peaks.
      for (int i1=0; i1<n1; ++i1) {
        float ai = ca[i1];
        float bi = cb[i1];
        float ci = cc[i1];
        if (bi>=ai && bi>=ci) {
          double q0 = bi;
          double q1 = 0.5*(ci-ai);
          double q2 = 0.5*(ci+ai)-bi;
          double up = (q2<0.0)?-0.5*q1/q2:0.0;
          double qp = q0+up*(q1+up*q2);
          if (qp>0.0) {
            float cp = (float)qp;
            if (cp>c1[i1]) {
              c2[i1] = c1[i1];
              c1[i1] = cp;
              u1[i1] = (float)(lag+up);
            } else if (cp>c2[i1]) {
              c2[i1] = cp;
            }
          }
        }
      }
    }

    // Use parameters for correlation peaks to construct an SPD 
    // tridiagonal system of equations to solve for smooth shifts.
    float s = _sigma1;
    float tiny = 0.001f;
    float cpow = 2.0f;
    float[] a = cabc[0];
    float[] b = cabc[1];
    float[] c = cabc[2];
    float[] r = u1;
    float cc = pow(c1[0]-c2[0],cpow);
    float sm = 0.0f;
    float sp = s*(1.0f-cc);
    b[0] = sp+cc+tiny;
    c[0] = -sp;
    r[0] = cc*u1[0];
    for (int i1=1; i1<n1-1; ++i1) {
      cc = pow(c1[i1]-c2[i1],cpow);
      sm = sp;
      sp = s*(1.0f-cc);
      a[i1] = -sm;
      b[i1] = sm+sp+cc+tiny;
      c[i1] = -sp;
      r[i1] = cc*u1[i1];
    }
    cc = pow(c1[n1-1]-c2[n1-1],cpow);
    sm = sp;
    a[n1-1] = -sm;
    b[n1-1] = sm+cc+tiny;
    r[n1-1] = cc*u1[n1-1];
    TridiagonalFMatrix tri = new TridiagonalFMatrix(n1,a,b,c);
    tri.solve(r,u);
  }

  private void findShiftsSmoothVersion1(
    int min, int max, float[] f, float[] g, float[] u) 
  {
    int n1 = f.length;

    // Parameters for highest correlation peaks.
    float[] b1 = new float[n1];
    float[] c1 = new float[n1];
    float[] u1 = new float[n1];

    // Parameters for second highest correlation peaks.
    float[] b2 = new float[n1];
    float[] c2 = new float[n1];
    float[] u2 = new float[n1];

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][] cabc = new float[3][n1];

    // Correlate for min lag.
    LocalCorrelationFilter lcf = _lcfSimple;
    lcf.setInputs(f,g);
    int lag1 = min;
    lcf.correlate(lag1,cabc[1]);
    lcf.normalize(lag1,cabc[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[] ca = (lag>min)?cabc[(i  )%3]:cabc[(i+2)%3];
      float[] cb =           cabc[(i+1)%3];
      float[] cc = (lag<max)?cabc[(i+2)%3]:cabc[(i  )%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = lag+1;
        lcf.correlate(lag1,cc);
        lcf.normalize(lag1,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, interpolate three correlation values with a quadratic
      // and then update parameters for the highest and 2nd highest peaks.
      for (int i1=0; i1<n1; ++i1) {
        float ai = ca[i1];
        float bi = cb[i1];
        float ci = cc[i1];
        if (bi>=ai && bi>=ci) {
          double q0 = bi;
          double q1 = 0.5*(ci-ai);
          double q2 = 0.5*(ci+ai)-bi;
          double up = (q2<0.0)?-0.5*q1/q2:0.0;
          double qp = q0+up*(q1+up*q2);
          if (qp>0.0) {
            float cp = (float)qp;
            if (cp>c1[i1]) {
              b2[i1] = b1[i1];
              c2[i1] = c1[i1];
              u2[i1] = u1[i1];
              b1[i1] = -(float)q2;
              c1[i1] = cp;
              u1[i1] = (float)(lag+up);
            } else if (cp>c2[i1]) {
              b2[i1] = -(float)q2;
              c2[i1] = cp;
              u2[i1] = (float)(lag+up);
            }
          }
        }
      }
    }

    // Use parameters for correlation peaks to construct an SPD 
    // tridiagonal system of equations to solve for smooth shifts.
    float tiny = 1.0e-5f;
    float cpow = 8.0f;
    float[] a = cabc[0];
    float[] b = cabc[1];
    float[] c = cabc[2];
    float[] r = u1;
    float a1 = pow(c1[0],cpow);
    float a2 = pow(c2[0],cpow);
    float w1 = a1*b1[0];
    float w2 = a2*b2[0];
    float cd = c1[0]-c2[0];
    //float cd = a1-a2;
    float sm = 0.0f;
    float sp = 1.0f-cd;
    b[0] = sp+cd*(w1+w2)+tiny;
    c[0] = -sp;
    r[0] = cd*(w1*u1[0]+w2*u2[0]);
    for (int i1=1; i1<n1-1; ++i1) {
      a1 = pow(c1[i1],cpow);
      a2 = pow(c2[i1],cpow);
      w1 = a1*b1[i1];
      w2 = a2*b2[i1];
      cd = c1[i1]-c2[i1];
      //cd = a1-a2;
      sm = sp;
      sp = 1.0f-cd;
      a[i1] = -sm;
      b[i1] = sm+sp+cd*(w1+w2)+tiny;
      c[i1] = -sp;
      r[i1] = cd*(w1*u1[i1]+w2*u2[i1]);
    }
    a1 = pow(c1[n1-1],cpow);
    a2 = pow(c2[n1-1],cpow);
    w1 = a1*b1[n1-1];
    w2 = a2*b2[n1-1];
    cd = c1[n1-1]-c2[n1-1];
    //cd = a1-a2;
    sm = sp;
    a[n1-1] = -sm;
    b[n1-1] = sm+cd*(w1+w2)+tiny;
    r[n1-1] = cd*(w1*u1[n1-1]+w2*u2[n1-1]);
    TridiagonalFMatrix tri = new TridiagonalFMatrix(n1,a,b,c);
    tri.solve(r,u);
  }

  private void findShifts(
    int min, int max, float[] f, float[] g, float[] u) 
  {
    int n1 = f.length;

    // Default shifts are zero.
    zero(u);

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][] c = new float[3][n1];

    // Array for current correlation maximum values.
    float[] cmax = new float[n1];

    // Correlate for min lag.
    LocalCorrelationFilter lcf = _lcfSimple;
    lcf.setInputs(f,g);
    int lag1 = min;
    lcf.correlate(lag1,c[1]);
    lcf.normalize(lag1,c[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[] ca = (lag>min)?c[(i  )%3]:c[(i+2)%3];
      float[] cb =           c[(i+1)%3];
      float[] cc = (lag<max)?c[(i+2)%3]:c[(i  )%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = lag+1;
        lcf.correlate(lag1,cc);
        lcf.normalize(lag1,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, update the correlation maximum value and displacement
      // using quadratic interpolation of three correlation values.
      for (int i1=0; i1<n1; ++i1) {
        float ai = ca[i1];
        float bi = cb[i1];
        float ci = cc[i1];
        if (bi>=ai && bi>=ci) {
          double c0 = bi;
          double c1 = 0.5*(ci-ai);
          double c2 = 0.5*(ci+ai)-bi;
          double up = (c2<0.0)?-0.5*c1/c2:0.0;
          double cp = c0+up*(c1+up*c2);
          if (cp>cmax[i1]) {
            cmax[i1] = (float)cp;
            u[i1] = (float)(lag+up);
          }
        }
      }
    }
  }

  private void findShifts(
    int dim, int min, int max, float[][] f, float[][] g, float[][] u) 
  {
    int n1 = f[0].length;
    int n2 = f.length;

    // Default shifts are zero.
    zero(u);

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][][] c = new float[3][n2][n1];

    // Array for current correlation maximum values.
    float[][] cmax = new float[n2][n1];

    // Correlate for min lag.
    LocalCorrelationFilter lcf = _lcfSimple;
    lcf.setInputs(f,g);
    int lag1 = (dim==1)?min:0;
    int lag2 = (dim==2)?min:0;
    lcf.correlate(lag1,lag2,c[1]);
    lcf.normalize(lag1,lag2,c[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[][] ca = (lag>min)?c[(i  )%3]:c[(i+2)%3];
      float[][] cb =           c[(i+1)%3];
      float[][] cc = (lag<max)?c[(i+2)%3]:c[(i  )%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = (dim==1)?lag+1:0;
        lag2 = (dim==2)?lag+1:0;
        lcf.correlate(lag1,lag2,cc);
        lcf.normalize(lag1,lag2,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, update the correlation maximum value and displacement
      // using quadratic interpolation of three correlation values.
      for (int i2=0; i2<n2; ++i2) {
        float[] ca2 = ca[i2];
        float[] cb2 = cb[i2];
        float[] cc2 = cc[i2];
        for (int i1=0; i1<n1; ++i1) {
          float ai = ca2[i1];
          float bi = cb2[i1];
          float ci = cc2[i1];
          if (bi>=ai && bi>=ci) {
            double c0 = bi;
            double c1 = 0.5*(ci-ai);
            double c2 = 0.5*(ci+ai)-bi;
            double up = (c2<0.0)?-0.5*c1/c2:0.0;
            double cp = c0+up*(c1+up*c2);
            if (cp>cmax[i2][i1]) {
              cmax[i2][i1] = (float)cp;
              u[i2][i1] = (float)(lag+up);
            }
          }
        }
      }
    }
  }

  private void findShifts(
    int dim, int min, int max, float[][][] f, float[][][] g, float[][][] u) 
  {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;

    // Default shifts are zero.
    zero(u);

    // Arrays to contain cross-correlations for three consecutive lags.
    float[][][][] c = new float[3][n3][n2][n1];

    // Array for current correlation maximum values.
    float[][][] cmax = new float[n3][n2][n1];

    // Correlate for min lag.
    LocalCorrelationFilter lcf = _lcfSimple;
    lcf.setInputs(f,g);
    int lag1 = (dim==1)?min:0;
    int lag2 = (dim==2)?min:0;
    int lag3 = (dim==3)?min:0;
    lcf.correlate(lag1,lag2,lag3,c[1]);
    lcf.normalize(lag1,lag2,lag3,c[1]);

    // For all lags in range [min,max], ...
    for (int lag=min; lag<=max; ++lag) {

      // Arrays ca, cb, and cc will contain three cross-correlations. For 
      // first and last lags, buffers a and c are the same. In other words, 
      // assume that correlation values are symmetric about the min and max 
      // lags scanned. This assumption enables local maxima to occur at the 
      // specified min and max lags, but forces displacements to lie within 
      // the range [min,max].
      int i = lag-min;
      float[][][] ca = (lag>min)?c[(i  )%3]:c[(i+2)%3];
      float[][][] cb =           c[(i+1)%3];
      float[][][] cc = (lag<max)?c[(i+2)%3]:c[(i  )%3];

      // Except for last lag, compute correlation for next lag in array cc.
      if (lag<max) {
        lag1 = (dim==1)?lag+1:0;
        lag2 = (dim==2)?lag+1:0;
        lag3 = (dim==3)?lag+1:0;
        lcf.correlate(lag1,lag2,lag3,cc);
        lcf.normalize(lag1,lag2,lag3,cc);
      }

      // For each sample, check for a local max correlation value. For each 
      // local max, update the correlation maximum value and displacement
      // using quadratic interpolation of three correlation values.
      for (int i3=0; i3<n3; ++i3) {
        for (int i2=0; i2<n2; ++i2) {
          float[] ca32 = ca[i3][i2];
          float[] cb32 = cb[i3][i2];
          float[] cc32 = cc[i3][i2];
          for (int i1=0; i1<n1; ++i1) {
            float ai = ca32[i1];
            float bi = cb32[i1];
            float ci = cc32[i1];
            if (bi>=ai && bi>=ci) {
              double c0 = bi;
              double c1 = 0.5*(ci-ai);
              double c2 = 0.5*(ci+ai)-bi;
              double up = (c2<0.0)?-0.5*c1/c2:0.0;
              double cp = c0+up*(c1+up*c2);
              if (cp>cmax[i3][i2][i1]) {
                cmax[i3][i2][i1] = (float)cp;
                u[i3][i2][i1] = (float)(lag+up);
              }
            }
          }
        }
      }
    }
  }
}
