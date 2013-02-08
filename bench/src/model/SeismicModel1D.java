/****************************************************************************
Copyright (c) 2012, Colorado School of Mines and others. All rights reserved.
This program and accompanying materials are made available under the terms of
available at http://www.eclipse.org/legal/cpl-v10.html
****************************************************************************/
package model;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.Cdouble.*;
import static edu.mines.jtk.util.ArrayMath.*;

// FOR TESTING ONLY
import javax.swing.*;
import edu.mines.jtk.mosaic.*;

/**
 * Synthetic seismograms for vertically propagating plane waves 
 * in layered media. Layer parameters include depth, velocity, 
 * density and Q.
 * @author Dave Hale, Colorado School of Mines
 * @version 2012.11.01
 */
public class SeismicModel1D {

  /**
   * Source type, either isotropic or anisotropic.
   * <p>
   * Isotropic sources produce upgoing and downgoing wavefields having 
   * the same sign. For example, marine airguns produce upgoing and 
   * downgoing pressure waves with the same sign, and a buried vibrator 
   * produces upgoing and downgoing waves with particle velocities having 
   * the same sign.
   * <p>
   * An anisotropic source produces upgoing and downgoing wavefields
   * with opposite sign. For example, dynamite buried on land produces
   * upgoing and downgoing waves with particle velocities having opposite 
   * signs.
   */
  public enum SourceType {
    /** Marine airgun (isotropic) */
    MARINE_AIRGUN,
    /** Vibroseis on land (isotropic) */
    LAND_VIBROSEIS,
    /** Dynamite on land (anisotropic) */
    LAND_DYNAMITE,
    /** Dynamite at sea (isotropic) */
    MARINE_DYNAMITE,
    /** Generic isotropic source */
    ISOTROPIC,
    /** Generic anisotropic source */
    ANISOTROPIC,
  };

  /**
   * Sensor type.
   * Sensors measure either pressure or particle velocity.
   */
  public enum SensorType {
    /** Geophone (measures particle velocity) */
    GEOPHONE,
    /** Hydrophone (measures pressure) */
    HYDROPHONE,
    /** Generic pressure sensor. */
    PRESSURE,
    /** Generic particle-velocity sensor. */
    VELOCITY,
  };

  /**
   * Sets the type of source.
   * The default source type is isotropic.
   * @param sourceType source type.
   */
  public void setSourceType(SourceType sourceType) {
    if (sourceType==SourceType.MARINE_AIRGUN) {
      sourceType = SourceType.ISOTROPIC;
    } else if (sourceType==SourceType.LAND_VIBROSEIS) {
      sourceType = SourceType.ISOTROPIC;
    } else if (sourceType==SourceType.MARINE_DYNAMITE) {
      sourceType = SourceType.ISOTROPIC;
    } else if (sourceType==SourceType.LAND_DYNAMITE) {
      sourceType = SourceType.ANISOTROPIC;
    }
    _sourceType = sourceType;
  }

  /**
   * Sets the type of sensor.
   * The default sensor type is pressure.
   * @param sensorType sensor type.
   */
  public void setSensorType(SensorType sensorType) {
    if (sensorType==SensorType.HYDROPHONE) {
      sensorType = SensorType.PRESSURE;
    } else if (sensorType==SensorType.GEOPHONE) {
      sensorType = SensorType.VELOCITY;
    }
    _sensorType = sensorType;
  }

  /**
   * Sets the absolute value of the surface reflection coefficient.
   * <p>
   * The top of the shallowest layer is assumed to be a free surface
   * with a specified reflection coefficient. For pressure seismograms,
   * the surface reflection coefficient will be negative for a wave 
   * incident from below. For particle-velocity seismograms, the surface 
   * reflection coefficient will be positive for a wave incident from
   * below.
   * <p>
   * The reflection coefficient must not be greater than 1.0.
   * The default coefficient is 0.9.
   * @param c reflection coefficient.
   */
  public void setSurfaceReflectionCoefficient(double c) {
    Check.argument(abs(c)<=1.0,"c <= 1.0");
    _surfaceRC = c;
  }

  /**
   * Sets the factor by which to oversample frequency. Oversampling 
   * reduces aliasing in time caused by frequency sampling. 
   * <p>
   * This factor must not be less than 1.0.
   * The default oversample factor is 2.0.
   * @param oversample the oversampling factor.
   */
  public void setOversample(double oversample) {
    Check.argument(oversample>=1.0,"oversample >= 1.0");
    _oversample = oversample;
  }

  /**
   * Sets factor by which ends of seismograms are exponentially decayed.
   * This decay reduces aliasing caused by sampling in frequency,
   * and will be removed after inverse Fourier transforming to time.
   * <p>
   * This factor must not be greater than 1.0.
   * The default decay factor is 0.1.
   */
  public void setDecay(double decay) {
    Check.argument(decay<=1.0,"decay <= 1.0");
    _decay = decay;
  }

  /**
   * Enables or disables use of a Ricker wavelet.
   * The default peak frequency is zero, for no Ricker wavelet.
   * @param fpeak if non-zero, peak frequency for Ricker wavelet.
   */
  public void setRickerWavelet(double fpeak) {
    _fpeak = fpeak;
  }

  /**
   * Enables or disables use of an anti-alias filter.
   * The anti-alias filter is a 12-pole causal Butterworth filter with 
   * amplitude that is 3dB down at half the Nyquist frequency.
   * @param aaf true, for anti-alias filter; false, otherwise.
   */
  public void setAntiAlias(boolean aaf) {
    _aaf = aaf;
  }

  /**
   * Adds a layer with specified depth, velocity, density and Q.
   * @param z depth of top of layer.
   * @param v layer velocity.
   * @param p layer density.
   * @param q quality factor Q.
   */
  public void addLayer(double z, double v, double p, double q) {
    // If layer already exists for the specified depth z, modify
    // the layer v, p, and q. Otherwise, insert a new layer.
    Layer layer = _layers.get(z);
    if (layer!=null) {
      layer.v = v;
      layer.p = p;
      layer.q = q;
      layer.fake = false;
    } else {
      layer = new Layer(z,v,p,q,0.0,false);
      _layers.put(z,layer);
    }
    // Update any fake layers within this layer.
    for (Double zk = _layers.higherKey(z); 
         zk!=null;
         zk = _layers.higherKey(zk)) {
      Layer lk = _layers.get(zk);
      if (lk.fake) {
        lk.v = v;
        lk.p = p;
        lk.q = q;
      } else {
        break;
      }
    }
  }

  /**
   * Adds a source with specified depth and relative amplitude.
   * @param z source depth.
   * @param a source amplitude.
   */
  public void addSource(double z, double a) {
    // If layer already exists for source at specified depth,
    // then simply modify its amplitude. Otherwise, insert a
    // new fake layer, with existing v, p, and q.
    Layer layer = _layers.get(z);
    if (layer!=null) {
      layer.a = a;
    } else {
      Double zk = _layers.floorKey(z);
      Layer lk = _layers.get(zk);
      layer = new Layer(z,lk.v,lk.p,lk.q,a,true);
      _layers.put(z,layer);
    }
  }

  /**
   * Adds a sensor with specified depth.
   * @param z sensor depth.
   */
  public void addSensor(double z) {
    _sensorDepths.add(z);
  }

  /**
   * Computes and returns seismograms for all sensors.
   * The number of seismograms (ns) equals the number of sensor depths.
   * Seismograms are ordered by sensor depths, from top to bottom.
   * The time of first sample is zero.
   * @param nt number of time samples.
   * @param dt time sampling interval.
   * @param fref reference frequency, in cycles per unit time.
   * @return array[ns][nt] of seismograms.
   */
  public float[][] makeSeismograms(int nt, double dt, double fref) {

    // Constants.
    Cdouble ci = new Cdouble(0.0,1.0);

    // Reference frequency in radians.
    double wref = 2.0*PI*fref;

    // Frequency sampling.
    int ntpad = (int)(_oversample*nt);
    int ntfft = FftReal.nfftFast(ntpad);
    int nw = ntfft/2+1;
    double dw = 2.0*PI/(ntfft*dt);

    // Imaginary part of complex frequency.
    double eps = -log(_decay)/(nt*dt);

    // Seismograms as functions of complex frequency.
    int nr = _sensorDepths.size();
    float[][] cs = new float[nr][2*nw];

    // Depth of surface.
    Double zs = _layers.firstKey();

    // For all complex frequencies, ...
    for (int iw=0,iwr=0,iwi=1; iw<nw; ++iw,iwr+=2,iwi+=2) {
      Cdouble cw = new Cdouble(iw*dw,eps);

      // Downgoing and upgoing waves for all layers.
      makeWaves(cw,wref);

      // At all receiver depths, evaluate wavefields.
      int ir = 0;
      for (Double zr:_sensorDepths) {

        boolean between = _layers.containsKey(zr);
        Double zabove = _layers.floorKey(zr);

        // If receiver is located exactly at a interface
        // (not including the surface interface), then
        // average wavefields above and below interface.
        if (!zr.equals(zs) && _layers.containsKey(zr)) {
          Layer l = _layers.get(zr);
          Layer u = _layers.get(_layers.lowerKey(zr));
          Cdouble cpower = ci.times(u.ck).times(l.z-u.z);
          Cdouble cexpd = exp(cpower);
          Cdouble cexpu = exp(neg(cpower));
          Cdouble csi = u.cd.times(cexpd).plus(u.cu.times(cexpu)
                       .plus(l.cd).plus(l.cu)).times(0.5);
          cs[ir][iwr] = (float)csi.r;
          cs[ir][iwi] = (float)csi.i;
        } 
        // Otherwise, if receiver is located exactly at the
        // surface or within a layer, then evaluate wavefield 
        // at receiver depth.
        else if (zr.equals(zs) || zabove!=null) {
          Layer u = _layers.get(zabove);
          Cdouble cpower = ci.times(u.ck).times(zr-u.z);
          Cdouble cexpd = exp(cpower);
          Cdouble cexpu = exp(neg(cpower));
          Cdouble csi = u.cd.times(cexpd).plus(u.cu.times(cexpu));
          cs[ir][iwr] = (float)csi.r;
          cs[ir][iwi] = (float)csi.i;
        }
        ++ir;
      }
    }

    // Butterworth anti-alias filter, if requested.
    if (_aaf) {
      int npole = 12;
      double pio2 = PI/2.0;
      double w3db = pio2/dt;
      for (int ipole=0; ipole<npole; ++ipole) {
        Cdouble cpower = new Cdouble(0.0,-(2*ipole+1)*pio2/npole);
        Cdouble cpole = exp(cpower).times(w3db);
        for (int iw=0,iwr=0,iwi=1; iw<nw; ++iw,iwr+=2,iwi+=2) {
          Cdouble cw = new Cdouble(iw*dw,eps);
          Cdouble ch = cpole.over(cpole.minus(cw));
          for (int ir=0; ir<nr; ++ir) {
            float csr = cs[ir][iwr];
            float csi = cs[ir][iwi];
            float chr = (float)ch.r;
            float chi = (float)ch.i;
            cs[ir][iwr] = csr*chr-csi*chi;
            cs[ir][iwi] = csr*chi+csi*chr;
          }
        }
      }
    }

    // Ricker wavelet, if requested.
    if (_fpeak>0.0) {
      double wpeak = 2.0*PI*_fpeak;
      double sigma = sqrt(2.0)/wpeak;
      float scale = (float)(4.0*sqrt(PI)/(dt*wpeak*wpeak*wpeak));
      for (int iw=0,iwr=0,iwi=1; iw<nw; ++iw,iwr+=2,iwi+=2) {
        Cdouble cw = new Cdouble(iw*dw,eps);
        Cdouble cws = cw.times(cw);
        Cdouble ch = cws.times(exp(neg(cws.over(wpeak*wpeak))));
        for (int ir=0; ir<nr; ++ir) {
          float csr = cs[ir][iwr];
          float csi = cs[ir][iwi];
          float chr = (float)ch.r;
          float chi = (float)ch.i;
          cs[ir][iwr] = scale*(csr*chr-csi*chi);
          cs[ir][iwi] = scale*(csr*chi+csi*chr);
        }
      }
    }

    // Transform frequency to time.
    float[][] s = new float[nr][nt];
    FftReal fft = new FftReal(ntfft);
    float[] sfft = new float[ntfft];
    float[] scal = new float[nt];
    for (int it=0; it<nt; ++it)
      scal[it] = (float)(exp(eps*it*dt)/ntfft);
    for (int ir=0; ir<nr; ++ir) {
      fft.complexToReal(-1,cs[ir],sfft);
      for (int it=0; it<nt; ++it)
        s[ir][it] = sfft[it]*scal[it];
    }
    return s;
  }

  /**
   * Prints parameters for all layers.
   */
  public void dumpLayers() {
    for (Double z:_layers.keySet()) {
      Layer l = _layers.get(z);
      System.out.println("Layer top      = "+l.z);
      System.out.println("      velocity = "+l.v);
      System.out.println("      density  = "+l.p);
      System.out.println("      Q        = "+l.q);
      System.out.println("      source   = "+l.a);
      if (l.fake) {
        System.out.println("      is fake");
      } else {
        System.out.println("      is not fake");
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // private

  private static class Layer {
    double z; // depth at top of layer
    double v; // layer velocity
    double p; // layer density
    double q; // layer quality factor
    double a; // source amplitude at top of layer; zero, if none
    boolean fake; // true, if this layer exists merely for source at top
    Cdouble cd; // downgoing wave at top of this layer
    Cdouble cu; // upgoing wave at top of this layer
    Cdouble ck; // complex wavenumber k = w/v within this layer
    Layer(double z, double v, double p, double q, double a, boolean fake) {
      this.z = z;
      this.v = v;
      this.p = p;
      this.q = q;
      this.a = a;
      this.fake = fake;
    }
  }
  private TreeMap<Double,Layer> _layers = new TreeMap<Double,Layer>();
  private TreeSet<Double> _sensorDepths = new TreeSet<Double>();
  private double _surfaceRC = 0.9f; // absolute value of refl coeff at surface
  private double _decay = 0.1f; // determines imag part of complex frequency
  private double _oversample = 2.0; // factor by which to oversample frequency
  private SourceType _sourceType = SourceType.MARINE_AIRGUN;
  private SensorType _sensorType = SensorType.HYDROPHONE;
  private boolean _aaf; // true, for anti-alias filter
  private double _fpeak; // non-zero, for Ricker wavelet

  private void makeWaves(Cdouble cw, double wref) {

    // Constants.
    Cdouble ciw = new Cdouble(-cw.i,cw.r);
    Cdouble cmiw = ciw.neg();

    // Lower and upper depths and layers. 
    Double zl,zu;
    Layer l,u;

    // Complex wavenumber for lower halfspace.
    zl = _layers.lastKey();
    l = _layers.get(zl);
    double gamma = atan(1.0/l.q)/PI;
    Cdouble cv = pow(cmiw.over(wref),gamma).times(l.v*cos(gamma*0.5*PI));
    l.ck = cw.over(cv);

    // Complex impedance of lower halfspace.
    Cdouble cil = cv.times(l.p);

    // Initially assume unit-amplitude downgoing wave in lower halfspace.
    Cdouble c1d = new Cdouble(1.0,0.0);
    Cdouble c1u = new Cdouble(0.0,0.0);

    // Initialize source terms to zero.
    Cdouble csd = new Cdouble(0.0,0.0);
    Cdouble csu = new Cdouble(0.0,0.0);

    // First loop over layers, bottom to top.
    for (zu=_layers.lowerKey(zl); zu!=null; zl=zu,zu=_layers.lowerKey(zl)) {
      l = _layers.get(zl);
      u = _layers.get(zu);

      // Complex velocity, wavenumber, impedance, and reflection coefficient.
      gamma = atan(1.0/l.q)/PI;
      cv = pow(cmiw.over(wref),gamma).times(u.v*cos(gamma*0.5*PI));
      u.ck = cw.over(cv);
      Cdouble ciu = cv.times(u.p);
      Cdouble cr = ciu.minus(cil).over(ciu.plus(cil));
      if (_sensorType==SensorType.PRESSURE) 
        cr.negEquals();

      // Complex inverse transmission coefficient.
      Cdouble cit = inv(cr.plus(1.0));

      // Complex propagation matrix.
      Cdouble cpower = ciw.over(cv).times(l.z-u.z);
      Cdouble cplus = cit.times(exp(cpower));
      Cdouble cminus = cit.times(exp(neg(cpower)));
      Cdouble ca11 = cminus;
      Cdouble ca12 = cr.times(cminus);
      Cdouble ca21 = cr.times(cplus);
      Cdouble ca22 = cplus;

      // Downgoing and upgoing waves, without source terms.
      Cdouble ctemp = c1d;
      c1d = ca11.times(ctemp).plus(ca12.times(c1u));
      c1u = ca21.times(ctemp).plus(ca22.times(c1u));

      // Complex source terms.
      ctemp = csd;
      csd = ca11.times(ctemp).plus(ca12.times(csu));
      csu = ca21.times(ctemp).plus(ca22.times(csu));
      if (l.a!=0.0) {
        csd.minusEquals(cminus.times(l.a));
        if (_sourceType==SourceType.ISOTROPIC) {
          csu.plusEquals(cplus.times(l.a));
        } else {
          csu.minusEquals(cplus.times(cr.times(2.0).plus(1.0)).times(l.a));
        }
      }

      // Upper layer becomes lower layer in next iteration.
      cil = ciu;
    }

    // Reflection coefficient at surface.
    double rs = abs(_surfaceRC);
    if (_sensorType==SensorType.PRESSURE)
      rs = -rs;

    // Amplitude of source at surface.
    l = _layers.get(zl);
    double as = l.a;

    // Downgoing and upgoing wave in lower halfspace.
    zl = _layers.lastKey();
    l = _layers.get(zl);
    l.cd = csu.times(rs).minus(csd).plus(as).over(c1d.minus(c1u.times(rs)));
    l.cu = new Cdouble(0.0,0.0);

    // Complex impedance of lower halfspace.
    gamma = atan(1.0/l.q)/PI;
    cil = pow(cmiw.over(wref),gamma).times(l.p*l.v*cos(gamma*0.5*PI));

    // Second loop over layers, bottom to top.
    for (zu=_layers.lowerKey(zl); zu!=null; zl=zu,zu=_layers.lowerKey(zl)) {
      l = _layers.get(zl);
      u = _layers.get(zu);

      // Complex velocity, impedance, and reflection coefficient.
      gamma = atan(1.0/l.q)/PI;
      cv = pow(cmiw.over(wref),gamma).times(u.v*cos(gamma*0.5*PI));
      Cdouble ciu = cv.times(u.p);
      Cdouble cr = ciu.minus(cil).over(ciu.plus(cil));
      if (_sensorType==SensorType.PRESSURE) 
        cr.negEquals();

      // Complex inverse transmission coefficient.
      Cdouble cit = inv(cr.plus(1.0));

      // Complex propagation matrix.
      Cdouble cpower = ciw.over(cv).times(l.z-u.z);
      Cdouble cplus = cit.times(exp(cpower));
      Cdouble cminus = cit.times(exp(neg(cpower)));
      Cdouble ca11 = cminus;
      Cdouble ca12 = cr.times(cminus);
      Cdouble ca21 = cr.times(cplus);
      Cdouble ca22 = cplus;

      // Downgoing and upgoing waves, without source terms.
      u.cd = ca11.times(l.cd).plus(ca12.times(l.cu));
      u.cu = ca21.times(l.cd).plus(ca22.times(l.cu));

      // Add complex source terms.
      if (l.a!=0.0) {
        u.cd.minusEquals(cminus.times(l.a));
        if (_sourceType==SourceType.ISOTROPIC) {
          u.cu.plusEquals(cplus.times(l.a));
        } else {
          u.cu.minusEquals(cplus.times(cr.times(2.0).plus(1.0)).times(l.a));
        }
      }

      // Upper layer becomes lower layer in next iteration.
      cil = ciu;
    }
  }

  ///////////////////////////////////////////////////////////////////////////
  // testing

  public static void main(String[] args) {
    SwingUtilities.invokeLater(new Runnable() {
      public void run() {
        //goTest0();
        goTest1();
        //goDumpLayers();
      }
    });
  }
  private static void goTest0() {
    SeismicModel1D sm = new SeismicModel1D();
    sm.setSourceType(SeismicModel1D.SourceType.MARINE_AIRGUN);
    sm.setSensorType(SeismicModel1D.SensorType.HYDROPHONE);
    sm.setSurfaceReflectionCoefficient(0.9);
    sm.addLayer(0.0,2.0,2.0,1.0e6);
    sm.addSource(0.02,1.0);
    sm.addSensor(2.02);
    sm.dumpLayers();
    int nt = 1101;
    double dt = 0.004;
    double fref = 0.5/dt;
    boolean aaf = false;
    Sampling st = new Sampling(nt,dt,0.0);
    float[][] s = sm.makeSeismograms(nt,dt,fref);
    System.out.println("min = "+min(s[0])+" max = "+max(s[0]));
    SimplePlot sp = new SimplePlot();
    sp.addPoints(st,s[0]);
  }
  private static void goTest1() {
    SeismicModel1D sm = new SeismicModel1D();
    sm.setSourceType(SeismicModel1D.SourceType.MARINE_AIRGUN);
    sm.setSensorType(SeismicModel1D.SensorType.HYDROPHONE);
    sm.setSurfaceReflectionCoefficient(0.9);
    sm.setDecay(0.9999);
    sm.addLayer(0.0,2.0,2.0,1.0e6);
    sm.addLayer(1.0,6.0,2.0,1.0e6);
    sm.addSource(0.02,2.0);
    sm.addSensor(0.02);
    sm.dumpLayers();
    int nt = 1101;
    double dt = 0.004;
    double fref = 0.5/dt;
    boolean aaf = false;
    Sampling st = new Sampling(nt,dt,0.0);
    float[][] s = sm.makeSeismograms(nt,dt,fref);
    System.out.println("min = "+min(s[0])+" max = "+max(s[0]));
    SimplePlot sp = new SimplePlot();
    sp.addPoints(st,s[0]);
  }
  private static void goDumpLayers() {
    SeismicModel1D sm = new SeismicModel1D();
    sm.addLayer(2.0,2.0,2.0,2.0);
    sm.addLayer(5.0,5.0,5.0,5.0);
    sm.addSource(3.0,3.0);
    sm.addSource(4.0,4.0);
    sm.addLayer(1.0,1.0,1.0,1.0);
    sm.addLayer(3.0,3.0,3.0,3.0);
    sm.dumpLayers();
  }
}
