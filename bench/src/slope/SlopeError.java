package slope;

import edu.mines.jtk.dsp.DynamicWarping;
import edu.mines.jtk.dsp.SincInterpolator;
import static edu.mines.jtk.util.ArrayMath.*;

public class SlopeError {

  public SlopeError(int shiftMin, int shiftMax, float d1) {
    _d1 = d1;
    _d1New = _d1/RESAMPLE;
    _lmin  = shiftMin;
    _lmax  = shiftMax;
    _limin = _lmin*RESAMPLE;
    _limax = _lmax*RESAMPLE;
    _nl    = 1+_lmax-_lmin;
    _nil   = 1+_limax-_limin;
    System.out.println("d1="+_d1+", d1New="+_d1New+", lmin/lmax"+
      _lmin+"/"+_lmax+", nl="+_nl);
  }
  
  public float[][] computeErrors(float[] f, float[] g) {
    float[][] e = zerofloat(_nil, f.length);
    computeErrors(f,g,e);
    return e;
  }
  
  public float[] backtrackReverse(float[][] d, float[][] e) {
    DynamicWarping dw = new DynamicWarping(_limin,_limax);
    float[] u = zerofloat(d.length);
    dw.backtrackReverse(d,e,u);
    float[] us = div(u, RESAMPLE);
    return us;
  }
  
  public void computeErrors(float[] f, float[] g, float[][] e) {
    int n1 = f.length;
    SincInterpolator si = new SincInterpolator();
    si.setUniformSampling(n1,_d1,0.0);
    si.setUniformSamples(g);
    int nx = n1*RESAMPLE;
    float[] x = rampfloat(0f, _d1New, nx);
    float[] y = zerofloat(nx);
    si.interpolate(nx, x, y);
    for (int i1=0; i1<n1; ++i1) {
      int illo = min(_nil-1,max(0,-_limin-i1*4));
      int ilhi = max(0,min(_nil,nx-_limin-i1*4));
      for (int il=0,j1=i1*RESAMPLE+il+_limin; il<illo; ++il,++j1)
        e[i1][il] = (j1>=0) ?
          error(f[i1],y[j1]) :
          error(f[-_lmin-(il/RESAMPLE)],y[0]);
      for (int il=illo,j1=i1*RESAMPLE+il+_limin; il<ilhi; ++il,++j1)
        e[i1][il] = error(f[i1],y[j1]);
      for (int il=ilhi,j1=i1*RESAMPLE+il+_limin; il<_nil; ++il,++j1)
        e[i1][il] = (j1<nx) ?
          error(f[i1],y[j1]) :
          error(f[n1-1-_lmin-(il/RESAMPLE)],y[nx-1]);
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // private
  
  private static final int RESAMPLE = 4;
  private int _nl;
  private int _nil;
  private int _lmax;
  private int _lmin;
  private int _limax;
  private int _limin;
  private float _d1;
  private float _d1New;
  
  private float error(float f, float g) {
    return pow(abs(f-g),2);
  }
}
