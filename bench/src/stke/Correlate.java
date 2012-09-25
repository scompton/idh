package stke;

import gp404.Sequence;
import gp404.SequencePlot;
import het.RecursiveExponentialFilter;
import edu.mines.jtk.dsp.RecursiveGaussianFilter;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Parallel;
import static edu.mines.jtk.util.ArrayMath.*;

public class Correlate {
  
  private RecursiveGaussianFilter _rgf1; 
  private RecursiveGaussianFilter _rgf2; 
  private RecursiveGaussianFilter _rgf3; 
  private RecursiveExponentialFilter _ref1;
  private RecursiveExponentialFilter _ref2;
  private RecursiveExponentialFilter _ref3;
  
  public Correlate(double sigma1, double sigma2, double sigma3) {
//    _rgf1 = new RecursiveGaussianFilter(sigma1);
//    _rgf2 = new RecursiveGaussianFilter(sigma2);
//    _rgf3 = new RecursiveGaussianFilter(sigma3);
    _ref1 = new RecursiveExponentialFilter(sigma1);
    _ref2 = new RecursiveExponentialFilter(sigma2);
    _ref3 = new RecursiveExponentialFilter(sigma3);
  }
  
  public float[][][] compute(final float[][][] g, final float[][][] s) {
    final int n1 = g[0][0].length;
    final int n2 = g[0].length;
    final int n3 = g.length;
    final float[][][] gg = new float[n3][n2][n1];
    final float[][][] ss = new float[n3][n2][n1];
    final float[][][] gs = new float[n3][n2][n1];
    final float[][][] c  = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          gg[i3][i2][i1] = g[i3][i2][i1]*g[i3][i2][i1];
          ss[i3][i2][i1] = s[i3][i2][i1]*s[i3][i2][i1];
          gs[i3][i2][i1] = g[i3][i2][i1]*s[i3][i2][i1];
        }
      }
    }});
    _ref1.apply1(gg,gg);
    _ref2.apply2(gg,gg);
    _ref3.apply3(gg,gg);
    _ref1.apply1(ss,ss);
    _ref2.apply2(ss,ss);
    _ref3.apply3(ss,ss);
    _ref1.apply1(gs,gs);
    _ref2.apply2(gs,gs);
    _ref3.apply3(gs,gs);
//    _rgf1.apply0XX(gg,gg);
//    _rgf2.applyX0X(gg,gg);
//    _rgf3.applyXX0(gg,gg);
//    _rgf1.apply0XX(ss,ss);
//    _rgf2.applyX0X(ss,ss);
//    _rgf3.applyXX0(ss,ss);
//    _rgf1.apply0XX(gs,gs);
//    _rgf2.applyX0X(gs,gs);
//    _rgf3.applyXX0(gs,gs);
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          c[i3][i2][i1] = gs[i3][i2][i1]/sqrt(gg[i3][i2][i1]*ss[i3][i2][i1]);
        }
      }
    }
    return c;
  }
  
  public float[][] compute(float[][] g, float[][] s) {
    int n1 = g[0].length;
    int n2 = g.length;
    float[][] gg = new float[n2][n1];
    float[][] ss = new float[n2][n1];
    float[][] gs = new float[n2][n1];
    float[][] c  = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        gg[i2][i1] = g[i2][i1]*g[i2][i1];
        ss[i2][i1] = s[i2][i1]*s[i2][i1];
        gs[i2][i1] = g[i2][i1]*s[i2][i1];
      }
    }
    _ref1.apply1(gg,gg);
    _ref2.apply2(gg,gg);
    _ref1.apply1(ss,ss);
    _ref2.apply2(ss,ss);
    _ref1.apply1(gs,gs);
    _ref2.apply2(gs,gs);
//    _rgf1.apply0X(gg,gg);
//    _rgf2.applyX0(gg,gg);
//    _rgf1.apply0X(ss,ss);
//    _rgf2.applyX0(ss,ss);
//    _rgf1.apply0X(gs,gs);
//    _rgf2.applyX0(gs,gs);
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        float cc = gs[i2][i1]/sqrt(gg[i2][i1]*ss[i2][i1]);
        if (Float.isInfinite(cc))
          throw new RuntimeException("Infinite value detected at i2,i1="
              +i2+","+i1);
        if (Float.isNaN(cc))
          throw new RuntimeException("NaN detected at i2,i1="+i2+","+i1);
        c[i2][i1] = cc;
      }
    }
    return c;
  }
  
  public static void weight(final float[][][] f, final float a, final float l) {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float aa  = a*a;
    final float aaa = aa*a;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          float v = f[i3][i2][i1];
          if (v<=0)
            f[i3][i2][i1] = 0.0f;
          else if (v<a) {
            f[i3][i2][i1] = pow(v,l)*(3.0f/aa*v*v-2.0f/aaa*v*v*v);
          } else {
            f[i3][i2][i1] = pow(v,l);
          }
        }
      }
    }});
  }
  
  public static void applyWeights(final float[][][] f, final float[][][] w) {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          f[i3][i2][i1] = f[i3][i2][i1]*w[i3][i2][i1];
        }
      }
    }});
  }
  
  public static Sequence plotWeights(float a, float l) {
    int nt = 201;
    double dt = 1.0/nt;
    double ft = 0;
    float aa  = a*a;
    float aaa = aa*a;
    Sampling s = new Sampling(nt,dt,ft);
    float[] x = new float[nt];
    x[0] = 0.0f;
    for (int n=1; n<nt; n++) {
      float v = (float)s.getValue(n);
//    if (v<a) {
//    x[n] = pow(v,l)*(3.0f/aa*v*v-2.0f/aaa*v*v*v);
//  }
//    if (v<a) {
//      x[n] = pow(v,l)*(3.0f*pow(v/a,2)-2.0f*pow(v/a,3));
//    }
      if ((v*5.0f)<a) {
        x[n] = pow(v,l)*(3.0f*(pow(5.0f*v/a,2))-2.0f*pow(5.0f*v/a,3));
      }
      else {
        x[n] = pow(v,l);
      }
    }
    return new Sequence(s,x);
  }
  
  public static void main(String[] args) {
    Check.argument(args.length==1,"args.length==1"); 
    float l = Float.parseFloat(args[0]);
    float[] aparms = new float[] {0.1f,0.4f,0.7f,1.0f};
//    float[] aparms = new float[] {1f,4f,7f,10f};
    int na = aparms.length;
    String[] names = new String[na];
    Sequence[] seqs = new Sequence[na];
    for (int ia=0; ia<na; ia++) {
      names[ia] = "Weights a="+aparms[ia];
      seqs[ia] = plotWeights(aparms[ia],l); 
    }
//    float a = Float.parseFloat(args[0]);
//    float[] lparms = new float[] {0.1f,0.4f,0.7f,1.0f};
//    int nl = lparms.length;
//    String[] names = new String[nl];
//    Sequence[] seqs = new Sequence[nl];
//    for (int il=0; il<nl; il++) {
//      names[il] = "Weights b="+lparms[il];
//      seqs[il] = plotWeights(a,lparms[il]); 
//    }
    new SequencePlot(names,seqs);
  }

}