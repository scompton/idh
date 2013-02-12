package warp;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import junit.framework.TestCase;
import junit.framework.TestSuite;

public class JTestSubsample extends TestCase {
  
  public void testSubsample1() {
    int nmax = 1500;
    int dmax = 250;
    for (int n=dmax; n<=nmax; n++) {
      for (int d=1; d<=dmax; d++) {
        try {
          Subsample.subsample(n,d);
        } catch (AssertionError e) {
          System.out.println("Failed test for n="+n+", d="+d);
          throw e;
        }
      }
    }
  }
  
  public void testSubsample2() {
    int n = 2000;
    int nm = n-1;
    float[] f = randfloat(n);
    for (int d=1; d<200; d++) {
      int g[] = Subsample.subsample(f,d);
      int ng = g.length;
      assert (g[0]==0);
      assert (g[ng-1]==nm);
      for (int ig=1; ig<ng; ig++) {
        assert((g[ig]-g[ig-1])>=d);
      }
    }
  }

  public void testSubsample3() {
    int n = 2000;
    int nm = n-1;
    int d = 200;
    float[] f = randfloat(n);
    for (int ng=2; ng<200; ng++) {
      int g[] = Subsample.subsample(f,d,ng);
      assert (g.length==ng);
      assert (g[0]==0);
      assert (g[ng-1]==nm);
    }
  }

  public void testMXB() {
    int n = 1013;
    int d = 14;
    int[] g = Subsample.subsample(n,d);
    
    float rmin = 0.5f;
    float rmax = 1.5f;
//    int kmin = (int)ceil( rmin*d);
//    int kmax = (int)floor(rmax*d);
    Map<Integer,float[][]> map = Subsample.getMXB(g,rmin,rmax,true);
    Set<Integer> keys = map.keySet();
    Iterator<Integer> it = keys.iterator();
    while (it.hasNext()) {
      float[][] mxb = map.get(it.next());
      for (float[] f : mxb) {
        int nf = f.length;
        assert (f[ 0]%1==0);
        assert (f[nf-1]==0);
      }
//      dump(mxb);
    }
  }
  
  /** 
   * This automatically generates a suite of all "test" methods
   * @return new Test
   */
  public static junit.framework.Test suite() {
    return new TestSuite(JTestSubsample.class);
  }

  /** 
   * Run all tests with text gui if this class main is invoked
   * @param args command line 
   */
  public static void main (String[] args) {
    junit.textui.TestRunner.run(suite());
  }

}
