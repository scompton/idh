package warp;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import junit.framework.TestCase;
import junit.framework.TestSuite;

public class JTestDecompose extends  TestCase {
  
  public void testDecompose() {
    int nmax = 1000;
    int dmax = 250;
    for (int n=1; n<=nmax; n++) {
      for (int d=1; d<=dmax; d++) {
        try {
          Decompose.decompose(n,d);
        } catch (AssertionError e) {
          System.out.println("Failed test for n="+n+", d="+d);
          throw e;
        }
      }
    }
  }

  public void testMXB() {
    int n = 1013;
    int d = 14;
    int[] g = Decompose.decompose(n,d);
    
    float rmin = 0.5f;
    float rmax = 1.5f;
//    int kmin = (int)ceil( rmin*d);
//    int kmax = (int)floor(rmax*d);
    Map<Integer,float[][]> map = Decompose.getMXB(g,rmin,rmax,true);

    System.out.println("size="+map.size());
    Set<Integer> keys = map.keySet();
    Iterator<Integer> it = keys.iterator();
    while (it.hasNext()) {
      float[][] mxb = map.get(it.next());
      for (float[] f : mxb) {
        int nf = f.length;
        assert (f[ 0]%1==0);
        assert (f[nf-1]==0);
      }
      dump(mxb);
    }
  }
  
  /** 
   * This automatically generates a suite of all "test" methods
   * @return new Test
   */
  public static junit.framework.Test suite() {
    return new TestSuite(JTestDecompose.class);
  }

  /** 
   * Run all tests with text gui if this class main is invoked
   * @param args command line 
   */
  public static void main (String[] args) {
    junit.textui.TestRunner.run(suite());
  }

}
