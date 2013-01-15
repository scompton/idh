package warp;

import junit.framework.TestCase;
import junit.framework.TestSuite;

public class JTestDecompose extends  TestCase {
  
  public void test() {
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
