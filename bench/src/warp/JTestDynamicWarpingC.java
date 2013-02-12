package warp;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.Map;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.mosaic.PlotPanel.Orientation;
import util.Viewer;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class JTestDynamicWarpingC extends TestCase {

  /*
   * Tests sparse accumulation for the simplest case when rmin and rmax
   * are equal to zero, and all errors are 1.0. This test verifies 
   * returned values for accumulating forward, accumulating reverse,
   * and smoothing (forward + reverse - duplicated error).
   */
  public void testAccumulateSparseSimple() {
    int nl = 50;
    int n1 = 100;
    int n1m1 = n1-1;
    float[][] e = new float[n1][nl];
    fill(1.0f,e);
    float rmin = 0.0f;
    float rmax = 0.0f;
    // Test Forward - After normalization, errors should equal 0,
    // at the zero index, and 1 at the last index. Errors should always
    // be increasing from left to right.
    for (int ig=1; ig<=n1; ig++) {
      int[] g = Subsample.subsample(n1,ig);
      int ng = g.length;
      int ngm1 = ng-1;
      float[][][] dmf = DynamicWarpingC.accumulateForwardSparse(e,rmin,rmax,g);
      float[][] df = dmf[0];
      float[][] mf = dmf[1];
      float lvf = 0.0f;
//      Viewer vef = new Viewer(DynamicWarpingC.transposeLag(df),
//          Orientation.X1RIGHT_X2UP);
//      vef.setTitle("Accumulated Errors Forward");
//      vef.setSize(900,600);
//      vef.setColorModel1(ColorMap.JET);
//      vef.setClips1(0.0f,1.0f);
//      vef.addColorBar("Error");
//      vef.show();
      for (int i1=0; i1<ng; i1++) {
        for (int il=0; il<nl; il++) {
          assert(mf[i1][il]==0.0f);
          if (i1==0)
            assert(df[i1][il]==0.0f);
          if (i1==n1m1)
            assert(df[i1][il]==1.0f);
          assert(df[i1][il]>=lvf);
        }
        lvf = df[i1][0];
      }
      // Test Reverse - After normalization, errors should equal 0,
      // at the last index, and 1 at the zero index. Errors should always
      // be increasing from right to left.
      float[][][] dmr = DynamicWarpingC.accumulateReverseSparse(e,rmin,rmax,g);
      float[][] dr = dmr[0];
      float[][] mr = dmr[1];
      float lvr = 0.0f;
//      Viewer ver = new Viewer(DynamicWarpingC.transposeLag(dr),
//          Orientation.X1RIGHT_X2UP);
//      ver.setTitle("Accumulated Errors Reverse");
//      ver.setSize(900,600);
//      ver.setColorModel1(ColorMap.JET);
//      ver.setClips1(0.0f,1.0f);
//      ver.addColorBar("Error");
//      ver.show();  
      for (int i1=ngm1; i1>=0; i1--) {
        for (int il=0; il<nl; il++) {
          assert(mr[i1][il]==0.0f);
          if (i1==0)
            assert(dr[i1][il]==1.0f);
          if (i1==n1m1)
            assert(dr[i1][il]==0.0f);
          assert(dr[i1][il]>=lvr);
        }
        lvr = dr[i1][0];
      }
      // Test Smoothing - For these simple alignment errors, the smoothed
      // errors should be constant everywhere.
      Map<Integer,float[][]> fmap = Subsample.getMXB(g,rmin,rmax,false);
      Map<Integer,float[][]> rmap = Subsample.getMXB(g,rmin,rmax,true);
      float[][] es = 
          DynamicWarpingC.smoothErrorsSparse(e,rmin,rmax,g,fmap,rmap);
//      Viewer ves = new Viewer(DynamicWarpingC.transposeLag(es),
//          Orientation.X1RIGHT_X2UP);
//      ves.setTitle("Smoothed Errors");
//      ves.setSize(900,600);
//      ves.setColorModel1(ColorMap.JET);
//      ves.setClips1(0.0f,1.0f);
//      ves.addColorBar("Error");
//      ves.show();
      float sc = es[0][0];
      System.out.println(sc);
      for (int i1=0; i1<ng; i1++) {
        for (int il=0; il<nl; il++) {
          assert(es[i1][il]==sc):
            "es["+i1+"]["+il+"]="+es[i1][il]+", expected "+sc;
        }
      }
    }
  }
  
  /** 
   * This automatically generates a suite of all "test" methods
   * @return new Test
   */
  public static junit.framework.Test suite() {
    return new TestSuite(JTestDynamicWarpingC.class);
  }

  /** 
   * Run all tests with text gui if this class main is invoked
   * @param args command line 
   */
  public static void main (String[] args) {
    junit.textui.TestRunner.run(suite());
  }

}
