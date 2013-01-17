package warp;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.HashMap;
import java.util.Map;

public class Decompose {

  /**
   * Ensures that the sparse grid decomposition falls on n-1.
   * The interval of the sparse grid points may be greater than
   * or equal to {@code d}, to ensure that the sparse grid
   * includes n-1.
   * @return an integer array of sparse grid points.
   */ 
  public static int[] decompose(int n, int d) {
    int[] g;
    int last = n-1;
    
    if (d>last) {
      g = new int[]{0,last};
    } else if (last%d == 0) { // n is perfectly divisible by d! 
      int ng = last/d + 1;
      g = new int[ng];
      ramp(0,d,g);
    } else {// Compute grid with d, ensuring that the grid includes last.
      int ng = (int)floor(last/d)+1;
      g = new int[ng];
      g[0] = 0;
      for (int ig=1; ig<ng; ig++) {
        int gm = g[ig-1];
        int id = round((float)(last-gm)/(ng-ig));
        assert (id>=d);
//        System.out.println("id="+id);
        g[ig] = gm+id;
      }
    }
    int ng = g.length;
    try {
      assert(g[ng-1]==last) : "Expected g[ng-1] to equal "+last+
      ", but g[ng-1]="+g[ng-1];  
    } catch (AssertionError e) {
      dump(g);
      throw e;
    }
    
    return g;
  }
  
  public static Map<Integer,float[][]> getMXB(
      int[] g, float rmin, float rmax, boolean up) 
  {
    Map<Integer,float[][]> mxbMap= new HashMap<Integer,float[][]>();
    int ks = up? 1:-1;
    int ms = up?-1: 1;
    int ng = g.length;
    for (int ig=1; ig<ng; ig++) {
      int dg = g[ig]-g[ig-1];
      if (mxbMap.containsKey(dg))
        continue;
      int kmin = (int)ceil( rmin*dg);
      int kmax = (int)floor(rmax*dg);
      int nk = kmax+1;
//      System.out.println("dg="+dg+", kmin="+kmin+", kmax="+kmax+", nk="+nk);
      float[][] mx = new float[nk][dg+1];
      for (int k=kmin; k<=kmax; k++) {
        float m = (float)k/dg*ms;
        for (int x=0; x<=dg; x++) {
          mx[k][x] = m*x+ks*k;
        }  
      }
      mxbMap.put(dg,mx);
    }
    return mxbMap;
  }
  
  public static void main(String[] args) {
    if (args.length != 2) {
      System.out.println("Usage: java Decompose numberOfSamples delta");
      System.exit(0);
    }
    int n = Integer.parseInt(args[0]);
    int d = Integer.parseInt(args[1]);
    int[] g = decompose(n,d); 
    dump(g);
  }

}
