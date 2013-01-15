package warp;

import edu.mines.jtk.util.ArrayMath;

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
      ArrayMath.ramp(0,d,g);
    } else {// Compute grid with d, ensuring that the grid includes last.
      int ng = (int)ArrayMath.floor(last/d)+1;
      g = new int[ng];
      g[0] = 0;
      for (int ig=1; ig<ng; ig++) {
        int gm = g[ig-1];
        int id = ArrayMath.round((float)(last-gm)/(ng-ig));
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
      ArrayMath.dump(g);
      throw e;
    }
    
    return g;
  }
  
  public static void main(String[] args) {
    if (args.length != 2) {
      System.out.println("Usage: java Decompose numberOfSamples delta");
      System.exit(0);
    }
    int n = Integer.parseInt(args[0]);
    int d = Integer.parseInt(args[1]);
    int[] g = decompose(n,d); 
    ArrayMath.dump(g);
  }

}
