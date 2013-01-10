package warp;

import edu.mines.jtk.util.ArrayMath;

public class Decompose {

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
  
  public static float[][] extrapolateErrors(int d, float[][] e) {
    int n1 = e.length;
    int nl = e[0].length;
    int last = n1-1;
    float[] elast = e[last];
    while (last%d != 0) last++;
    int n1new = last+1;
    float[][] e2 = new float[n1new][nl];
    for (int i1=0; i1<n1; i1++) {
      System.arraycopy(e[i1],0,e2[i1],0,nl);
    }
    for (int i1=n1; i1<n1new; i1++) {
      System.arraycopy(elast,0,e2[i1],0,nl);
    }
    return e2;    
  }
  
  public static float[][][] extrapolateErrors1(int d1, float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    int n1m = n1-1;
    int ne1m = n1m;
    while (ne1m%d1 != 0) ne1m++;
    int ne1 = ne1m+1;
    float[][][] ee = new float[n2][ne1][nl];
    for (int i2=0; i2<n2; i2++) {
      float[] e1Last = e[i2][n1m];
      for (int i1=0; i1<n1; i1++) {
        System.arraycopy(e[i2][i1],0,ee[i2][i1],0,nl);
      }
      for (int i1=n1; i1<ne1; i1++) {
        System.arraycopy(e1Last,0,ee[i2][i1],0,nl);
      }
    }
    
    return ee;    
  }
  
  public static float[][][][] extrapolateErrors1(int d1, float[][][][] e) {
    int nl = e[0][0][0].length;
    int n1 = e[0][0].length;
    int n2 = e[0].length;
    int n3 = e.length;
    int n1m = n1-1;
    int ne1m = n1m;
    while (ne1m%d1 != 0) ne1m++;
    int ne1 = ne1m+1;
    float[][][][] ee = new float[n3][n2][ne1][nl];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        float[] e1Last = e[i3][i2][n1m];
        for (int i1=0; i1<n1; i1++) {
          System.arraycopy(e[i2][i1],0,ee[i2][i1],0,nl);
        }
        for (int i1=n1; i1<ne1; i1++) {
          System.arraycopy(e1Last,0,ee[i2][i1],0,nl);
        }
      }  
    }
    return ee;    
  }
  
  public static float[][][] extrapolateErrors2(int d2, float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    int n2m = n2-1;
    int ne2m = n2m;
    while (ne2m%d2 != 0) ne2m++;
    int ne2= ne2m+1;
    float[][][] ee = new float[ne2][n1][nl];
    float[][] e2Last = e[n2m];
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        for (int il=0; il<nl; il++) {
          ee[i2][i1][il] = e[i2][i1][il];
        }
      }
    }  
    for (int i2=n2; i2<ne2; i2++) {
      ee[i2] = e2Last;
    }
    
    return ee;    
  }
  
  public static float[][][][] extrapolateErrors2(int d2, float[][][][] e) {
    int nl = e[0][0][0].length;
    int n1 = e[0][0].length;
    int n2 = e[0].length;
    int n3 = e.length;
    int n2m = n2-1;
    int ne2m = n2m;
    while (ne2m%d2 != 0) ne2m++;
    int ne2= ne2m+1;
    float[][][][] ee = new float[n3][ne2][n1][nl];
    for (int i3=0; i3<n3; i3++) {
      float[][] e2Last = e[i3][n2m];
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          for (int il=0; il<nl; il++) {
            ee[i3][i2][i1][il] = e[i3][i2][i1][il];
          }
        }
      }  
      for (int i2=n2; i2<ne2; i2++) {
        ee[i3][i2] = e2Last;
      }
    }
    
    return ee;    
  }
  
//  public static double[] interpValues(int d) {
//    for (int k=0; k<d; k++) {
//      double slope = (double)k/d;
//    }
//  }
  
  /**
   * Expands the length of the sparse grid array so that
   * the sparse grid has uniform spacing of d. 
   * @param n the length that will be coarsely sampled by d.
   * @param d the course sampling delta
   * @return an array of sparse grid values. This array may contain
   *  values greater than n in order to maintain uniform spacing.
   */
  public static int[] decomposeUniform(int n, int d) {
    int last = n-1;
    while (last%d != 0) last++;
    int ng = last/d + 1;
    int[] g = new int[ng];
    ArrayMath.ramp(0,d,g);
    return g;
  }
  
  // Ensures that the sparse grid decomposition falls on 
  // n-1. The sparse grid has mostly d spacing, except near
  // the end.
  public static int[] decompose(int n, int d) {
    int[] g;
    int last = n-1;
    
    // n is perfectly divisible by d! 
    if (last%d == 0) {
      int ng = last/d + 1;
      g = new int[ng];
      ArrayMath.ramp(0,d,g);
    } else {// Compute grid with d, ensuring that the grid includes last.  
      int ndm = (int)ArrayMath.ceil((float)last/d);
      int dmax = (ndm-1)*d;
      int diff = last-dmax;
//      System.out.println("Start ndm="+ndm);
//      System.out.println("Start dmax="+dmax);
//      System.out.println("Start diff="+diff);
      if (diff!=(d+1) && diff!=(d-1)) {
        if (diff<d && dmax>d) {
          dmax -= d;
          diff = last-dmax;
        } else {
          diff = last-dmax;
          ndm += 1;
        }  
      } else {
        ndm += 1;
      }
      g = new int[ndm];  
      int max = ndm-1;
      for (int i=0; i<max; i++) {
        g[i] = i*d;
      }
      g[max] = g[max-1] + diff;
//      System.out.println("End ndm="+ndm);
//      System.out.println("End dmax="+dmax);
//      System.out.println("End diff:"+diff);
    }

    int ng = g.length;
    assert(g[ng-1]==last) : "Expected g[ng-1] to equal "+last+
      ", but g[ng-1]="+g[ng-1];
//    ArrayMath.dump(g);
    return g;
  }
  
}
