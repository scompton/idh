package util;

import edu.mines.jtk.util.Check;

public class Mask {

  /**
   * Returns a mask for the input array f. The mask has zero
   * values at indices that correspond to zero amplitudes in f,
   * and values of one everywhere else. This mask can be
   * applied later to any processed version of f where computed
   * values at originally zero amplitude locations should be 
   * replaced.
   * @param f the input array.
   * @return an mask array with size equal to f, containing
   *  zero values for zero amplitudes in f and values of one
   *  everywhere else.
   */
  public static float[][] getMask(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] m = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        if (f[i2][i1]==0.0f)
          m[i2][i1] = 0.0f;
        else
          m[i2][i1] = 1.0f;
      }
    }
    return m;
  }
  
  /**
   * Applies a precomputed mask m, to an input array f.
   * Amplitudes in f will zeroed at indices i2,i1 where
   * m[i2][i1]==0.
   * @param f the input array
   * @param m the mask array
   */
  public static void applyMask(float[][] f, float[][] m) {
    int n2 = f.length;
    int n1 = f[0].length;
    Check.argument(n2==m.length,"f.length==m.length");
    Check.argument(n1==m[0].length,"f[0].length==m[0].length");
    for (int i2=0; i2<n2; i2++)
      for (int i1=0; i1<n1; i1++)
        f[i2][i1] = f[i2][i1]*m[i2][i1];
  }
  
  /**
   * Uses a precomputed mask to extrapolate masked values in
   * f along the fast dimension. Where usually a mask would
   * zero the amplitudes in f at these locations, this method
   * will replace those amplitudes with the nearest unmasked 
   * value in f.  
   * @param f the input array
   * @param m the mask array
   */
  public static void extrapolate1(float[][] f, float[][] m) {
    int n2 = f.length;
    int n1 = f[0].length;
    Check.argument(n2==m.length,"f.length==m.length");
    Check.argument(n1==m[0].length,"f[0].length==m[0].length");
    for (int i2=0; i2<n2; i2++) {
      for (int i1=0; i1<n1; i1++) {
        if (m[i2][i1]==0) {
          int i1n = findIndex(m[i2],i1);
          f[i2][i1] = (i1n==Integer.MAX_VALUE)?0.0f:f[i2][i1n];
        }
      }
    }
  }
  
  /**
   * Uses a precomputed mask to extrapolate masked values in
   * f along the slow dimension. Where usually a mask would
   * zero the amplitudes in f at these locations, this method
   * will replace those amplitudes with the nearest unmasked 
   * value in f.
   * @param f the input array
   * @param m the mask array
   */
  public static void extrapolate2(float[][] f, float[][] m) {
    int n2 = f.length;
    int n1 = f[0].length;
    Check.argument(n2==m.length,"f.length==m.length");
    Check.argument(n1==m[0].length,"f[0].length==m[0].length");
    float[][] mt = transpose12(m);
    for (int i1=0; i1<n1; i1++) {
      for (int i2=0; i2<n2; i2++) {
        if (mt[i1][i2]==0) {
          int i2n = findIndex(mt[i1],i2);
          f[i2][i1] = (i2n==Integer.MAX_VALUE)?0.0f:f[i2n][i1];
        }
      }
    }
  }
  
  ///////////////////////////////////////////////////////////////////////////
  // Private
  
  /**
   * Finds the index in a mask array </code>m</code> that is
   * nearest to index </code>is</code> corresponding to an 
   * unmasked value.
   * @param m a mask array
   * @param is the starting index for the search
   * @return the index of the unmasked value closest to is,
   *  or {@link Integer#MAX_VALUE} if all of m contains no 
   *  unmasked values.
   */
  private static int findIndex(float[] m, int is) {
    int n1 = m.length;
    int ibw = -1;
    for (int i1=is-1; i1>=0; i1--) {
      if (m[i1]==1.0f) {
        ibw = i1;
        break;
      }
    }
    int ifw = -1;
    for (int i1=is+1; i1<n1; i1++) {
      if (m[i1]==1.0f) {
        ifw = i1;
        break;
      }
    }
    int dbw = (ibw<0)?Integer.MAX_VALUE:is-ibw; 
    int dfw = (ifw<0)?Integer.MAX_VALUE:ifw-is;
    return (dbw<dfw)?ibw:ifw;
  }
  
  /**
   * Transpose fast and slow dimensions of an array f.
   * @param f the input array
   * @return a transposed array
   */
  private static float[][] transpose12(float[][] f) {
    int n2 = f.length;
    int n1 = f[0].length;
    float[][] ft = new float[n1][n2];
    for (int i1=0; i1<n1; i1++)
      for (int i2=0; i2<n2; i2++)
        ft[i1][i2] = f[i2][i1];
    return ft;
  }

}
