package warp;

import edu.mines.jtk.util.Almost;

public class SlopeTest {

  public static void main(String[] args) {
    oneSided();
    twoSided();
  }
  
  private static void oneSided() {
    int is,ic;
    int iemin = 100;
    for (int dg=1; dg<61; dg++) {
      for (int dir=1; dir<2; dir++) {
        if (dir==0) {
          is =  1;
          ic = -1;
        }
        else {
          is = -1;
          ic =  1;
        }
        int iemax = iemin+is*dg;
        for (int il=0; il<100; il++) {
          int kmax = 30;
          int kmin = 0;
          for (int k=kmin; k<=kmax; k++) {
            int rk = k*ic;
            computeSlopes(rk,is,dg,il,iemin,iemax);
          }  
        }
      }  
    }
  }
  
  private static void twoSided() {
    int is;
    int iemin = 100;
    for (int dg=1; dg<61; dg++) {
      for (int dir=0; dir<2; dir++) {
        if (dir==0)
          is=1;
        else
          is=-1;
        int iemax = iemin+is*dg;
        for (int il=0; il<100; il++) {
          int kmax = 30;
          int kmin = -kmax;
          while ((il+kmin)<0) kmin++;
          for (int k=kmin; k<=kmax; k++) {
            computeSlopes(k,is,dg,il,iemin,iemax);
          }  
        }
      }  
    }
  }
  
  private static void computeSlopes (
      int k, int is, int dg, int il, int iemin, int iemax) 
  {
    double slope = (double)-k*is/dg;
    double ylast = il+k;
    for (int j=iemin+is; j!=iemax+is; j+=is) {
      double y = il+(slope*(j-iemin)+k);
      float x = j;
//      System.out.format("slope=%5.4f, k=%d, y=%8.4f, x=%.1f\n",
//          slope,k,y,x);
      if (k<0)
        assert (y>ylast);
      if (k>0)
        assert (y<ylast);
      ylast = y;
//      if (j==iemax) {
//        double yl = il+(((double)-k*is/dg)*(j-iemin)+k);
//        assert(a.equal(yl,il)): "yl-il="+(yl-il)+
//        ", il="+il+", j-iemin="+(j-iemin)+", yl="+yl+", dg="+dg+
//        ", k="+k+", slope="+slope;
//      }
      
    }
  }
  
}
