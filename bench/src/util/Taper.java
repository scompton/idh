package util;

import static edu.mines.jtk.util.ArrayMath.*;

import javax.swing.SwingUtilities;

import viewer.Viewer2D;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.mosaic.PixelsView;

public class Taper {

  public static void taperConstant(float[][] f, int wl, int wr) {
    int n2 = f.length;
    
    while ((wl+1)>(n2-1)) wl--;
    wl = max(0,wl);
    float[] tl = copy(f[wl]);
    for (int i2=0; i2<=wl; i2++)
      f[i2] = copy(tl);
    
    while ((n2-1-wr-1)<0) wr++;
    wr = min(n2-1,wr);
    float[] tr = copy(f[n2-1-wr]);
    for (int i2=n2-1; i2>=n2-1-wr; i2--)
        f[i2] = copy(tr);
  }
  
  public static void taperConstant(float[][][] f, int wl, int wr) {
    int n3 = f.length;
    int n2 = f[0].length;
    
    while ((wl+1)>(n2-1)) wl--;
    wl = max(0,wl);
    for (int i3=0; i3<n3; i3++) {
      float[] tl = copy(f[i3][wl]);
      for (int i2=0; i2<=wl; i2++)
        f[i3][i2] = copy(tl);  
    }
    
    while ((n2-1-wr-1)<0) wr++;
    wr = min(n2-1,wr);
    for (int i3=0; i3<n3; i3++) {
      float[] tr = copy(f[i3][n2-1-wr]);
      for (int i2=n2-1; i2>=n2-1-wr; i2--)
          f[i3][i2] = copy(tr);  
    }
  }
  
  public static void main(String[] args) {
    int n2 = 50;
    int n1=100;
    final float[][] f = new float[n2][n1];
    for (int i2=0; i2<n2; i2++) {
      if (i2<7 || i2>47)
        fill(0.0f,f[i2]);
      else
        fill(i2,f[i2]);
    }
    final float[][] ft = copy(f);
    Taper.taperConstant(ft,8,3);
    SwingUtilities.invokeLater(new Runnable() {
      @Override
      public void run() {
        Viewer2D v = new Viewer2D();
        v.setTitle("Input");
        PixelsView pv1 = v.addPixels(f,"f");
        pv1.setColorModel(ColorMap.JET);
        v.addColorBar(null);
        v.show();
        Viewer2D vt = new Viewer2D();
        vt.setTitle("Tapered");
        PixelsView pv2 = vt.addPixels(ft,"ft");
        pv2.setColorModel(ColorMap.JET);
        vt.addColorBar(null);
        vt.show();
      }
    });
  }
}
