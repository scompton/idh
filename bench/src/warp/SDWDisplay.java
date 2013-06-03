package warp;

import javax.swing.JFrame;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotPanel.Orientation;
import edu.mines.jtk.mosaic.PointsView;
import static edu.mines.jtk.util.ArrayMath.*;
import viewer.Viewer2D;

public class SDWDisplay {
  
  private static final int WIDTH = 650;
  private static final int HEIGHT = 900;

  public static void display2(
      float[][] pp, float[][] ps, float[][] u, float[][] psw, float[][] vpvs,
      float[][][] e, int[][] g1, int[]g2, String desc) 
  {
    if (desc==null) desc = "";

    if (pp!=null) {
      Viewer2D v1 = new Viewer2D();
      PixelsView pv1 = v1.addPixels(pp,"pp");
      pv1.setPercentiles(1.0f,98.0f);
      if (g1!=null && g2!=null) {
        float[][][] x1x2 = getSparseGridCoords(g1,g2);
        PointsView pt = v1.addPoints2(x1x2[0],x1x2[1],"sc");
        pt.setStyle("rO");
        pt.setMarkSize(8.0f);
        Sampling s1 = v1.getSampling1();
        Sampling s2 = v1.getSampling2();
        v1.setVLimits(s1.getFirst(),s1.getLast());
        v1.setHLimits(s2.getFirst(),s2.getLast());
      }
      v1.setTitle("PP "+desc);
      v1.setHLabel("Inline");
      v1.setVLabel("Time sampling");
      v1.addColorBar("Amplitude");
      v1.setSize(WIDTH,HEIGHT);
      v1.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
      v1.show();
    }

    if (ps!=null) {
      Viewer2D v2 = new Viewer2D();
      PixelsView pv2 = v2.addPixels(ps,"ps");
      pv2.setPercentiles(1.0f,98.0f);
      v2.setTitle("PS "+desc);
      v2.setHLabel("Inline");
      v2.setVLabel("Time sampling");
      v2.addColorBar("Amplitude");
      v2.setSize(WIDTH,HEIGHT);
      v2.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
      v2.show();
    }

    if (u!=null) {
      Viewer2D v3 = new Viewer2D();
      PixelsView pv3 = v3.addPixels(u,"u");
      pv3.setColorModel(ColorMap.JET);
      v3.setHLabel("Inline");
      v3.setVLabel("Time sampling");
      v3.addColorBar("Time shift");
      v3.setTitle("U "+desc);
      v3.setSize(WIDTH,HEIGHT);
      v3.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
      v3.show();  
    }

    if (psw!=null) {
      Viewer2D v4 = new Viewer2D();
      PixelsView pv4 = v4.addPixels(psw,"psw");
      pv4.setPercentiles(1.0f,98.0f);
      v4.setTitle("PS warped "+desc);
      v4.setHLabel("Inline");
      v4.setVLabel("Time sampling");
      v4.addColorBar("Amplitude");
      v4.setSize(WIDTH,HEIGHT);
      v4.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
      v4.show();  
    }

    if (vpvs!=null) {
      Viewer2D v5 = new Viewer2D();
      PixelsView pv5 = v5.addPixels(vpvs,"vpvs");
      pv5.setColorModel(ColorMap.JET);
      v5.setTitle("Vp/Vs "+desc);
      v5.setHLabel("Inline");
      v5.setVLabel("Time sampling");
      v5.addColorBar("Interval Vp/Vs");
      v5.setSize(WIDTH,HEIGHT);
      v5.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
      v5.show();  
    }

    if (e!=null) {
      Viewer2D v6 = new Viewer2D(Orientation.X1RIGHT_X2UP);
      print("e[0][0].length="+e[0][0].length+", u[0].length="+u[0].length);
      PixelsView pv6 = v6.addPixels(e,"e");
      pv6.setClips(0.0f,0.2f);
      if (u!=null && g1!=null && g2!=null) {
        float[][][] x1x2 = getShiftCoords(u,g1,g2);
        PointsView pt1 = v6.addPoints(x1x2[0],x1x2[1],"uc");
        pt1.setStyle("gO");
        PointsView pt2 = v6.addPoints(u,"u");
        pt2.setStyle("w-");
        Sampling s1 = v6.getSampling1();
        Sampling s2 = v6.getSampling2();
        v6.setHLimits(s1.getFirst(),s1.getLast());
        v6.setVLimits(s2.getFirst(),s2.getLast());
      }
      v6.setTitle("Alignment errors "+desc);
      v6.setHLabel("Time sampling");
      v6.setVLabel("Time shift");
      v6.addColorBar("Time shift");
      v6.setSize(HEIGHT,WIDTH);
      v6.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
      v6.show();
    }

  }

  public static void displayFlat(float[][] flat) {
    Viewer2D v1 = new Viewer2D();
    PixelsView pv1 = v1.addPixels(flat,"flat");
    pv1.setPercentiles(1.0f,98.0f);
    v1.setTitle("Flat");
    v1.setHLabel("Inline");
    v1.setVLabel("Time sampling");
    v1.addColorBar("Amplitude");
    v1.setSize(WIDTH,HEIGHT);
    v1.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    v1.show();
  }

  public static void displayLSF2(float[][] p2, float[][] el) {
    Viewer2D v1 = new Viewer2D();
    PixelsView pv1 = v1.addPixels(p2,"p2");
    pv1.setPercentiles(1.0f,98.0f);
    pv1.setColorModel(ColorMap.BLUE_WHITE_RED);
    v1.setTitle("P2");
    v1.setHLabel("Inline");
    v1.setVLabel("Time sampling");
    v1.addColorBar("Slope");
    v1.setSize(WIDTH,HEIGHT);
    v1.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    v1.show();

    Viewer2D v2 = new Viewer2D();
    PixelsView pv2 = v2.addPixels(el,"el");
    pv2.setPercentiles(1.0f,98.0f);
    pv2.setColorModel(ColorMap.JET);
    v2.setTitle("Linearity");
    v2.setHLabel("Inline");
    v2.setVLabel("Time sampling");
    v2.addColorBar("Linearity");
    v2.setSize(WIDTH,HEIGHT);
    v2.getViewerFrame().setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    v2.show();
  }

  private static float[][][] getShiftCoords(float[][] u, int[][] g1, int[] g2) {
    int n2 = g1.length;
    int ng1 = g1[0].length;
    int ng2 = g2.length;
    float[][] x1 = fillfloat(-10,ng1,n2);
    float[][] x2 = fillfloat(-10,ng1,n2);
    for (int i2=0; i2<ng2; i2++) {
      int g2i = g2[i2];
      for (int i1=0; i1<ng1; i1++) {
        x1[g2i][i1] = g1[g2i][i1];
        x2[g2i][i1] = u[g2i][g1[g2i][i1]];
      }
    }
    return new float[][][] {x1,x2};
  }

  private static float[][][] getSparseGridCoords(int[][] g1, int[] g2) {
    int ng1 = g1[0].length;
    int ng2 = g2.length;
    float[][] x1 = new float[ng2][ng1];
    float[][] x2 = new float[ng2][ng1];
    for (int i2=0; i2<ng2; i2++) {
      for (int i1=0; i1<ng1; i1++) {
        x1[i2][i1] = g1[g2[i2]][i1];
        x2[i2][i1] = g2[i2];
      }
    }
    return new float[][][] {x1,x2};
  }

  private static void print(String s) {
    System.out.println(s);
  }
}
