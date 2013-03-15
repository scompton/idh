package viewer;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.image.IndexColorModel;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.dsp.EigenTensors2;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.ColorBar;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.TensorsView;
import edu.mines.jtk.mosaic.PixelsView.Interpolation;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PlotPanel.Orientation;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Clips;

/**
 * Wraps a PlotPanel with one Tile for 2D or 3D pixels and 
 * includes convenient options for interactively editing the plot.
 * This class inlcludes methods to overlay points or pixels on
 * a base PixelsView.
 * </p>
 * For 3D arrays a slider allows panning through the volume.
 */
public class Viewer2D {

  /**
   * Constructs a pixels view of {@code f}.
   * @param f the pixel values.
   */
  public Viewer2D(float[][] f) {
    this(f,null);
  }
  
  /**
   * Constructs a pixels view of {@code f}, with specified 
   * orientation.
   * @param f the pixel values.
   * @param orientation the orientation for the pixels.
   */
  public Viewer2D(float[][] f, Orientation orientation) {
    this(null,null,f,orientation);
  }
  
  /**
   * Constructs a pixels view of {@code f} with given samplings.
   * @param f the pixel values.
   * @parama s1 the sampling of the first dimension of f.
   * @parama s2 the sampling of the second dimension of f.
   */
  public Viewer2D(Sampling s1, Sampling s2, float[][] f) {
    this(s1,s2,f,null);
  }
  
  /**
   * Constructs a pixels view of {@code f} with given samplings 
   * and with specified orientation. 
   * @param f the pixel values.
   * @parama s1 the sampling of the first dimension of f.
   * @parama s2 the sampling of the second dimension of f.
   * @param orientation the orientation for the pixels.
   */
  public Viewer2D(
      Sampling s1, Sampling s2, float[][] f, Orientation orientation)
  {
    _s1 = (s1==null)?new Sampling(f[0].length):s1;
    _s2 = (s2==null)?new Sampling(f.length   ):s2;
    _orientation = (orientation==null )?Orientation.X1DOWN_X2RIGHT:orientation;
    
    // Make initial panel.
    _pp = new PlotPanel(_orientation);
    _pv1 = _pp.addPixels(_s1,_s2,f);
    _pf = new ViewerFrame(_pp,new PixelsView[]{_pv1});
  }
  
  /**
   * Constructs a pixels view of {@code f}.
   * @param f the pixel values.
   */
  public Viewer2D(double[][] f) {
    this(f,null);
  }
  
  /**
   * Constructs a pixels view of {@code f}, with specified 
   * orientation.
   * @param f the pixel values.
   * @param orientation the orientation for the pixels.
   */
  public Viewer2D(double[][] f, Orientation orientation) {
    this(null,null,f,orientation);
  }
  
  /**
   * Constructs a pixels view of {@code f} with given samplings.
   * @param f the pixel values.
   * @parama s1 the sampling of the first dimension of f.
   * @parama s2 the sampling of the second dimension of f.
   */
  public Viewer2D(Sampling s1, Sampling s2, double[][] f) {
    this(s1,s2,f,null);
  }
  
  /**
   * Constructs a pixels view of {@code f} with given samplings 
   * and with specified orientation. 
   * @param f the pixel values.
   * @parama s1 the sampling of the first dimension of f.
   * @parama s2 the sampling of the second dimension of f.
   * @param orientation the orientation for the pixels.
   */
  public Viewer2D(
      Sampling s1, Sampling s2, double[][] f, Orientation orientation) 
  {
    this(s1,s2,convertToFloat(f),orientation);
  }
  
  /**
   * Constructs a 2D pixels view of {@code f} with a slider
   * bar that allows slicing through the full 3D volume.
   * @param f the pixel values.
   */
  public Viewer2D(float[][][] f) {
    this(f,null);
  }
  
  /**
   * Constructs a 2D pixels view of {@code f} with a slider
   * bar that allows slicing through the full 3D volume with
   * the specified orientation.
   * @param f the pixel values.
   * @param orientation the orientation for the pixels.
   */
  public Viewer2D(float[][][] f, Orientation orientation) {
    this(null,null,null,f,orientation);
  }
  
  /**
   * Constructs a 2D pixels view of {@code f} with a slider
   * bar that allows slicing through the full 3D volume with
   * the specified samplings.
   * @param s1 sampling of the first dimension of f.
   * @param s2 sampling of the second dimension of f.
   * @param s3 sampling of the third dimension of f.
   * @param f the pixel values.
   */
  public Viewer2D(Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    this(s1,s2,s3,f,null);
  }
  
  /**
   * Constructs a 2D pixels view of {@code f} with a slider
   * bar that allows slicing through the full 3D volume with
   * the specified samplings and orientation.
   * @param s1 sampling of the first dimension of f.
   * @param s2 sampling of the second dimension of f.
   * @param s3 sampling of the third dimension of f.
   * @param f the pixel values.
   * @param orientation the orientation for the pixels.
   */
  public Viewer2D(
      Sampling s1, Sampling s2, Sampling s3, float[][][] f, 
      Orientation orientation)
  {
    _f = f;
    _s1 = (s1==null)?new Sampling(f[0][0].length):s1;
    _s2 = (s2==null)?new Sampling(f[0].length   ):s2;
    _s3 = (s3==null)?new Sampling(f.length      ):s3;
    _orientation  = (orientation==null )?Orientation.X1DOWN_X2RIGHT:orientation;
    
    // Make initial panel, displaying the middle frame.
    int n3 = _s3.getCount();
    _i3 = n3/2;
    int r3 = (int)_s3.getValue(_i3);
    _pp = new PlotPanel(_orientation);
    _title = String.valueOf(_i3);
    _pp.setTitle(_title);
    _pv1 = _pp.addPixels(_s1,_s2,f[_i3]);
    Clips clips = new Clips(f);
    _pv1.setClips(clips.getClipMin(),clips.getClipMax());
    
    SliderListener sl = new SliderListener();
    DefaultBoundedRangeModel brm = new DefaultBoundedRangeModel(
        r3,0,0,n3-1);
    JSlider slider = new JSlider(brm);
    slider.setMajorTickSpacing(n3/10);
    slider.setMinorTickSpacing(n3/150);
    slider.setPaintLabels(true);
    slider.setPaintTicks(true);
    slider.addChangeListener(sl);
      
    _pf = new ViewerFrame(_pp,new PixelsView[]{_pv1});
    _pf.add(slider,BorderLayout.SOUTH);
  }
  
  public void addPixels(float[][] f) {
    Check.argument(f.length==_s2.getCount(),
        "f.length is not consistent with sampling");
    Check.argument(f[0].length==_s1.getCount(),
        "f[0].length is not consistent with sampling");
    _pv2 = _pp.addPixels(_s1,_s2,f);
    _pf.addOptions(new PixelsView[]{_pv2},"2");
  }
  
  public void addPixels(float[][][] f) {
    Check.argument(f.length==_s3.getCount(),
        "f.length is not consistent with sampling");
    Check.argument(f[0].length==_s2.getCount(),
        "f[0].length is not consistent with sampling");
    Check.argument(f[0][0].length==_s1.getCount(),
        "f[0][0].length is not consistent with sampling");
    _g = f;
    _pv2 = _pp.addPixels(_s1,_s2,f[_i3]);
    _pf.addOptions(new PixelsView[]{_pv2},"2");
  }
  
  /**
   * Adds a view of points (x1,x2) for a sampled function x2(x1).
   * @param x2 array of x2 coordinates.
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints(float[] x2, String label) {
    Check.argument(_s1.getCount()==x2.length,
        "x2.length is not consistend with sampling");
    int size = _ptMap.size();
    PointsView pt = _pp.addPoints(_s1,x2);
    _pf.addOptions(pt,label,Integer.toString(size+1));
    return pt;
  }
  
  /**
   * Adds a view of points (x1,x2) for a sampled function x2(x1),
   * for all x3.
   * @param x2 array of x2 coordinates for each x3.
   * @param label the label for these points in the options menu,
   *  or {@code null} for no label.
   * @return the points view.
   */
  public PointsView addPoints(float[][] x2, String label) {
    Check.argument(_s1.getCount()==x2[0].length,
        "x2.length is not consistend with sampling");
    int size = _ptMap.size();
    PointsView pt = _pp.addPoints(_s1,x2[_i3]);
    _ptMap.put(pt,x2);
    updatePoints();
    _pf.addOptions(pt,label,Integer.toString(size+1));
    return pt;
  }
  
  public void addPoints(float[] x1, float[] x2) {
    addPoints(x1,x2,"rO",14.0f);
  }
  
  public void addPoints(
      float[] x1, float[] x2, String style, float size) 
  {
    _pt2 = _pp.addPoints(x1,x2);
    _pt2.setStyle(style);
    _pt2.setMarkSize(size);
    _pf.addOptions(_pt2,"pt2","2");
  }
  
  public void addPoints2(float[][] x1, float[][] x2) {
    addPoints2(x1, x2,"rO",14.0f);
  }
  
  public void addPoints2(float[][] x1, float[][] x2, String style, float size) {
    _pt2 = _pp.addPoints(x1,x2);
    _pt2.setStyle(style);
    _pt2.setMarkSize(size);
    _pf.addOptions(_pt2,"pt2","2");
  }
  
  public void addPoints3(float[][][] x1, float[][][] x2) {
    _x13 = x1;
    _x23 = x2;
    _pt2 = _pp.addPoints(x1[_i3],x2[_i3]);
    _pt2.setStyle("rO");
    _pt2.setMarkSize(14.0f);
    _pf.addOptions(_pt2,"pt2","2");
  }
  
  public void addPoints(float[][] x1, float[][] x2) {
    addPoints(x1,x2,"rO",14.0f);
  }
  
  public void addPoints(float[][] x1, float[][] x2, String style, float size) {
    Check.argument(x1.length==_f.length,"x1.length==f.length");
    Check.argument(x2.length==_f.length,"x2.length==f.length");
    _x1 = x1;
    _x2 = x2;
    _pt2 = _pp.addPoints(x1[_i3],x2[_i3]);
    _pt2.setStyle(style);
    _pt2.setMarkSize(size);
    _pf.addOptions(_pt2,"pt2","2");
  }

  public void addTensors(EigenTensors2 et2) {
    addTensors(et2,_s2.getCount()/10,Color.YELLOW,1.0f);
  }
  
  public void addTensors(EigenTensors2 et2, int ne, Color color, float width) {
    TensorsView tensorView = new TensorsView(_s1,_s2,et2);
    TensorsView.Orientation orientation;
    if (_orientation.equals(PlotPanel.Orientation.X1DOWN_X2RIGHT))
      orientation = TensorsView.Orientation.X1DOWN_X2RIGHT;
    else
      orientation = TensorsView.Orientation.X1RIGHT_X2UP;
    tensorView.setOrientation(orientation);
    tensorView.setEllipsesDisplayed(ne);
    tensorView.setLineColor(color);
    Sampling e1 = new Sampling(12,0.25,0.25);
    Sampling e2 = new Sampling(10,0.5,0.0);
    tensorView.setEllipsesDisplayed(e1,e2);
    tensorView.setLineWidth(width);
    tensorView.setScale(1);
    _pp.addTiledView(tensorView);
  }
  
  public void setInterpolation(Interpolation interpolation) {
    _pv1.setInterpolation(interpolation);
  }
  
  public void setTitle(String title) {
    _title = title;
    if (_i3!=Integer.MIN_VALUE)
      _pp.setTitle(_title+" "+_i3);
    else
      _pp.setTitle(_title);
  }

  public void setHLabel(String label) {
    _pp.setHLabel(label);
  }
  
  public void setVLabel(String label) {
    _pp.setVLabel(label);
  }
  
  public void setHLimits(double hmin, double hmax) {
    _pp.setHLimits(hmin,hmax);
  }
  
  public void setVLimits(double vmin, double vmax) {
    _pp.setVLimits(vmin,vmax);
  }
  
  public void setVFormat(String format) {
    _pp.setVFormat(format);
  }
  
  public void setHFormat(String format) {
    _pp.setHFormat(format);
  }
  
  public void setHInterval(double interval) {
    _pp.setHInterval(interval);
  }

  public void setVInterval(double interval) {
    _pp.setVInterval(interval);
  }
  
  public void setClips1(float clipMin, float clipMax) {
    _pv1.setClips(clipMin,clipMax);
  }
  
  public void setClips2(float clipMin, float clipMax) {
    if (_pv2==null)
      throw new IllegalStateException(
          "Second PixelsView has not been added to plot panel.");
    _pv2.setClips(clipMin,clipMax);
  }

  public void setColorModel1(IndexColorModel colorModel) {
    _pv1.setColorModel(colorModel);
  }
  
  public void setColorModel2(IndexColorModel colorModel) {
    _pv2.setColorModel(colorModel);
  }
  
  public ColorBar addColorBar(String label) {
    return _pp.addColorBar(label);
  }
  
  public void setColorBarWidthMinimum(int widthMinimum) {
    _pp.setColorBarWidthMinimum(widthMinimum);
  }
  
  public void setSize(int width, int height) {
    _pf.setSize(width,height);
  }
  
  public void setFontSizeForPrint(double fontSize, double plotWidth) {
    _pf.setFontSizeForPrint(fontSize,plotWidth);
  }
  
  public void setFontSizeForSlide(double fracWidth, double fracHeight) {
    _pf.setFontSizeForSlide(fracWidth, fracHeight);
  }
  
  public void paintToPng(double dpi, double win, String fileName) {
    _pf.paintToPng(dpi,win,fileName);
  }
  
  public void show() {
    _pf.setVisible(true);
  }
  
  public static void main(String[] args) throws IOException {
    if (args.length != 4) {
      System.out.println("usage: java Viewer datasetPath n1 n2 n3");
      System.exit(0);
    }
    ArrayInputStream ais = new ArrayInputStream(args[0]);
    int n1 = Integer.parseInt(args[1]);
    int n2 = Integer.parseInt(args[2]);
    int n3 = Integer.parseInt(args[3]);
    float[][][] f = new float[n3][n2][n1]; 
    ais.readFloats(f);
    ais.close();
    Viewer2D v = new Viewer2D(f);
    v.show();
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private
  
  
  private ViewerFrame _pf;
  private PlotPanel _pp;
  private PixelsView _pv1;
  private PixelsView _pv2;
  private Map<PointsView,float[][]> _ptMap = 
      new HashMap<PointsView,float[][]>(); // Map for all added points views.
  PointsView[] _pts = new PointsView[0]; // Array of points views for updating.
  private PointsView _pt2;
  private Orientation _orientation;
  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
  private String _title = "";
  private float[][][] _f;
  private float[][][] _g;
  private float[][] _x1;
  private float[][] _x2;
  private float[][][] _x13;
  private float[][][] _x23;
  private int _i3 = Integer.MIN_VALUE;
  
  private void updatePoints() {
    _pts = _ptMap.keySet().toArray(new PointsView[0]);
  }

  private static float[][] convertToFloat(double[][] a) {
    int n2 = a.length;
    float[][] b = new float[n2][];
    for (int i2=0; i2<n2; ++i2) {
      int n1 = a[i2].length;
      b[i2] = new float[n1];
      for (int i1=0; i1<n1; ++i1) {
        b[i2][i1] = (float)a[i2][i1];
      }
    }
    return b;
  }
  
  private class SliderListener implements ChangeListener {
    @Override
    public void stateChanged(ChangeEvent e) {
      JSlider source = (JSlider)e.getSource();
      int i3 = (int)source.getValue();
      _pv1.set(_f[i3]);
      if (_g!=null)
        _pv2.set(_g[i3]);
      for (PointsView pt : _pts) {
        pt.set(_s1,_ptMap.get(pt)[i3]);
      }
//      if (_p!=null)
//        _pt1.set(_s1,_p[i3]);
      if (_x1!=null && _x2!=null)
        _pt2.set(_x1[i3],_x2[i3]);
      if (_x13!=null && _x23!=null)
        _pt2.set(_x13[i3],_x23[i3]);
      _i3 = i3;
      _pp.removeTitle();
      _pp.setTitle(_title+" "+_i3);
      _pf.repaint();
    }
  }
  
}
