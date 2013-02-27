package viewer;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.image.IndexColorModel;
import java.io.IOException;

import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JSlider;
import javax.swing.KeyStroke;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PlotPanel.Orientation;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.TiledView;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Clips;

/**
 * Wraps a PlotPanel with one Tile for 2D or 3D pixels and 
 * includes convenient options for interactively editing the plot.
 * This class inlcludes methods to overlay points or pixels on
 * a base PixelsView.
 * </p>
 * For 3D arrays the view is 2-Dimensional, but a slider is added 
 * to the frame to allow panning through the volume.
 */
public class Viewer {

  public Viewer(float[][] f) {
    this(f,null);
  }
  
  public Viewer(float[][] f, Orientation o) {
    this(null,null,f,o);
  }
  
  public Viewer(Sampling s1, Sampling s2, float[][] f) {
    this(s1,s2,f,null);
  }
  
  public Viewer(
      Sampling s1, Sampling s2, float[][] f, Orientation orientation)
  {
    _s1 = (s1==null)?new Sampling(f[0].length):s1;
    _s2 = (s2==null)?new Sampling(f.length   ):s2;
    orientation  = (orientation==null )?Orientation.X1DOWN_X2RIGHT:orientation;
    
    // Make initial panel.
    _pp = new PlotPanel(orientation);
    _pv1 = _pp.addPixels(_s1,_s2,f);
    _pf = new ViewerFrame(_pp,new PixelsView[]{_pv1});
  }
  
  public Viewer(double[][] f) {
    this(f,null);
  }
  
  public Viewer(double[][] f, Orientation o) {
    this(null,null,f,o);
  }
  
  public Viewer(Sampling s1, Sampling s2, double[][] f) {
    this(s1,s2,f,null);
  }
  
  public Viewer(Sampling s1, Sampling s2, double[][] f, Orientation o) {
    this(s1,s2,convertToFloat(f),o);
  }
  
  public Viewer(float[][][] f) {
    this(f,null);
  }
  
  public Viewer(float[][][] f, Orientation o) {
    this(null,null,null,f,o);
  }
  
  public Viewer(Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    this(s1,s2,s3,f,null);
  }
  
  public Viewer(
      Sampling s1, Sampling s2, Sampling s3, float[][][] f, Orientation o)
  {
    _f = f;
    _s1 = (s1==null)?new Sampling(f[0][0].length):s1;
    _s2 = (s2==null)?new Sampling(f[0].length   ):s2;
    _s3 = (s3==null)?new Sampling(f.length      ):s3;
    o  = (o==null )?Orientation.X1DOWN_X2RIGHT:o;
    
    // Make initial panel, displaying the middle frame.
    int n3 = _s3.getCount();
    _i3 = n3/2;
    int r3 = (int)_s3.getValue(_i3);
    _pp = new PlotPanel(o);
    _title = String.valueOf(_i3);
    _pp.setTitle(_title);
    _pv1 = _pp.addPixels(_s1,_s2,f[_i3]);
    Clips clips = new Clips(f);
    _pv1.setClips(clips.getClipMin(),clips.getClipMax());
    
    SliderListener sl = new SliderListener();
    //TODO Fix sampling usage here. The slider need an int so it 
    // doesn't make sense to use sampling values.
    DefaultBoundedRangeModel brm = new DefaultBoundedRangeModel(
        r3,0,(int)_s3.getFirst(),(int)_s3.getLast());
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
    _pf.addPixelOverlay(new PixelsView[]{_pv2},"2");
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
    _pf.addPixelOverlay(new PixelsView[]{_pv2},"2");
  }
  
  public void addPoints(float[] x2) {
    Check.argument(_s1.getCount()==x2.length,
        "x2.length is not consistend with sampling");
    _pt1 = _pp.addPoints(_s1,x2);
    _pt1.setLineColor(Color.WHITE);
    _pf.addPointsOverlay(_pt1,"pt1","1");
  }
  
  public void addPoints(float[][] x2) {
    Check.argument(_s1.getCount()==x2[0].length,
        "x2.length is not consistend with sampling");
    _p = x2;
    _pt1 = _pp.addPoints(_s1,x2[_i3]);
    _pt1.setLineColor(Color.WHITE);
    _pf.addPointsOverlay(_pt1,"pt1","1");
  }
  
  public void addPoints(float[] x1, float[] x2) {
    _pt2 = _pp.addPoints(x1,x2);
    _pt2.setStyle("rO");
    _pf.addPointsOverlay(_pt2,"pt2","2");
  }
  
  public void addPoints2(float[][] x1, float[][] x2) {
    _pt2 = _pp.addPoints(x1,x2);
    _pt2.setStyle("rO");
    _pf.addPointsOverlay(_pt2,"pt2","2");
  }
  
  public void addPoints3(float[][][] x1, float[][][] x2) {
    _x13 = x1;
    _x23 = x2;
    _pt2 = _pp.addPoints(x1[_i3],x2[_i3]);
    _pt2.setStyle("rO");
    _pf.addPointsOverlay(_pt2,"pt2","2");
  }
  
  public void addPoints(float[][] x1, float[][] x2) {
    Check.argument(x1.length==_f.length,"x1.length==f.length");
    Check.argument(x2.length==_f.length,"x2.length==f.length");
    _x1 = x1;
    _x2 = x2;
    _pt2 = _pp.addPoints(x1[_i3],x2[_i3]);
    _pt2.setStyle("rO");
    _pf.addPointsOverlay(_pt2,"pt2","2");
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
  
  public void addColorBar(String label) {
    _pp.addColorBar(label);
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
    Viewer v = new Viewer(f);
    v.show();
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private
  
  private ViewerFrame _pf;
  private PlotPanel _pp;
  private PixelsView _pv1;
  private PixelsView _pv2;
  private PointsView _pt1;
  private PointsView _pt2;
  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
  private JMenu _options;
  private String _title = "";
  private float[][][] _f;
  private float[][][] _g;
  private float[][] _p;
  private float[][] _x1;
  private float[][] _x2;
  private float[][][] _x13;
  private float[][][] _x23;
  private int _i3 = Integer.MIN_VALUE;

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
  
  private static void addColorKeyListener(PointsView ptv, PlotPanel pp) {
    pp.addKeyListener(new ColorKeyListener(ptv));
  }
  private static class ColorKeyListener implements KeyListener {
    ColorKeyListener(PointsView ptv) {
      _ptv = ptv;
    }

    static enum PtColor {
      R("r",Color.RED),
      B("b",Color.BLUE),
      G("g",Color.GREEN);
      
      PtColor(String k, Color c) {
        _k = k;
        _c = c;
      }
      
      public static PtColor fromString(String k) {
        for (PtColor pc : values()) {
          if (k.equals(pc._k))
            return pc;
        }
        return null;
      }
      
      String _k;
      Color _c;
    }
    @Override
    public void keyPressed(KeyEvent e) {
      String k = KeyEvent.getKeyText(e.getKeyCode());
      PtColor pc = PtColor.fromString(k); 
      _ptv.setLineColor(pc._c);
      _ptv.setMarkColor(pc._c);
      System.out.println("Changing colo to "+pc._k);
    }
    @Override
    public void keyReleased(KeyEvent e) {}
    @Override
    public void keyTyped(KeyEvent e) {}
    
    private PointsView _ptv;
  }
  
//  private void addRemoveOptions(
//      JMenu options, final TiledView tv, String label, String key)
//  {
//    String name = "Add/Remove "+label;
//    JMenuItem addRemove = new JMenuItem(name);
//    addRemove.addActionListener(new AddRemoveListener(tv));
//    addRemove.setAccelerator(KeyStroke.getKeyStroke(key));
//    options.add(addRemove);
//  }
//  
//  private class AddRemoveListener implements ActionListener {
//    public AddRemoveListener(TiledView tv) {
//      _tv = tv;
//      _tvLive = true;
//    }
//    @Override
//    public void actionPerformed(ActionEvent e) {
//      if (_tvLive) {
//        _pp.remove(_tv);
//        _tvLive = false;
//      } else {
//        _pp.addTiledView(_tv);
//        _tvLive = true;
//      }
//    }
//    private TiledView _tv;
//    private boolean _tvLive;
//  }
  
  private class SliderListener implements ChangeListener {
    @Override
    public void stateChanged(ChangeEvent e) {
      JSlider source = (JSlider)e.getSource();
      int i3 = (int)source.getValue();
      _pv1.set(_f[i3]);
      if (_g!=null)
        _pv2.set(_g[i3]);
      if (_p!=null)
        _pt1.set(_s1,_p[i3]);
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
