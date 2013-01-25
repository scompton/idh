package util;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.IndexColorModel;
import java.io.IOException;

import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PlotPanel.Orientation;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Clips;

public class ViewerOld {
  
  public enum CMaps {
    GRAY("gray",ColorMap.GRAY),
    JET("jet",ColorMap.JET),
    HUE("hue",ColorMap.HUE),
    BWR("blue white red",ColorMap.BLUE_WHITE_RED),
    PRISM("prism",ColorMap.PRISM),
    AT("alpha test",ColorMap.JET);
    
    private CMaps(String s, IndexColorModel cim) {
      _s = s;
      _cim = cim;
    }
    
    @Override
    public String toString() {
      return _s;
    }
    
    public static CMaps fromString(String map) {
      for (CMaps cmaps : values()) {
        if (map.equals(cmaps._s))
          return cmaps;
      }
      return null;
    }
    
    private String _s;
    private IndexColorModel _cim;
  }
  
  public ViewerOld(float[][] f) {
    this(f,/*g*/null,/*p*/null,/*x1*/null,/*x2*/null,null,null,
        null,null,null,null,
        null,0.0f,0.0f,0,0,null);
  }
  
  public ViewerOld(float[][] f, float[][] g, float[] p, 
      String title, Sampling s1, Sampling s2, String s1Label, String s2Label,
      String cbar, float clipMin, float clipMax, int width, int height,
      Orientation orientation) 
  {
    this(f,g,p,null,null,title,s1,s2,s1Label,s2Label,cbar,null,clipMin,clipMax,
        width,height,orientation);
  }
  
  public ViewerOld(float[][] f, float[][] g, float[] p,
      String title, Sampling s1, Sampling s2, String s1Label, String s2Label,
      String cbar, IndexColorModel cmap, float clipMin, float clipMax,
      int width, int height, Orientation orientation)
  {
    this(f,g,p,null,null,title,s1,s2,s1Label,s2Label,cbar,cmap,clipMin,clipMax,
        width,height,orientation);
  }
  
  public ViewerOld(float[][] f, float[][] g, float[] p, float[] x1, float[] x2,
      String title, Sampling s1, Sampling s2, String s1Label, String s2Label,
      String cbar, IndexColorModel cmap, float clipMin, float clipMax,
      int width, int height, Orientation orientation) 
  {
    this(f,g,p,x1,x2,title,s1,s2,s1Label,s2Label,cbar,cmap,clipMin,clipMax,
        width,height,null,null,orientation);
  }
  
  public ViewerOld(float[][] f, float[][] g, float[] p, float[] x1, float[] x2,
      String title, Sampling s1, Sampling s2, String s1Label, String s2Label,
      String cbar, IndexColorModel cmap, float clipMin, float clipMax,
      int width, int height, double[] s1Limits, double[] s2Limits,
      Orientation orientation) 
  {
    title = (title==null)?"":title;
    s1 = (s1==null)?new Sampling(f[0].length):s1;
    s2 = (s2==null)?new Sampling(f.length   ):s2;
    s1Label = (s1Label==null)?"":s1Label;
    s2Label = (s2Label==null)?"":s2Label;
    double s1min = (s1Limits==null)?s1.getFirst():s1Limits[0];
    double s1max = (s1Limits==null)?s1.getLast( ):s1Limits[1];
    double s2min = (s2Limits==null)?s2.getFirst():s2Limits[0];
    double s2max = (s2Limits==null)?s2.getLast( ):s2Limits[1];
    orientation = (orientation==null)?Orientation.X1DOWN_X2RIGHT:orientation;
    
    init2D(
        f,g,p,x1,x2,title,s1,s2,s1Label,s2Label,cbar,cmap,clipMin,clipMax,
        width,height,s1min,s1max,s2min,s2max,orientation);
  }
  
  public ViewerOld(float[][][] f) {
    this(f,null,null,null,null,null,null,null,null,null,0.0f,0.0f,0,0,null);
  }
  
  public ViewerOld(float[][][] f, float[][] points, String title, 
      Sampling s1, Sampling s2, Sampling s3,
      String s1Label, String s2Label, String s3Label, String cbar,
      float clipMin, float clipMax, int width, int height,
      Orientation orientation) 
  {
    this(f,points,title,s1,s2,s3,s1Label,s2Label,s3Label,cbar,null,
        clipMin,clipMax,width,height,orientation);
  }
  
  public ViewerOld(float[][][] f, float[][] points, String title, 
      Sampling s1, Sampling s2, Sampling s3,
      String s1Label, String s2Label, String s3Label, String cbar,
      IndexColorModel cmap, float clipMin, float clipMax, int width, int height,
      Orientation orientation) 
  {
    this(f,points,title,s1,s2,s3,s1Label,s2Label,s3Label,cbar,null,
        clipMin,clipMax,width,height,null,null,orientation);
  }
  
  public ViewerOld(float[][][] f, float[][] points, String title, 
      Sampling s1, Sampling s2, Sampling s3,
      String s1Label, String s2Label, String s3Label, String cbar,
      IndexColorModel cmap, float clipMin, float clipMax, int width, int height,
      double[] s1Limits, double[] s2Limits, Orientation orientation) 
  {
    this(f,null,points,null,null,title,s1,s2,s3,s1Label,s2Label,s3Label,cbar,
        null,clipMin,clipMax,width,height,null,null,orientation);
  }
  
  public ViewerOld(float[][][] f, float[][][] g, float[][] p, float[][] x1,
      float[][] x2, String title, Sampling s1, Sampling s2, Sampling s3,
      String s1Label, String s2Label, String s3Label, String cbar,
      IndexColorModel cmap, float clipMin, float clipMax, int width, int height,
      double[] s1Limits, double[] s2Limits, Orientation orientation) 
  {
    title = (title==null)?"":title;
    s1 = (s1==null)?new Sampling(f[0][0].length):s1;
    s2 = (s2==null)?new Sampling(f[0].length   ):s2;
    s3 = (s3==null)?new Sampling(f.length      ):s3;
    s1Label = (s1Label==null)?"":s1Label;
    s2Label = (s2Label==null)?"":s2Label;
    s3Label = (s3Label==null)?"":s3Label;
    double s1min = (s1Limits==null)?s1.getFirst():s1Limits[0];
    double s1max = (s1Limits==null)?s1.getLast( ):s1Limits[1];
    double s2min = (s2Limits==null)?s2.getFirst():s2Limits[0];
    double s2max = (s2Limits==null)?s2.getLast( ):s2Limits[1];
    orientation = (orientation==null)?Orientation.X1DOWN_X2RIGHT:orientation;
    
    init3D(f,g,p,x1,x2,title,s1,s2,s3,s1Label,s2Label,s3Label,cbar,cmap,
        clipMin,clipMax,width,height,s1min,s1max,s2min,s2max,orientation);
  }
  
  private void init2D(float[][] f, Sampling s1, Sampling s2, Orientation o) {
    // Make initial panel.
    PlotPanel pp = new PlotPanel(o);
    PixelsView pv1 = pp.addPixels(s1,s2,f);
    
    JMenuBar menuBar = new JMenuBar();
    JMenu options = new JMenu("Options");
    addClipOptions(options,pv1,null);
    addColorOptions(options,pv1,null);
    menuBar.add(options);
    
    // Add everything to the PlotFrame, and display.
    _pf = new PlotFrame(pp);
    _pf.setJMenuBar(menuBar);
    _pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  }
  
  private PlotFrame _pf;
  
  private void init2D(
      float[][] f, float[][] g, float[] p, float[] x1, float[] x2, String title,
      Sampling s1, Sampling s2, String s1Label, String s2Label, String cbar, 
      IndexColorModel cmap, float clipMin, float clipMax, int width, int height,
      double s1min, double s1max, double s2min, double s2max,
      Orientation orientation)
  {
 // Make initial panel.
    PlotPanel pp = new PlotPanel(orientation);
    pp.setTitle(title);
    final PixelsView pv1 = pp.addPixels(s1,s2,f);
    if (clipMin!=0.0f || clipMax!=0.0f) {
      pv1.setClips(clipMin,clipMax);
    }
    if (cmap!=null) {
      pv1.setColorModel(cmap);  
    }
    pp.addColorBar(cbar);
    setHV(pp,orientation,s1Label,s2Label,s1min,s1max,s2min,s2max);
    
    // Possibly add second image and menu options for each PixelView.
    JMenuBar menuBar = new JMenuBar();
    JMenu options = new JMenu("Options");
    if (g!=null) {
      PixelsView pv2 = pp.addPixels(s1,s2,g);
      addClipOptions(options,pv1,"f");
      addColorOptions(options,pv1,"f");
      addClipOptions(options,pv2,"g");
      addColorOptions(options,pv2,"g");
      addAlphaOptions(options,pv2,"g");
    } else {
      addClipOptions(options,pv1,null);
      addColorOptions(options,pv1,null);
    }
    menuBar.add(options);
    
    // Possibly add points.
    if (p!=null) {
      PointsView pt1 = pp.addPoints(s1,p);
      pt1.setLineColor(Color.WHITE);
    }
    if (x1!=null && x2!=null) {
      PointsView pt2 = pp.addPoints(x1,x2);
      pt2.setStyle("rO");
    }
    
    // Add everything to the PlotFrame, and display.
    PlotFrame pf = new PlotFrame(pp);
    pf.setJMenuBar(menuBar);
    if (width!=0 && height!=0)
      pf.setSize(width,height);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
  }
  
  private void init3D(float[][][] f, float[][][] g, float[][] p, float[][] x1,
      float[][] x2, String title, Sampling s1, Sampling s2, Sampling s3,
      String s1Label, String s2Label, String s3Label, String cbar,
      IndexColorModel cmap, float clipMin, float clipMax, int width, int height,
      double s1min, double s1max, double s2min, double s2max,
      Orientation orientation)
  {
    // Make initial panel, displaying the middle frame.
    int n3 = s3.getCount();
    int i3 = n3/2;
    int r3 = (int)s3.getValue(i3);
    PlotPanel pp = new PlotPanel(orientation);
    final PixelsView pv1 = pp.addPixels(s1,s2,f[i3]);
    if (clipMin==0.0f && clipMax==0.0f) {
      Clips clips = new Clips(f);
      pv1.setClips(clips.getClipMin(),clips.getClipMax());
    } else {
      pv1.setClips(clipMin,clipMax);
    }
    if (cmap!=null) {
      pv1.setColorModel(cmap);  
    }
    title = title+" "+s3Label;
    pp.setTitle(title+" "+r3);
    pp.addColorBar(cbar);
    setHV(pp,orientation,s1Label,s2Label,s1min,s1max,s2min,s2max);
    
    JMenuBar menuBar = new JMenuBar();
    JMenu options = new JMenu("Options");
    PixelsView pv2 = null;
    if (g!=null) {
      pv2 = pp.addPixels(s1,s2,g[i3]);
      addClipOptions(options,pv1,"f");
      addColorOptions(options,pv1,"f");
      addClipOptions(options,pv2,"g");
      addColorOptions(options,pv2,"g");
      addAlphaOptions(options,pv2,"g");
    } else {
      addClipOptions(options,pv1,null);
      addColorOptions(options,pv1,null);
    }
    menuBar.add(options);
    
    PointsView pt1=null,pt2=null;
    if (p!=null) {
      Check.argument(p.length==f.length,"p.length==f.length");
      pt1 = pp.addPoints(s1,p[i3]);
      pt1.setLineColor(Color.WHITE);
    }
    if (x1!=null && x2!=null) {
      Check.argument(x1.length==f.length,"x1.length==f.length");
      Check.argument(x2.length==f.length,"x2.length==f.length");
      pt2 = pp.addPoints(x1[i3],x2[i3]);
      pt2.setStyle("rO"); 
    }
    
    SliderListener sl = new SliderListener(
        pp,pv1,pv2,pt1,pt2,s1,f,g,p,x1,x2,title);
    DefaultBoundedRangeModel brm = 
        new DefaultBoundedRangeModel(r3,0,(int)s3.getFirst(),(int)s3.getLast());
    JSlider slider = new JSlider(brm);
    slider.setMajorTickSpacing(n3/10);
    slider.setMinorTickSpacing(n3/150);
    slider.setPaintLabels(true);
    slider.setPaintTicks(true);
    slider.addChangeListener(sl);
    
    PlotFrame pf = new PlotFrame(pp);
    pf.add(slider,BorderLayout.SOUTH);
    pf.setJMenuBar(menuBar);
    if (width!=0 && height!=0)
      pf.setSize(width,height);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
  }
  
  private void setHV(
      PlotPanel pp, Orientation o, String s1Label, String s2Label,
      double s1min, double s1max, double s2min, double s2max)
  {
    if (o.equals(Orientation.X1DOWN_X2RIGHT)) {
      pp.setHLabel(s2Label);
      pp.setVLabel(s1Label);
      pp.setVLimits(s1min,s1max);
      pp.setHLimits(s2min,s2max);  
    } else {
      pp.setHLabel(s1Label);
      pp.setVLabel(s2Label);
      pp.setHLimits(s1min,s1max);
      pp.setVLimits(s2min,s2max);
    }
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
    
    new ViewerOld(f);
  }
  
  private void addClipOptions(
      JMenu options, final PixelsView pv, String label)
  {
    String name = (label==null)?"Change Clips":"Change Clips ("+label+")";
    JMenuItem changeClips = new JMenuItem(name);
    changeClips.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        float clipMin = pv.getClipMin();
        float clipMax = pv.getClipMax();
        new ClipFrame(clipMin,clipMax,pv);
      }
    });
    options.add(changeClips);
  }
  
  private void addAlphaOptions(
      JMenu options, final PixelsView pv, String label)
  {
    String name = (label==null)?"Change Alpha":"Change Alpha("+label+")";
    JMenuItem changeAlpha = new JMenuItem(name);
    changeAlpha.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        float a = pv.getColorModel().getAlpha(0)/255.0f;
        new AlphaFrame(a,pv);
      }
    });
    options.add(changeAlpha);
  }
  
  private void addColorOptions(
      JMenu options, final PixelsView pv, String label)
  {
    String name = (label==null)?"Change Colormap":"Change Colormap ("+label+")";
    JMenu changeCmap = new JMenu(name);
    JMenuItem gray = new JMenuItem(CMaps.GRAY.toString());
    JMenuItem jet = new JMenuItem(CMaps.JET.toString());
    JMenuItem bwr = new JMenuItem(CMaps.BWR.toString());
    JMenuItem hue = new JMenuItem(CMaps.HUE.toString());
    JMenuItem prism = new JMenuItem(CMaps.PRISM.toString());
    ChangeColorMapListener ccml = new ChangeColorMapListener(pv);
    gray.addActionListener(ccml);
    jet.addActionListener(ccml);
    bwr.addActionListener(ccml);
    hue.addActionListener(ccml);
    prism.addActionListener(ccml);
    changeCmap.add(gray);
    changeCmap.add(jet);
    changeCmap.add(bwr);
    changeCmap.add(hue);
    changeCmap.add(prism);
    options.add(changeCmap);
  }
  
  private class SliderListener implements ChangeListener {

    public SliderListener(
        PlotPanel pp, PixelsView pv1, PixelsView pv2, PointsView pt1, 
        PointsView pt2, Sampling s1, float[][][] f, float[][][] g, float[][] p,
        float[][] x1, float[][] x2, String title)
    {
      _pp = pp;
      _pv1 = pv1;
      _pv2 = pv2;
      _pt1 = pt1;
      _pt2 = pt2;
      _s1 = s1;
      _f = f;
      _g = g;
      _p = p;
      _x1 = x1;
      _x2 = x2;
      _title = title;
    }
    
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
//      _pp.removeTitle();
      _pp.setTitle(_title+" "+i3);
    }
    
    private PlotPanel _pp;
    private PixelsView _pv1;
    private PixelsView _pv2;
    private PointsView _pt1;
    private PointsView _pt2;
    private Sampling _s1;
    private float[][][] _f;
    private float[][][] _g;
    private float[][] _p;
    private float[][] _x1;
    private float[][] _x2;
    private String _title;
  }
  
  private static class ChangeColorMapListener implements ActionListener {

    public ChangeColorMapListener(PixelsView pv) {
      _pv = pv;
    }
    
    @Override
    public void actionPerformed(ActionEvent e) {
      CMaps source = CMaps.fromString(e.getActionCommand());
      float a = _pv.getColorModel().getAlpha(0)/255.0f;
      IndexColorModel icm;
      switch (source) {
        case GRAY:  icm = ColorMap.setAlpha(CMaps.GRAY._cim, a); break;
        case JET:   icm = ColorMap.setAlpha(CMaps.JET._cim,  a); break;
        case BWR:   icm = ColorMap.setAlpha(CMaps.BWR._cim,  a); break;
        case HUE:   icm = ColorMap.setAlpha(CMaps.HUE._cim,  a); break;
        case PRISM: icm = ColorMap.setAlpha(CMaps.PRISM._cim,a); break;
        default: throw new IllegalArgumentException(
            source+" is not a valid color map");
      }
      _pv.setColorModel(icm);
    }
    
    private PixelsView _pv;
  }

}
