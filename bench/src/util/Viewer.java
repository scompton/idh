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

public class Viewer {
  
  public enum CMaps {
    GRAY("gray",ColorMap.GRAY),
    JET("jet",ColorMap.JET),
    HUE("hue",ColorMap.HUE),
    BWR("blue white red",ColorMap.BLUE_WHITE_RED),
    PRISM("prism",ColorMap.PRISM);
    
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
  
  public Viewer(float[][] f) {
    this(f,null,null,null,null,null,null,null,null,null,0.0f,0.0f,0,0,null);
  }
  
  public Viewer(float[][] f, float[][] g, float[] points, String title, 
      Sampling s1, Sampling s2, String s1Label, String s2Label, String cbar, 
      float clipMin, float clipMax, int width, int height,
      Orientation orientation) 
  {
    this(f,g,points,title,s1,s2,s1Label,s2Label,cbar,null,clipMin,clipMax,
        width,height,orientation);
  }
  
  public Viewer(float[][] f, float[][] g, float[] points, String title, 
      Sampling s1, Sampling s2, String s1Label, String s2Label, String cbar, 
      IndexColorModel cmap, float clipMin, float clipMax, int width, int height,
      Orientation orientation) 
  {
    title = (title==null)?"":title;
    s1 = (s1==null)?new Sampling(f[0].length):s1;
    s2 = (s2==null)?new Sampling(f.length   ):s2;
    s1Label = (s1Label==null)?"":s1Label;
    s2Label = (s2Label==null)?"":s2Label;
    orientation = (orientation==null)?Orientation.X1DOWN_X2RIGHT:orientation;
    
    init2D(
        f,g,points,title,s1,s2,s1Label,s2Label,cbar,cmap,clipMin,clipMax,
        width,height,orientation);
  }
  
  public Viewer(float[][][] f) {
    this(f,null,null,null,null,null,null,null,null,null,0.0f,0.0f,0,0,null);
  }
  
  public Viewer(float[][][] f, float[][] points, String title, 
      Sampling s1, Sampling s2, Sampling s3,
      String s1Label, String s2Label, String s3Label, String cbar,
      float clipMin, float clipMax, int width, int height,
      Orientation orientation) 
  {
    title = (title==null)?"":title;
    s1 = (s1==null)?new Sampling(f[0][0].length):s1;
    s2 = (s2==null)?new Sampling(f[0].length   ):s2;
    s3 = (s3==null)?new Sampling(f.length      ):s3;
    s1Label = (s1Label==null)?"":s1Label;
    s2Label = (s2Label==null)?"":s2Label;
    s3Label = (s3Label==null)?"":s3Label;
    orientation = (orientation==null)?Orientation.X1DOWN_X2RIGHT:orientation;
    
    init3D(f,points,title,s1,s2,s3,s1Label,s2Label,s3Label,cbar,clipMin,clipMax,
        width,height,orientation);
  }
  
  private void init2D(
      float[][] f, float[][] g, float[] points, String title,
      Sampling s1, Sampling s2, String s1Label, String s2Label, String cbar, 
      IndexColorModel cmap, float clipMin, float clipMax, int width, int height,
      Orientation orientation)
  {
    
    PlotPanel pp = new PlotPanel(orientation);
    pp.setTitle(title);
    final PixelsView pv = pp.addPixels(s1,s2,f);
    if (clipMin!=0.0f || clipMax!=0.0f) {
      pv.setClips(clipMin,clipMax);
    }
    if (cmap!=null) {
      pv.setColorModel(cmap);  
    }
    if (orientation.equals(Orientation.X1DOWN_X2RIGHT)) {
      pp.setHLabel(s2Label);
      pp.setVLabel(s1Label);
    } else {
      pp.setHLabel(s1Label);
      pp.setVLabel(s2Label);
    }
    pp.addColorBar(cbar);
    if (points!=null) {
      PointsView ptv = pp.addPoints(s1,points);
      ptv.setLineColor(Color.WHITE);
    }
    
    JMenuBar menuBar = new JMenuBar();
    JMenu options = new JMenu("Options");
    addClipOptions(options,pv);
    addColorOptions(options,pv);
    menuBar.add(options);

    PlotFrame pf = new PlotFrame(pp);
    pf.setJMenuBar(menuBar);
    if (width!=0 && height!=0)
      pf.setSize(width,height);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
  }
  
  private void init3D(float[][][] f, float[][] points, String title, 
      Sampling s1, Sampling s2, Sampling s3,
      String s1Label, String s2Label, String s3Label, String cbar,
      float clipMin, float clipMax, int width, int height,
      Orientation orientation)
  {
    int n3 = s3.getCount();
    int i3 = n3/2;
    int r3 = (int)s3.getValue(i3);
    PlotPanel pp = new PlotPanel(orientation);
    final PixelsView pv = pp.addPixels(s1,s2,f[i3]);
    if (clipMin==0.0f && clipMax==0.0f) {
      Clips clips = new Clips(f);
      pv.setClips(clips.getClipMin(),clips.getClipMax());
    } else {
      pv.setClips(clipMin,clipMax);
    }
    title = title+" "+s3Label;
    pp.setTitle(title+" "+r3);
    if (orientation.equals(Orientation.X1DOWN_X2RIGHT)) {
      pp.setHLabel(s2Label);
      pp.setVLabel(s1Label);
    } else {
      pp.setHLabel(s1Label);
      pp.setVLabel(s2Label);
    }
    pp.addColorBar(cbar);
    PointsView ptv = null;
    if (points!=null) {
      Check.argument(points.length==f.length,"points.length==f.length");
      ptv = pp.addPoints(s1,points[i3]);
      ptv.setLineColor(Color.WHITE);
    }
    
    SliderListener sl = new SliderListener(pp,pv,ptv,s1,f,points,title);
    DefaultBoundedRangeModel brm = 
        new DefaultBoundedRangeModel(r3,0,(int)s3.getFirst(),(int)s3.getLast());
    JSlider slider = new JSlider(brm);
    slider.setMajorTickSpacing(n3/10);
    slider.setMinorTickSpacing(n3/150);
    slider.setPaintLabels(true);
    slider.setPaintTicks(true);
    slider.addChangeListener(sl);
    
    JMenuBar menuBar = new JMenuBar();
    JMenu options = new JMenu("Options");
    addClipOptions(options,pv);
    addColorOptions(options,pv);
    menuBar.add(options);

    PlotFrame pf = new PlotFrame(pp);
    pf.add(slider,BorderLayout.SOUTH);
    pf.setJMenuBar(menuBar);
    if (width!=0 && height!=0)
      pf.setSize(width,height);
    pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    pf.setVisible(true);
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
    
    new Viewer(f);
  }
  
  private void addClipOptions(JMenu options, final PixelsView pv) {
    JMenuItem changeClips = new JMenuItem("Change Clips");
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
  
  private void addColorOptions(JMenu options, final PixelsView pv) {
    JMenu changeCmap = new JMenu("Change Colormap");
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
        PlotPanel pp, PixelsView pv, PointsView ptv, Sampling s1, 
        float[][][] f, float[][] points, String title)
    {
      _pp = pp;
      _pv = pv;
      _ptv = ptv;
      _s1 = s1;
      _f = f;
      _points = points;
      _title = title;
    }
    
    @Override
    public void stateChanged(ChangeEvent e) {
      JSlider source = (JSlider)e.getSource();
      int i3 = (int)source.getValue();
      _pv.set(_f[i3]);
      if (_points!=null) {
        _ptv.set(_s1,_points[i3]);
      }
//      _pp.removeTitle();
      _pp.setTitle(_title+" "+i3);
    }
    
    private PlotPanel _pp;
    private PixelsView _pv;
    private PointsView _ptv;
    private Sampling _s1;
    private float[][][] _f;
    private float[][] _points;
    private String _title;
  }
  
  private static class ChangeColorMapListener implements ActionListener {

    public ChangeColorMapListener(PixelsView pv) {
      _pv = pv;
    }
    
    @Override
    public void actionPerformed(ActionEvent e) {
      CMaps source = CMaps.fromString(e.getActionCommand());
      switch (source) {
        case GRAY: _pv.setColorModel(CMaps.GRAY._cim); break;
        case JET: _pv.setColorModel(CMaps.JET._cim); break;
        case BWR: _pv.setColorModel(CMaps.BWR._cim); break;
        case HUE: _pv.setColorModel(CMaps.HUE._cim); break;
        case PRISM: _pv.setColorModel(CMaps.PRISM._cim); break;
      }
    }
    
    private PixelsView _pv;
  }

}
