package viewer;

import java.awt.Component;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionAdapter;
import java.awt.event.MouseMotionListener;
import java.awt.image.IndexColorModel;
import java.io.IOException;

import javax.swing.KeyStroke;

import edu.mines.jtk.awt.Mode;
import edu.mines.jtk.awt.ModeManager;
import edu.mines.jtk.awt.ModeMenuItem;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PlotPanelPixels3.Orientation;
import edu.mines.jtk.mosaic.Tile;

import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.mosaic.PlotPanelPixels3.AxesPlacement;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.util.Check;

/**
 * Wraps PlotPanelPixels3 with some convenient options. 
 */
public class Viewer3 {

  /**
   * Construct a 3 panel view of f, with default orientation,
   * {@link Orientation#X1DOWN_X2RIGHT} and default Samplings.
   * @param f
   */
  public Viewer3(float[][][] f) {
    this(f,null);
  }
  
  /**
   * Construct a 3 panel view of f, with specified orientation
   * and default Samplings.
   * @param f
   */
  public Viewer3(float[][][] f, Orientation o) {
    this(null,null,null,f,o);
  }
  
  /**
   * Construct a 3 panel view of f, with default orientation,
   * {@link Orientation#X1DOWN_X2RIGHT} and specified Samplings.
   * @param f
   */
  public Viewer3(Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    this(s1,s2,s3,f,null);
  }
  
  /**
   * Construct a 3 panel view of f, with specified orientation
   * and Samplings.
   * @param f
   */
  public Viewer3(
      Sampling s1, Sampling s2, Sampling s3, float[][][] f, 
      Orientation orientation)
  {
    _f = f;
    _s1 = (s1==null)?new Sampling(f[0][0].length):s1;
    _s2 = (s2==null)?new Sampling(f[0].length   ):s2;
    _s3 = (s3==null)?new Sampling(f.length      ):s3;
    _orientation = (orientation==null )?Orientation.X1DOWN_X2RIGHT:orientation;
    
    // Make initial panel, displaying the middle frame.
    _pp = new PlotPanelPixels3(_orientation,AxesPlacement.LEFT_BOTTOM,
        _s1,_s2,_s3,f);
    _pp.getMosaic().setHeightElastic(0,100);
    _pp.getMosaic().setHeightElastic(1,200);
    _pv1 = new PixelsView[]{
        _pp.getPixelsView12(),
        _pp.getPixelsView13(),
        _pp.getPixelsView23()
    };
    _vf = new ViewerFrame(_pp,_pv1);
    SliceMode sm = new SliceMode(_vf.getModeManager());
    _vf.addToMenu(new ModeMenuItem(sm));
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
    //TODO Fix
//    updateOptions(_pv2,"2");
  }
  
  //TODO Fix
//  public void addPoints(float[] x2) {
//    Check.argument(_s1.getCount()==x2.length,
//        "x2.length is not consistend with sampling");
//    _pt1 = _pp.addPoints(_s1,x2);
//    _pt1.setLineColor(Color.WHITE);
//    addRemoveOptions(_options,_pt1,"pt1","1");
//  }
//  
//  //TODO Fix
//  public void addPoints(float[][] x2) {
//    Check.argument(_s1.getCount()==x2[0].length,
//        "x2.length is not consistend with sampling");
//    _p = x2;
//    _pt1 = _pp.addPoints(_s1,x2[_i3]);
//    _pt1.setLineColor(Color.WHITE);
//    addRemoveOptions(_options,_pt1,"pt1","1");
//  }
//  
//  //TODO Fix
//  public void addPoints(float[] x1, float[] x2) {
//    _pt2 = _pp.addPoints(x1,x2);
//    _pt2.setStyle("rO");
//    addRemoveOptions(_options,_pt2,"pt2","2");
//  }
//  
//  //TODO Fix
//  public void addPoints2(float[][] x1, float[][] x2) {
//    _pt2 = _pp.addPoints(x1,x2);
//    _pt2.setStyle("rO");
//    addRemoveOptions(_options,_pt2,"pt2","2");
//  }
//  
//  //TODO Fix
//  public void addPoints3(float[][][] x1, float[][][] x2) {
//    _x13 = x1;
//    _x23 = x2;
//    int[] slices = _pp.getSlices();
//    int k3 = slices[2];
//    _pt2 = _pp.addPoints(1,0,x1[k3],x2[k3]);
//    _pt2.setStyle("rO");
//    addRemoveOptions(_options,_pt2,"pt2","2");
//  }
//  
//  //TODO Fix
//  public void addPoints(float[][] x1, float[][] x2) {
//    Check.argument(x1.length==_f.length,"x1.length==f.length");
//    Check.argument(x2.length==_f.length,"x2.length==f.length");
//    _x1 = x1;
//    _x2 = x2;
//    _pt2 = _pp.addPoints(x1[_i3],x2[_i3]);
//    _pt2.setStyle("rO");
//    addRemoveOptions(_options,_pt2,"pt2","2");
//  }

  public void setTitle(String title) {
    _title = title;
    if (_i3!=Integer.MIN_VALUE)
      _pp.setTitle(_title+" "+_i3);
    else
      _pp.setTitle(_title);
  }

  public void setLabel1(String label) {
    _pp.setLabel1(label);
  }
  
  public void setLabel2(String label) {
    _pp.setLabel2(label);
  }
  
  public void setLabel3(String label) {
    _pp.setLabel3(label);
  }
  
  public void setLimits1(double min, double max) {
    if (_orientation==Orientation.X1DOWN_X2RIGHT) {
      _pp.setVLimits(1,min,max);
    } else if (_orientation==Orientation.X1DOWN_X3RIGHT) {
      _pp.setVLimits(1,min,max);
    } else if (_orientation==Orientation.X1RIGHT_X2UP) {
      _pp.setHLimits(0,min,max);
    } else if (_orientation==Orientation.X1RIGHT_X3UP) {
      _pp.setHLimits(0,min,max);
    }
  }
  
  public void setLimits2(double min2, double max2) {
    if (_orientation==Orientation.X1DOWN_X2RIGHT) {
      _pp.setHLimits(0,min2,max2);
    } else if (_orientation==Orientation.X1DOWN_X3RIGHT) {
      _pp.setHLimits(1,min2,max2);
      _pp.setVLimits(0,min2,max2);
    } else if (_orientation==Orientation.X1RIGHT_X2UP) {
      _pp.setVLimits(1,min2,max2);
    } else if (_orientation==Orientation.X1RIGHT_X3UP) {
      _pp.setHLimits(1,min2,max2);
      _pp.setVLimits(0,min2,max2);
    }
  }
  
  public void setLimits3(double min3, double max3) {
    if (_orientation==Orientation.X1DOWN_X2RIGHT) {
      _pp.setHLimits(1,min3,max3);
      _pp.setVLimits(0,min3,max3);
    } else if (_orientation==Orientation.X1DOWN_X3RIGHT) {
      _pp.setHLimits(0,min3,max3);
    } else if (_orientation==Orientation.X1RIGHT_X2UP) {
      _pp.setHLimits(1,min3,max3);
      _pp.setVLimits(0,min3,max3);
    } else if (_orientation==Orientation.X1RIGHT_X3UP) {
      _pp.setVLimits(1,min3,max3);
    }
  }
  
  public void setHInterval(double interval) {
    _pp.setHInterval(interval);
  }

  public void setVInterval(double interval) {
    _pp.setVInterval(interval);
  }
  
  public void setClips1(float clipMin, float clipMax) {
    _pp.setClips(clipMin,clipMax);
  }
  
  public void setClips2(float clipMin, float clipMax) {
    if (_pv2==null)
      throw new IllegalStateException(
          "Second PixelsView has not been added to plot panel.");
    _pv2.setClips(clipMin,clipMax);
  }

  public void setColorModel1(IndexColorModel colorModel) {
    _pp.setColorModel(colorModel);
  }
  
  public void setColorModel2(IndexColorModel colorModel) {
    _pv2.setColorModel(colorModel);
  }
  
  public void addColorBar(String label) {
    _pp.addColorBar(label);
//    _pp.setColorBarWidthMinimum(100);
  }
  
  public void setSize(int width, int height) {
    _vf.setSize(width,height);
  }
  
  public void setFontSizeForPrint(double fontSize, double plotWidth) {
    _vf.setFontSizeForPrint(fontSize,plotWidth);
  }
  
  public void setFontSizeForSlide(double fracWidth, double fracHeight) {
    _vf.setFontSizeForSlide(fracWidth, fracHeight);
  }
  
  public void paintToPng(double dpi, double win, String fileName) {
    _vf.paintToPng(dpi,win,fileName);
  }
  
  public void show() {
    _vf.setVisible(true);
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
    Viewer3 v = new Viewer3(f);
    v.show();
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private
  
  private ViewerFrame _vf;
  private PlotPanelPixels3 _pp;
  private PixelsView[] _pv1;
  private PixelsView _pv2;
  private PointsView _pt1;
  private PointsView _pt2;
  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
  private Orientation _orientation;
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
  
//  private void updateOptions(PixelsView[] pv, String label) {
//    addClipOptions(_options,pv,label);
//    addColorOptions(_options,pv,label);
//    addAlphaOptions(_options,pv,label);
//  }
  
  private class SliceMode extends Mode {

    protected SliceMode(ModeManager manager) {
      super(manager);
      setName("Change slices");
      setMnemonicKey(KeyEvent.VK_S);
      setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_S,0));
      setShortDescription("Click/drag in any tile to change slices");
    }

    @Override
    protected void setActive(Component component, boolean active) {
      if (component instanceof Tile) {
        if (active) {
          component.addMouseListener(_ml);
          component.addMouseMotionListener(_mml);
        } else {
          component.removeMouseListener(_ml);
          component.removeMouseMotionListener(_mml);
        }
      }
    }
    
    private static final long serialVersionUID = 1L;
    
    // Handles mouse pressed and released events.
    private MouseListener _ml = new MouseAdapter() {
      public void mousePressed(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        Object source = e.getSource();
        if (source instanceof Tile) {
          Tile tile = (Tile)source;
          setSlices(tile,x,y);
        }
      }
      public void mouseReleased(MouseEvent e) {} // Do nothing.
    };

    // Handles mouse dragged events.
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        int x = e.getX();
        int y = e.getY();
        Object source = e.getSource();
        if (source instanceof Tile) {
          Tile tile = (Tile)source;
          setSlices(tile,x,y);
        }
      }
    };

    private void setSlices(Tile tile, int x, int y) {
      int row = tile.getRowIndex();
      int col = tile.getColumnIndex();
//      System.out.format("row=%d, col=%d\n",row,col);
      double worldX = tile.pixelToWorldHorizontal(x);
      double worldY = tile.pixelToWorldVertical(y);
      int[] slc = _pp.getSlices();
      int k1=slc[0],k2=slc[1],k3=slc[2];
      if (_orientation==Orientation.X1DOWN_X2RIGHT) {
        if (row==0) {
          k2 = _s2.indexOfNearest(worldX);
          k3 = _s3.indexOfNearest(worldY);
        } else if (col==0) {
          k1 = _s1.indexOfNearest(worldY);
          k2 = _s2.indexOfNearest(worldX);
        } else {
          k1 = _s1.indexOfNearest(worldY);
          k3 = _s3.indexOfNearest(worldX);
        }
      }
//      System.out.format(
//          "x=%d, y=%d, worldX=%f, worldY=%f, k1=%d, k2=%d, k3=%d\n",
//          x,y,worldX,worldY,k1,k2,k3);
      _pp.setSlices(k1,k2,k3);
    }
    
  }
  
}
