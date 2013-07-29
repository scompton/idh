package warp;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.util.ArrayMath;

public class ErrorMapUtil {

  /**
   * 
   * @param sl
   * @param s1
   * @param s2
   */
  public ErrorMapUtil(Sampling sl, Sampling s1, Sampling s2) {
    this(sl,s1,s2,null);
  }

  /**
   * 
   * @param sl
   * @param s1
   * @param s2
   * @param s3
   */
  public ErrorMapUtil(Sampling sl, Sampling s1, Sampling s2, Sampling s3) {
    _errorMap = new HashMap<>();
    _sl = sl;
    _s1 = s1;
    _s2 = s2;
    _s3 = s3;
  }

  /**
   * 
   * @param xl
   * @param x1
   * @param x2
   */
  public void add(double xl, double x1, double x2) {
    int il = _sl.indexOfNearest(xl);
    int i1 = _s1.indexOfNearest(x1);
    int i2 = _s2.indexOfNearest(x2);
    List<Integer> list = _errorMap.get(i2);
    if (list==null)
      list = new LinkedList<>();
    list.add(il);
    list.add(i1);
    _errorMap.put(i2,list);
  }

  /**
   * 
   * @param xl
   * @param x1
   * @param x2
   * @param x3
   */
  public void add(double xl, double x1, double x2, double x3) {
    int il = _sl.indexOfNearest(xl);
    int i1 = _s1.indexOfNearest(x1);
    int i2 = _s2.indexOfNearest(x2);
    int i3 = _s3.indexOfNearest(x3);
    int i23 = i2*i3;
    List<Integer> list = _errorMap.get(i23);
    if (list==null)
      list = new LinkedList<>();
    list.add(il);
    list.add(i1);
    _errorMap.put(i23,list);
  }

  /**
   * 
   * @return
   */
  public Map<Integer,int[]> getMap() {
    Map<Integer,int[]> map = new HashMap<>();
    for (Integer integer : _errorMap.keySet()) {
      List<Integer> list = _errorMap.get(integer);
      int n = list.size();
      int[] l1 = new int[n];
      for (int i=0; i<n; i++)
        l1[i] = list.get(i);
      map.put(integer,l1);
    }
    return map;
  }

  public static void main(String[] args) {
    int nl = 169;
    int n1 = 1501;
    int n2 = 1211;
    int n3 = 826;
    double dt = 2.0;
    double dx = 22.5;
    double dy = 12.5;
    Sampling sl = new Sampling(nl,dt,-200.0);
    Sampling s1 = new Sampling(n1,dt,0.0);
    Sampling s2 = new Sampling(n2,dx,1001);
    Sampling s3 = new Sampling(n3,dy,566);

    // 2D
    System.out.println("2D");
    ErrorMapUtil emu = new ErrorMapUtil(sl,s1,s2);
    emu.add(-100.0,5.0,1200.0);
    emu.add(-75.0,12.0,1200.0);
    emu.add(-150.0,300.0,1200.0);
    emu.add(0.0,0.0,1800.0);
    Map<Integer,int[]> map = emu.getMap();
    for (Integer integer : map.keySet()) {
      int[] l1 = map.get(integer);
      System.out.println("il coordinate array for index "+integer+":");
      ArrayMath.dump(l1);
    }
    
    // 3D
    System.out.println("3D");
    emu = new ErrorMapUtil(sl,s1,s2,s3);
    emu.add(-100.0,5.0,1200.0,566.0);
    emu.add(-75.0,12.0,1200.0,650.0);
    emu.add(-150.0,300.0,1200.0,650.0);
    emu.add(0.0,0.0,1800.0,800.0);
    map = emu.getMap();
    for (Integer integer : map.keySet()) {
      int[] l1 = map.get(integer);
      System.out.println("il coordinate array for index "+integer+":");
      ArrayMath.dump(l1);
    }
  }

  private Map<Integer,List<Integer>> _errorMap;
  private Sampling _sl, _s1, _s2, _s3;
}
