package util;

import edu.mines.jtk.util.ArrayMath;
import edu.mines.jtk.util.Check;

public class LagIndexTest {

  public LagIndexTest(int shiftMin, int shiftMax) {
    Check.argument(shiftMax-shiftMin>1, "shiftMax-shiftMin>1");
    _lmin = shiftMin;
    _nl = 1+shiftMax-shiftMin;
  }
  
  public int getLagIndex() {
    int il = -1;
    int nlm1 = _nl-1;
    il = (_lmin+nlm1<=0)?nlm1:(_lmin<=0)?ArrayMath.abs(_lmin):0;
    return il;
  }
  
  public static void main(String[] args) {
    LagIndexTest lit = new LagIndexTest(Integer.parseInt(args[0]),
        Integer.parseInt(args[1]));
    int il = lit.getLagIndex();
    System.out.println("index="+il);
  }

  private final int _lmin, _nl;
}
