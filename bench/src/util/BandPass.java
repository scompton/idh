package util;
import edu.mines.jtk.dsp.BandPassFilter;
import edu.mines.jtk.dsp.BandPassFilter.Extrapolation;
import edu.mines.jtk.util.Parallel;
import edu.mines.jtk.util.Stopwatch;

public class BandPass {

  public static float[][][] goBandPass(
      final double klower, final double kupper, 
      final double kwidth, final double aerror,
      final float[][][] f) 
  {
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n1 = f[0][0].length;
    final float[][][] ff = new float[n3][n2][n1];
    final BandPassFilter bpf = new BandPassFilter(klower,kupper,kwidth,aerror);
    bpf.setExtrapolation(Extrapolation.ZERO_SLOPE);
    bpf.setFilterCaching(true);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          bpf.apply(f[i3][i2],ff[i3][i2]);
        }
      }
    }});
    return ff;
  }
  
  public static float[][] goBandPass(
      final double klower, final double kupper, 
      final double kwidth, final double aerror,
      final float[][] f) 
  {
    Stopwatch s = new Stopwatch();
    s.start();
    final int n2 = f.length;
    final int n1 = f[0].length;
    final float[][] ff = new float[n2][n1];
    final BandPassFilter bpf = new BandPassFilter(klower,kupper,kwidth,aerror);
    bpf.setExtrapolation(Extrapolation.ZERO_SLOPE);
    bpf.setFilterCaching(true);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; i1++)
        bpf.apply(f[i2],ff[i2]);
    }});
    System.out.println("Finished in "+s.time()+" seconds");
    return ff;
  }
  
}
