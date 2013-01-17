#############################################################################
# Print figures for fault displacements in 2D images

import sys
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from fault import FaultScanner2,FaultSemblance

#############################################################################

pngDir = "./png"
#pngDir = None

sigmaTheta = 20
smoother = FaultScanner2.Smoother.SHEAR
#smoother = FaultScanner2.Smoother.FFT

def main(args):
  #goSlopes()
  #goAlign()
  #goSemblance()
  #goScan()
  #goThin()
  goShifts()

def getImage():
  #return imageSyn()
  return imageF3d()
  #return imageTpd()

def goShifts():
  s1,s2,g = getImage()
  g = slog(g)
  #plot2(s1,s2,g,title="log input")
  fse = FaultSemblance()
  g = fse.taper(10,g)
  p = fse.slopes(g)
  sn,sd = fse.semblanceNumDen(p,g)
  fsc = FaultScanner2(sigmaTheta,[sn,sd],smoother)
  f,t = fsc.scan(-15,15)
  shiftMin,shiftMax = -20,20
  faults = fsc.findFaults([f,t],shiftMax-shiftMin);
  ff = faults.getLikelihoods()
  plot2(s1,s2,g,ff,gmin=0,gmax=1,gmap=jetr,label="Fault likelihood",png="flg")
  g = fsc.smooth(4,p,ff,g)
  plot2(s1,s2,g,ff,gmin=0,gmax=1,gmap=jetr,label="Fault likelihood",png="flgs")
  plot2(s1,s2,g,label="Log amplitude",png="gs")
  p = fse.slopes(g)
  faults.findShifts(g,p,shiftMin,shiftMax)
  faults.clean()
  s = faults.getShifts()
  s = mul(s1.delta*1000.0,s)
  s = neg(s)
  print "s min =",min(s)," max =",max(s)
  plot2(s1,s2,g,s,gmin=0,gmax=28,gmap=jetr,
        label="Vertical component of throw (ms)",png="fs")
  plot2(s1,s2,g,s,gmin=0,gmax=15,gmap=jetr,
        label="Vertical component of throw (ms)",png="fs15")

def goThin():
  s1,s2,g = getImage()
  g = slog(g)
  #plot2(s1,s2,g,title="log input")
  fse = FaultSemblance()
  g = fse.taper(10,g)
  for iter in range(1):
    p = fse.slopes(g)
    #p = zerofloat(len(p[0]),len(p))
    sn,sd = fse.semblanceNumDen(p,g)
    fsc = FaultScanner2(sigmaTheta,[sn,sd],smoother)
    f,t = fsc.scan(-15,15)
    fs = copy(f); RecursiveGaussianFilter(1.0).apply00(fs,fs)
    plot2(s1,s2,g,fs,gmin=0.5,gmax=1,gmap=jetr,
          label="Fault likelihood",png="fls")
    ft,tt = fsc.thin([f,t])
    g = fsc.smooth(8,p,ft,g)
    plot2(s1,s2,g,ft,gmin=0.5,gmax=1,gmap=jetr,
          label="Fault likelihood",png="flt")
    plot2(s1,s2,g,label="Log amplitude",png="gs")

def goScan():
  s1,s2,g = getImage()
  g = slog(g)
  fse = FaultSemblance()
  g = fse.taper(10,g)
  p = fse.slopes(g)
  sn,sd = fse.semblanceNumDen(p,g)
  fsc = FaultScanner2(sigmaTheta,[sn,sd])
  st = Sampling(31,1.0,-15.0)
  for theta in st.values:
    f = fsc.likelihood(theta)
    png = "fl"+str(int(theta))
    plot2(s1,s2,g,f,gmin=0.5,gmax=1,gmap=jetr,label="Fault likelihood",png=png)
  tmin,tmax = st.first,st.last
  f,t = fsc.scan(tmin,tmax)
  plot2(s1,s2,g,f,gmin=0.5,gmax=1,gmap=jetr,label="Fault likelihood",png="fl")

def goSemblance():
  s1,s2,g = getImage()
  g = slog(g)
  fse = FaultSemblance()
  g = fse.taper(10,g)
  p = fse.slopes(g)
  sn0,sd0 = fse.semblanceNumDen(p,g)
  print "semblances for different vertical smoothings:"
  for sigma in [0,5,10,20]:
    ref = RecursiveExponentialFilter(sigma)
    sn = copy(sn0)
    sd = copy(sd0)
    ref.apply1(sn,sn)
    ref.apply1(sd,sd)
    s = fse.semblanceFromNumDen(sn,sd)
    print "sigma =",sigma," s min =",min(s)," max =",max(s)
    title = "semblance: sigma = "+str(sigma)
    png="s"+str(sigma)
    #plot2(s1,s2,g,s,gmin=0,gmax=1,title=title)
    plot2(s1,s2,g,s,gmin=0.5,gmax=1,label="Semblance",png=png)

def goAlign():
  s1,s2,g = getImage()
  g = slog(g)
  n1,n2 = len(g[0]),len(g)
  fse = FaultSemblance()
  p = fse.slopes(g)
  ref = RecursiveExponentialFilter(5)
  sn,sd = fse.semblanceNumDen(p,g)
  ref.apply1(sn,sn)
  ref.apply1(sd,sd)
  s = fse.semblanceFromNumDen(sn,sd)
  #plot2(s1,s2,g,s,gmin=0,gmax=1,title="semblance with alignment")
  plot2(s1,s2,g,s,gmin=0,gmax=1,gmap=jetf,label="Semblance",png="s5wa")
  p = zerofloat(n1,n2) # semblance with zero slopes
  sn,sd = fse.semblanceNumDen(p,g)
  ref.apply1(sn,sn)
  ref.apply1(sd,sd)
  s = fse.semblanceFromNumDen(sn,sd)
  #plot2(s1,s2,g,s,gmin=0,gmax=1,title="semblance without alignment")
  plot2(s1,s2,g,s,gmin=0,gmax=1,gmap=jetf,label="Semblance",png="s5woa")

def goSlopes():
  s1,s2,g = getImage()
  #plot2(s1,s2,g,title="input")
  plot2(s1,s2,g,label="Amplitude",png="ga")
  g = slog(g)
  #plot2(s1,s2,g,title="log input")
  plot2(s1,s2,g,label="Log amplitude",png="g")
  fse = FaultSemblance()
  p = fse.slopes(g)
  p = mul(s1.delta/s2.delta,p)
  #plot2(s1,s2,g,p,gmin=-0.9,gmax=0.9,title="slopes")
  plot2(s1,s2,g,p,gmin=-0.07,gmax=0.07,gmap=bwrf,label="Slope (s/km)",png="p")

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))
 
#############################################################################
# data read/write

def readImage(n1,n2,fileName):
  f = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(f)
  ais.close()
  return f

def writeImage(f,fileName):
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(f)
  aos.close()

def imageTpd():
  #ft,fx = 0.500,0.000
  #dt,dx = 0.004,0.025
  f1,f2 = 0.0,0.0
  d1,d2 = 1.0,1.0
  n1,n2 = 251,357
  fileName = "/data/seis/tp/csm/oldslices/tp73.dat"
  x = readImage(n1,n2,fileName)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return s1,s2,x

def imageF3d():
  n1,n2 = 462,951
  d1,d2 = 0.004,0.025
  f1,f2 = 0.004,0.000
  fileName = "/data/seis/f3d/f3d75.dat"
  x = readImage(n1,n2,fileName)
  subset = True
  if subset:
    j1,j2 = 240,0
    n1,n2 = n1-j1,440
    f1,f2 = f1+j1*d1,f2+j2*d2
    x = copy(n1,n2,j1,j2,x)
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  return s1,s2,x

def imageSyn():
  s1,s2,f = imageTpd()
  n1,n2 = s1.count,s2.count
  g = zerofloat(n1,n2)
  f1 = copy(f[n2/8])
  f2 = shiftRamp(f1)
  for i2 in range(n2):
    for i1 in range(n1):
      if i2<n2/4:
        copy(f1,g[i2])
      elif i2<n2/3:
        copy(f2,g[i2])
      elif i2<2*n2/3:
        copy(f1,g[i2])
      else:
        copy(f2,g[i2])
  #zero(g[1*n2/3])
  #zero(g[2*n2/3])
  #r = randomNoise(1.0,n1,n2)
  #g = add(g,r)
  #rgf = RecursiveGaussianFilter(2.0)
  #rgf.applyX0(g,g)
  return s1,s2,g

def shiftRamp(f):
  n = len(f)
  g = copy(f)
  t = rampfloat(0.0,1.0-8.0/(n-1),n)
  si = SincInterp()
  si.interpolate(n,1.0,0.0,f,n,t,g)
  return g

def randomNoise(a,n1,n2):
  ran = Random(3)
  r = mul(2.0*a,sub(randfloat(ran,n1,n2),0.5))
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply00(r,r)
  return r
 
#############################################################################
# plotting

def plot2(s1,s2,f,g=None,gmin=None,gmax=None,gmap=None,
          label=None,title=None,png=None):
  n1 = len(f[0])
  n2 = len(f)
  panel = panel2()
  if label:
    panel.addColorBar(label)
  else:
    panel.addColorBar()
  panel.setColorBarWidthMinimum(175)
  panel.setHLabel("Inline (km)")
  panel.setVLabel("Time (s)")
  pv = panel.addPixels(s1,s2,f)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.GRAY)
  if g:
    pv = panel.addPixels(s1,s2,g)
    pv.setInterpolation(PixelsView.Interpolation.NEAREST)
    if gmin==None: gmin = min(g)
    if gmax==None: gmax = max(g)
    if gmap==None: gmap = jetf
    pv.setClips(gmin,gmax)
    pv.setColorModel(gmap)
  frame2(panel,png)

def jetFill(alpha):
  return ColorMap.setAlpha(ColorMap.JET,alpha)
def bwrFill(alpha):
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,alpha)
def jetRamp(alpha):
  return ColorMap.setAlpha(ColorMap.JET,rampfloat(0.0,alpha/256,256))
def bwrNotch(alpha):
  a = zerofloat(256)
  for i in range(len(a)):
    if i<128:
      a[i] = alpha*(128.0-i)/128.0
    else:
      a[i] = alpha*(i-127.0)/128.0
  return ColorMap.setAlpha(ColorMap.BLUE_WHITE_RED,a)
jetf = jetFill(0.5)
jetr = jetRamp(1.0)
bwrf = bwrFill(0.5)
bwrn = bwrNotch(1.0)


def panel2():
  #panel = PlotPanel(1,1,
  #  PlotPanel.Orientation.X1DOWN_X2RIGHT,
  #  PlotPanel.AxesPlacement.NONE)
  panel = PlotPanel(1,1,
    PlotPanel.Orientation.X1DOWN_X2RIGHT)
  return panel

def frame2(panel,png=None):
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  #frame.setBackground(Color(0xfd,0xfe,0xff)) # easy to make transparent
  frame.setFontSizeForPrint(8,240) # single column
  frame.setSize(1240,770)
  frame.setVisible(True)
  if png and pngDir:
    frame.paintToPng(720,3.3,pngDir+"/"+png+".png")
  return frame

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

"""
Outline:
find slopes p
use p to align minus-plus images fm and fp
for all thetas
  shear fm and fp
  compute correlation coefficients c
  unshear c
  compute fault likelihood d
  remember dmax and corresponding theta

smooth dmax laterally and pick peaks
for all fault curves
  gather samples on both sides of fault
  cross-correlate to find displacements
"""
