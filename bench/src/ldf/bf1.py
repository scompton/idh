# Test 1D bilateral filter.
import sys
from math import *
from java.awt import *
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

from ldf import *

pngDir = None
#pngDir = "png/blf/"

#############################################################################
# functions

def main(args):
  goRandomBlocks()
  #goShowRangeFunctions()

gauss = BilateralFilter.Type.GAUSS
huber = BilateralFilter.Type.HUBER
tukey = BilateralFilter.Type.TUKEY

def goRandomBlocks():
  n1 = 801
  x = makeBlocks(n1)
  y = add(x,makeRandom(3141,8.0,3.0,n1))
  yrms,ymad,yqqd = rms(y),mad(y),qqd(y)
  print "y rms =",yrms," mad =",ymad," qqd =",yqqd
  plot2(x,y,-12,12,"rbxy")
  sigmaS = 20.0
  z = zerofloat(n1)
  for scale in [0.01,0.5,1.0,1.5,10,100.0]:
    sigmaX = scale*yqqd
    print "sigmaS =",sigmaS," sigmaX =",sigmaX
    bf = BilateralFilter(sigmaS,sigmaX)
    bf.setType(tukey)
    bf.apply(y,z)
    png = "rbxz"+str(int(scale*10+0.5))
    if len(png)==1: png = "0"+png
    plot2(x,z,-12,12,png=png)

def goShowRangeFunctions():
  xmin,xmax,sigma = -3.5,3.5,1.0
  nx = 351
  dx = (xmax-xmin)/(nx-1)
  fx = xmin
  sx = Sampling(nx,dx,fx)
  yg = BilateralFilter.sampleRangeFunction(gauss,sigma,sx)
  yt = BilateralFilter.sampleRangeFunction(tukey,sigma,sx)
  #yh = BilateralFilter.sampleRangeFunction(huber,sigma,sx)
  #x = rampfloat(fx,dx,nx)
  #yg = mul(x,yg)
  #yh = mul(x,yh)
  #yt = mul(x,yt)
  sp = SimplePlot()
  sp.setFontSizeForPrint(8,240)
  sp.setHLabel("x")
  sp.setVLabel("r(x)")
  sp.setSize(700,500)
  solid = PointsView.Line.SOLID
  dash = PointsView.Line.DASH
  dot = PointsView.Line.DOT
  #pv = sp.addPoints(sx,yh); pv.setLineStyle(dot); pv.setLineWidth(4)
  pv = sp.addPoints(sx,yg); pv.setLineStyle(dash); pv.setLineWidth(4)
  pv = sp.addPoints(sx,yt); pv.setLineStyle(solid); pv.setLineWidth(4)
  if pngDir:
    sp.paintToPng(720,3,pngDir+"fx.png")

def plot(x,ymin=0.0,ymax=0.0):
  sp = SimplePlot.asPoints(x)
  if ymin<ymax:
    sp.setVLimits(ymin,ymax)

def plot2(x,y,ymin=0.0,ymax=0.0,png=None):
  sp = SimplePlot()
  sp.setFontSizeForPrint(8,240)
  sp.setHLabel("Sample index")
  sp.setVLabel("Amplitude")
  sp.setSize(720,300)
  pv = sp.addPoints(x) 
  pv.setLineWidth(3)
  pv.setLineStyle(PointsView.Line.DOT)
  pv = sp.addPoints(y) 
  pv.setLineWidth(3)
  #pv.setLineColor(Color.RED)
  if ymin<ymax:
    sp.setVLimits(ymin,ymax)
  if png and pngDir:
    sp.paintToPng(720,3,pngDir+png+".png")

def makeBlocks(n1):
  nb = 17
  db = 1.0
  pb = 8.0
  sb = 1.0
  xb = pb
  m1 = 1+n1/nb
  x = zerofloat(n1)
  for i1 in range(n1):
    if (i1+1)%m1==0:
      pb = pb-db
      sb = -sb
    x[i1] = sb*pb
  return x

def makeRandom(seed,scale,sigma,n1):
  r = Random(seed)
  x = mul(2.0*scale,sub(randfloat(r,n1),0.5))
  rgf = RecursiveGaussianFilter(sigma)
  rgf.apply0(x,x)
  return x

def makeImpulse(n1):
  x = zerofloat(n1)
  x[n1/2] = 1.0
  return x

def rms(x):
  n = len(x)
  m = sum(x)/n
  x = sub(x,m)
  return sqrt(sum(mul(x,x))/n)

def mad(x):
  m = Quantiler.estimate(0.5,x)
  y = abs(sub(x,m))
  return Quantiler.estimate(0.5,y)

def qqd(x):
  return 0.5*(Quantiler.estimate(0.75,x)-Quantiler.estimate(0.25,x))

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
