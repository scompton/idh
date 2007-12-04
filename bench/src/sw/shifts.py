import sys
from math import *
from java.awt import *
from java.lang import *
from javax.swing import *
from lcc import *
from sw import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *

True = 1
False = 0

n1 = 1501  
d1 = 0.004
f1 = 0.0
s1 = Sampling(n1,d1,f1)

n2 = 623  
d2 = 0.025
f2 = 0
s2 = Sampling(n2,d2,f2)

n3 = 367
d3 = 0.025
f3 = 0
s3 = Sampling(n3,d3,f3)

datadir = "/data/seis/sw/all/"
datbdir = "/datb/seis/sw/all/"
datcdir = "/datb/seis/sw/all/"

##############################################################################
# Read/write

def readFloats3(file):
  f = Array.zerofloat(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readFloats(f)
  af.close()
  return f

def writeFloats3(file,f):
  af = ArrayFile(datadir+file,"rw")
  af.writeFloats(f)
  af.close()
  return f

def readBytes3(file):
  b = Array.zerobyte(n1,n2,n3)
  af = ArrayFile(datadir+file,"r")
  af.readBytes(b)
  af.close()
  return b

def writeBytes3(file,lag):
  af = ArrayFile(datadir+file,"rw")
  af.writeBytes(lag)
  af.close()

def readFloats12(file,i3):
  f = readFloats3(file)
  return slice12(f,i3)

def readFloats13(file,i2):
  f = readFloats3(file)
  return slice13(f,i2)

def readFloats23(file,i1):
  f = readFloats3(file)
  return slice23(f,i1)

def slice12(f,i3):
  f12 = Array.copy(f[i3])
  return f12

def slice13(f,i2):
  f13 = Array.zerofloat(n1,n3)
  for i3 in range(n3):
    Array.copy(n1,f[i3][i2],f13[i3])
  return f13

def slice23(f,i1):
  f23 = Array.zerofloat(n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      f23[i3][i2] = f[i3][i2][i1];
  return f23

def copyFile(x,y):
  x = ArrayFile(x,"r")
  y = ArrayFile(y,"r")
  n = 65536
  b = Array.zerobyte(n)
  size = x.length()
  while size>0:
    m = min(n,size)
    x.readBytes(b,0,m)
    y.writeBytes(b,0,m)
    size -= m
  x.close()
  y.close()

#############################################################################
# Research

class WhitenFilter(FileFloat3Chunks.Filter):
  def __init__(self,sigma1,sigma2,sigma3):
    self.sf = ShiftFinder(sigma1,sigma2,sigma3)
  def apply(self,i1,i2,i3,x,y):
    x,y = x[0],y[0]
    n1,n2,n3 = len(x[0][0]),len(x[0]),len(x)
    xmin,xmax = Array.min(x),Array.max(x)
    tiny = 1.0e-5*xmax
    print "WhitenFilter.apply"
    print "  min =",xmin," max =",xmax
    print "  i1 =",i1," i2 =",i2," i3 =",i3
    print "  n1 =",n1," n2 =",n2," n3 =",n3
    r = Array.randfloat(n1,n2,n3)
    Array.mul(tiny,r,r)
    Array.add(r,x,x)
    r = None
    self.sf.whiten(x,y)

def whitenArrayFile(afx,afy):
  sigma1,sigma2,sigma3 = 12,6,6
  wf = WhitenFilter(sigma1,sigma2,sigma3)
  #mc = 50000000 # 50 Mfloats
  mc = 120000000 # 120 Mfloats
  l1,l2,l3 = 3*sigma1,3*sigma2,3*sigma3
  ff3c = FileFloat3Chunks(mc,n1,l1,l1,n2,l2,l2,n3,l3,l3)
  print "whitenArrayFile: chunk size =",ff3c.getChunkSize()
  ff3c.apply(wf,[afx],[afy])

def whiten():
  #fs = ["s02","s04"]
  #fw = ["w02","w04"]
  fs = ["s04"]
  fw = ["w04"]
  for i in range(len(fs)):
    sname = datadir+fs[i]+".dat"
    wname = datadir+fw[i]+".dat"
    print "whitening",sname,"..."
    afs = ArrayFile(sname,"r")
    afw = ArrayFile(wname,"rw")
    whitenArrayFile(afs,afw)
    afs.close()
    afw.close()
    print "... done"

class ShiftFinderFilter(FileFloat3Chunks.Filter):
  def __init__(self,niter,lmin,lmax,sigma1,sigma2,sigma3):
    self.sf = ShiftFinder(sigma1,sigma2,sigma3)
    self.niter = niter
    self.lmin = lmin
    self.lmax = lmax
  def apply(self,i1,i2,i3,x,y):
    print "ShiftFinderFilter.apply:"
    print "  i1 =",i1," i2 =",i2," i3 =",i3
    print "  n1 =",n1," n2 =",n2," n3 =",n3
    f = x[0] # input array, not modified
    g,u1,u2,u3 = y[0],y[1],y[2],y[3] # input/output arrays
    n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
    du = Array.zerofloat(n1,n2,n3)
    gs = Array.zerofloat(n1,n2,n3)
    print "  shift1"
    sf.find1(lmin,lmax,f,g,du)
    print "    du min =",Array.min(du)," max =",Array.max(du)
    sf.shift1(du,g,gs,u1,u2,u3)
    print "    shift1 complete"
    Array.copy(gs,g)
    print "  shift3"
    sf.find3(lmin,lmax,f,g,du)
    print "    du min =",Array.min(du)," max =",Array.max(du)
    sf.shift3(du,g,gs,u1,u2,u3)
    print "    shift3 complete"
    Array.copy(gs,g)
    print "  shift2"
    sf.find2(lmin,lmax,f,g,du)
    print "    du min =",Array.min(du)," max =",Array.max(du)
    sf.shift2(du,g,gs,u1,u2,u3)
    print "    shift2 complete"
    Array.copy(gs,g)

def findShifts():
  nshift = 3
  lmin,lmax = -2,2
  sigma1,sigma2,sigma3 = 12,12,12
  sff = ShiftFinderFilter(lmin,lmax,sigma1,sigma2,sigma3)
  fname = datadir+"w02.dat"
  gname = datadir+"w04.dat"
  faf = ArrayFile(fname,"r")
  for ishift in range(1,nshift+1):
    print "findShifts: iteration =",ishift
    gsname = datadir+"gs"+ishift+".dat"
    u1name = datadir+"u1s"+ishift+".dat"
    u2name = datadir+"u2s"+ishift+".dat"
    u3name = datadir+"u3s"+ishift+".dat"
    copyFile(gname,gsname)
    gname = gsname
    gaf = ArrayFile(gname,"rw")
    u1f = ArrayFile(u1name,"rw")
    u2f = ArrayFile(u2name,"rw")
    u3f = ArrayFile(u3name,"rw")
    x = [faf,gaf,u1f,u2f,u3f]
    y = [    gaf,u1f,u2f,u3f]
    sff.apply(sff,x,y)
    gaf.close()
    u1f.close()
    u2f.close()
    u3f.close()

def main(args):
  #whiten()
  findShifts()
  return

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
