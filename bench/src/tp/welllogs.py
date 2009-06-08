import sys
from math import *

from java.awt import *
from java.io import *
from java.lang import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.sgl.test import *

from tp import *

# Directories and files for well logs, headers, directional surveys
tpDir = "/data/seis/tp/"
doeWellLogsDir = tpDir+"doe/WellLogs/"
csmWellLogsDir = tpDir+"csm/welllogs/"
doeWellHeaders = doeWellLogsDir+"WellHeaders.txt"
doeDirectionalSurveys = doeWellLogsDir+"DirectionalSurveys.txt"
doeWellLogs = ""
csmWellLogs = ""

# File for 3D seismic depth image to view with wells.
csmSeismiczDir = tpDir+"csm/seismicz/"
csmSeismic = csmSeismiczDir+"tpsz.dat"
s1c = Sampling(2762,0.002,0.000)
s2c = Sampling(161,0.025,0.000)
s3c = Sampling(357,0.025,0.000)

def main(args):
  #makeBinaryWellLogs("test")
  #makeBinaryWellLogs("shallow")
  #makeBinaryWellLogs("deep")
  #makeBinaryWellLogs("all")
  #viewWellCoordinates("deep")
  #viewWellsWithSeismic("deep","velocity")
  viewWellsWithSeismic("deep","gamma")

def setGlobals(what):
  global csmWellLogs,doeWellLogs
  if what=="all":
    doeWellLogs = doeWellLogsDir+"LAS_log_files/"
    csmWellLogs = csmWellLogsDir+"tpwall.dat"
  elif what=="deep":
    doeWellLogs = doeWellLogsDir+"LAS_log_files/Deeper_LAS_files/"
    csmWellLogs = csmWellLogsDir+"tpwdeep.dat"
  elif what=="shallow":
    doeWellLogs = doeWellLogsDir+"LAS_log_files/Shallow_LAS_files/"
    csmWellLogs = csmWellLogsDir+"tpwshallow.dat"
  elif what=="test":
    doeWellLogs = doeWellLogsDir+"LAS_log_files/test/"
    csmWellLogs = csmWellLogsDir+"tpwtest.dat"

def makeBinaryWellLogs(what):
  wldata = WellLog.Data(doeWellLogs,doeWellHeaders,doeDirectionalSurveys)
  wldata.printInfo()
  wldata.writeBinary(csmWellLogs)

def dumpWellHeaders(what):
  setGlobals(what)
  whdata = WellHeader.Data(doeWellHeaders)
  whdata.printInfo()

def dumpDirectionalSurveys(what):
  setGlobals(what)
  dsdata = DirectionalSurvey.Data(doeDirectionalSurveys)
  dsdata.printInfo()

def viewWellCoordinates(what):
  setGlobals(what)
  wldata = WellLog.Data.readBinary(csmWellLogs)
  wldata.printInfo()
  spz1 = SimplePlot()
  sp23 = SimplePlot()
  zlist = []
  x1list = []
  x2list = []
  x3list = []
  for log in wldata.getAll():
    zlist.append(log.z)
    x1list.append(log.dx1dz())
    x2list.append(log.x2)
    x3list.append(log.x3)
  tv = PointsView(zlist,x1list)
  spz1.add(tv)
  tv = PointsView(x2list,x3list)
  sp23.add(tv)

def viewWellsWithSeismic(what,curve):
  setGlobals(what)
  wdata = WellLog.Data.readBinary(csmWellLogs)
  n1c,n2c,n3c = s1c.count,s2c.count,s3c.count
  ais = ArrayInputStream(csmSeismic)
  x = Array.zerofloat(n1c,n2c,n3c)
  ais.readFloats(x)
  ais.close()
  print "x min =",Array.min(x)," max =",Array.max(x)
  ipg = ImagePanelGroup(s1c,s2c,s3c,x)
  world = World()
  world.addChild(ipg)
  addWellGroups(world,wdata,curve)
  frame = TestFrame(world)
  frame.setVisible(True)

def addWellGroups(world,wdata,curve):
  for log in wdata.getLogsWith(curve):
    pg = makePointGroup(log)
    world.addChild(pg)

def makePointGroup(log):
  n = log.n
  xyz = Array.zerofloat(3*n)
  Array.copy(n,0,1,log.x3,0,3,xyz)
  Array.copy(n,0,1,log.x2,1,3,xyz)
  Array.copy(n,0,1,log.x1,2,3,xyz)
  states = StateSet()
  cs = ColorState()
  cs.setColor(Color.YELLOW)
  states.add(cs)
  pg = PointGroup(xyz)
  pg.setStates(states)
  return pg

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
