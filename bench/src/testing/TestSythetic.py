from imports import *
from synthetic import Synthetic

from data import SUDataGrabber
from edu.mines.jtk.awt import ColorMap;

from random import *


def main(args):
  st,sx,f = getVikingGrabenCDP(500)
  plotGather(st,sx,f)
  """
  nref = 100
  vel = 2
  nx,dx,fx = 100,.005,0
  nt,dt,ft = 100,.004,0

  sx = Sampling(nx,dx,fx)
  st = Sampling(nt,dt,ft)
  p = Synthetic.makeCMPReflections(vel,nref,st,sx)
  SimplePlot.asPixels(p)
  """

def getVikingGrabenCDP(icdp):
  if 145 <= icdp | icdp <= 147:
    fileName = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp145_147.dat"
    print "a"
    return getCDP(icdp, 145, fileName)
  elif 400 <= icdp | icdp <= 800:
    fileName = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp400_800.dat"
    print "b"
    return getCDP(icdp, 400, fileName)
  elif 1100 <= icdp | icdp <= 1500:
    fileName = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp1100_1500.dat"
    print "c"
    return getCDP(icdp, 1100, fileName)
  elif 1700 <= icdp | icdp <= 2000:
    fileName = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp1700_2000.dat"
    print "d"
    return getCDP(icdp, 1700, fileName)
  else:
    print "Available CDPs 145-147, 400-800, 1100-1500, and 1700-2000"
    print " "
    print " "
    return None

def getCDP(icdp, startCDP, fileName):
  st = Sampling(1500,0.004,0.004)
  if icdp%2==0:
    sx = Sampling(60,0.050,0.262) # sampling for even cdps
  else:
    sx = Sampling(60,0.050,0.287) # sampling for odd cdps
  nt,dt,ft = st.count,st.delta,st.first
  nx,dx,fx = sx.count,sx.delta,sx.first
  f = zerofloat(nt,nx)
  af = ArrayFile(fileName, "r", ByteOrder.LITTLE_ENDIAN, ByteOrder.LITTLE_ENDIAN)
  af.skipBytes(4*nt*nx*(icdp-startCDP))
  af.readFloats(f)
  af.close()
  for ix in range(nx/2):
    fi = f[ix]
    f[ix] = f[nx-1-ix]
    f[nx-1-ix] = fi
  return st,sx,f

def plotGather(st,sx,p,tmin=None,tmax=None,perc=None,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.setHLabel("Offset (km)")
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  if tmin:
    sp.setVLimits(tmin,tmax)
  pv = sp.addPixels(st,sx,p)
  if perc:
    pv.setPercentiles(100-perc,perc)
  sp.getPlotPanel()


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())



