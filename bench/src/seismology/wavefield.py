#############################################################################
# Demo of 2 wavelet estimations from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.awt.ColorMap import *
from edu.mines.jtk.lapack import *
from seismology import Wavefield
from java.util import Random

############################################################################
pngDir = "./png/figures/"
#pngDir = None

def main(args):
  goSimpleTest()

def goSimpleTest():
  #parameters
  dens = 2500 #(kg/m^3)
  pwv = 2000.0 #m/s
  swv = 1000.0 #m/s
  freq = 25.0 #Hz
  period = 1.0/freq

  theta = 00.0 #(degrees)
  rnf = 80.0 #(meters)
  rff = 1000.0 #(meters)
  wf = Wavefield()
  wf.setMediumProp(dens,pwv,swv)
  wf.setFreq(freq)

  ft = 0.00
  dt = .001
  nx0 = (int)(period/dt+1.0)
  fcomp = 1
  nthetas = 10
  dthetas = 10
  ucomp = 1
  nt = (int) (rff/swv/dt+1.25*period/dt)
  st = Sampling(nt,dt,ft)

  ucomp=1
  r = rnf
  u1rnfthetas = wf.calcUThetas(ucomp,fcomp,nthetas,nt,ft,dt,r)
  ucomp=1
  r = rff
  u1rffthetas = wf.calcUThetas(ucomp,fcomp,nthetas,nt,ft,dt,r)
  ucomp=3
  r = rnf
  u3rnfthetas = wf.calcUThetas(ucomp,fcomp,nthetas,nt,ft,dt,r)
  ucomp=3
  r = rff
  u3rffthetas = wf.calcUThetas(ucomp,fcomp,nthetas,nt,ft,dt,r)
  u = [u1rnfthetas,u1rffthetas,u3rnfthetas,u3rffthetas]
  umax = max(abs(u))


  plotHolo = False 
  plotSeis = False 
  if plotHolo:
    r = rnf
    plotHodograms(nthetas,u1rnfthetas,u3rnfthetas,title="r = "+str(r)+" (m)",pngDir=pngDir)
    r = rff
    plotHodograms(nthetas,u1rffthetas,u3rffthetas,title="r = "+str(r)+" (m)",pngDir=pngDir)


  if plotSeis:
    ucomp=1
    r = rnf
    umax = max(abs(u1rnfthetas))
    for i in range(0,nthetas):
      label = "u"+str(ucomp)+" (m)" 
      title = str(i*dthetas)+" degrees from the horizontal (x1) and r = "+str(r)+" (m)"
      plotSequences(st,[u1rnfthetas[i]],amax=umax,labels=[label],title=title,pngDir=pngDir)

    ucomp=1
    r = rff
    umax = max(abs(u1rffthetas))
    for i in range(0,nthetas):
      label = "u"+str(ucomp)+" (m)" 
      title = str(i*dthetas)+" degrees from the horizontal (x1) and r = "+str(r)+" (m)"
      plotSequences(st,[u1rffthetas[i]],amax=umax,labels=[label],title=title,pngDir=pngDir)

    ucomp=3
    r = rnf
    umax = max(abs(u3rnfthetas))
    for i in range(0,nthetas):
      label = "u"+str(ucomp)+" (m)" 
      title = str(i*dthetas)+" degrees from the horizontal (x1) and r = "+str(r)+" (m)"
      plotSequences(st,[u3rnfthetas[i]],amax=umax,labels=[label],title=title,pngDir=pngDir)

    ucomp=3
    r = rff
    umax = max(abs(u3rffthetas))
    for i in range(0,nthetas):
      label = "u"+str(ucomp)+" (m)" 
      title = str(i*dthetas)+" degrees from the horizontal (x1) and r = "+str(r)+" (m)"
      plotSequences(st,[u3rffthetas[i]],amax=umax,labels=[label],title=title,pngDir=pngDir)

  x0 = zerofloat(nx0)
  for i in range(0,nx0):
    t = i*dt
    x0[i] = sin(2.0*PI*freq*t)
  plotSequences(Sampling(nx0,dt,ft),[x0],amax=1.0,labels=["x0"],title="Source Wavelet",pngDir=pngDir)


def plotSequences(st,xs,amax=None,labels=None,title=None,pngDir=None):
  nx = len(xs)
  pp = PlotPanel(nx,1)
  dt = st.getDelta()
  for ix,xi in enumerate(xs):
    pv = pp.addPoints(ix,0,st,xi)
    if labels:
      pp.setVLabel(ix,labels[ix])
    if amax:
      pp.setVLimits(ix,-amax,amax)
      pp.setVInterval(ix,amax/5)
  pp.setHLabel("Time (s)")
  pp.setHInterval(ix,100*dt)
  pf = PlotFrame(pp)
  pf.setVisible(True)
  if title:
    pp.setTitle(title)
  if pngDir:
    pf.paintToPng(360,3.08,pngDir+labels[0]+" "+title+".png")

def plotHodogram(xs1,xs2,min1=None,min2=None,max1=None,max2=None,theta=None,title=None,pngDir=None):
  nx = len(xs1)
  pp = PlotPanel(1,1)
  pv = pp.addPoints(0,0,xs1,xs2)
  pp.setVLabel(0,"U3(t) (m)")
  if min1:
    pp.setVLimits(0,min2,max2)
    pp.setHLimits(0,min1,max1)
  pp.setHLabel("U1(t) (m)")
  pf = PlotFrame(pp)
  pf.setVisible(True)
  if title:
    pp.setTitle(title)
  if pngDir:
    pf.paintToPng(360,3.08,pngDir+title+".png")

def plotHodograms(nthetas,u1s,u2s,amax=None,labels=None,title=None,pngDir=None):
  umax1 = max(u1s)
  umax2 = max(u2s)
  umin1 = min(u1s)
  umin2 = min(u2s)
  umax = 0.0
  umin = 0.0
  if umax1>umax2:
    umax = umax1
  else:
    umax = umax2
  if umin1>umin2:
    umin = umin2
  else:
    umin = umin1


  for i in range(0,nthetas):
    titler = str(10*i) + title
    plotHodogram(u1s[i],u2s[i],umin,umin,umax,umax,(i*10),titler,pngDir)



#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

