#############################################################################
#builds synthetic figures for CWP presentation

from imports import *

from edu.mines.jtk.dsp.Conv import *
from wwarp import WaveletWarping
from wwarp import ShapingFilter
from java.util import Random

############################################################################

#pngDir = "./pres14/figures/"
#pngDir = "./png/figuresAbs/"
pngDir = "./pres14/pqfgdh/"
#pngDir = None

def main(args):
  goSimpleTest()

def goSimpleTest():
  nt,ni = 481,2 # number of time samples; number of impulses
  freq,decay = 0.08,0.05 # peak frequency and decay for wavelet
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  amax = 5
  st = Sampling(nt,dt,ft)
  for mp in [False]: # True, for minimum-phase; False for other
    hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
    for r in [2.0]: # for stretch and squeeze, ...
      fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
      u = rampfloat(0.0,r,nt)
      if r<=1.0:
        u = add((1.0-r)*(nt-1),u)
      ww = WaveletWarping()
      ww.setFrequencyRange(fmin,fmax)
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)

      p,q = makeImpulses(r,nt,ni)

      lq = ww.applyL(u,q)
      slq = ww.applyS(u,lq)

      f = addWavelet(freq,decay,p,mp)
      g = addWavelet(freq,decay,q,mp)

      lg = ww.applyL(u,g)
      slg = ww.applyS(u,lg)
      slgdiv2 = div(slg,2)

      hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
      ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
      aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet

      ag = ww.applyA(na,ka,ak,g)
      lag = ww.applyL(u,ag) # lowpass, if squeezing
      slag = ww.applyS(u,lag)

      aslg = ww.applyA(na,ka,ak,slg)

      hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet

      hslag = ww.applyH(nh,kh,hw,slag)

      undelay = False 
      d = ww.displayDifferencesShifts(na,ka,u,f,g,undelay)

      undelay = True 
      dud = ww.displayDifferencesShifts(na,ka,u,f,g,undelay)

      nhw = normalize(hw)
      nhk = normalize(hk)
      naw = normalize(aw)
      nak = normalize(ak)

      title="pq"
      plot2TracesSideBySide(st,p,q,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="pSLq"
      plot2TracesSideBySide(st,p,slq,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="fg"
      plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="fSg"
      plot2TracesSideBySide(st,f,slg,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="fSgdiv2"
      plot2TracesSideBySide(st,f,slgdiv2,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="fAg"
      plot2TracesSideBySide(st,f,ag,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="fSLAg"
      plot2TracesSideBySide(st,f,slag,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="SLAgASLg"
      plot2TracesSideBySide(st,f,aslg,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="fHSLAg"
      plot2TracesSideBySide(st,f,hslag,tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="3d"
      plot3TracesSideBySide(st,d[0],d[1],d[2],tmin*dt,tmax*dt,5,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="0dud"
      plot1Trace(st,d[0],tmin*dt,250*dt,6,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="1dud"
      plot1Trace(st,d[1],tmin*dt,250*dt,6,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)
      title="2dud"
      plot1Trace(st,d[2],tmin*dt,250*dt,6,0.90,0.80,16.0/9.0,
      title=title,pngDir=pngDir,twocol=True)

 
      title="hwaveletsfg"
      plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],0.90,0.80,16.0/9.0,title=title,pngDir=pngDir,
      twocol=True)
      title="awaveletsfg"
      plotWavelets(Sampling(na,dt,ka*dt),[naw,nak],0.90,0.80,16.0/9.0,title=title,pngDir=pngDir,
      twocol=True)
      title="hKnownwaveletsfg"
      plotWavelet(Sampling(nh,dt,kh*dt),nhk,0.90,0.70,16.0/9.0,title=title,pngDir=pngDir,
      twocol=True)

def normalize(h):
  #return div(h,max(max(h),-min(h)))
  return div(h,rms(h))

def rms(h):
  return sqrt(sum(mul(h,h))/len(h))


def makeImpulses(r,nt,ni):
  p = zerofloat(nt)
  q = zerofloat(nt)
  tmax = nt-1
  if r<=1.0:
    ts = rampfloat(tmax/(ni+1),tmax/(ni+1),ni)
  else:
    ts = rampfloat(tmax/(ni+1)/r,tmax/(ni+1)/r,ni)
  si = SincInterp.fromErrorAndFrequency(0.01,0.45)
  rj = -1.0
  for ji in range(ni):
    tj = ts[ji]
    tp = ts[ji]
    tq = r*tp
    if r<=1.0:
      tq += (1.0-r)*tmax
    #print "tp =",tp," tq =",tq
    rj = -rj
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)
  return p,q

def addWavelet(fpeak,decay,p,mp=False):
  w = 2.0*PI*fpeak
  print "w = "+str(w)
  if not mp:
    decay *= 2.0
    print "decay = "+str(decay)
    w -= 2.0*PI*0.04
    print "w= "+str(w)
  r = exp(-decay)
  print "r= "+str(r)
  a1,a2 = -2.0*r*cos(w),r*r
  print "a =",[1,a1,a2]
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  print "pole1 = "+poles[0].toString()
  print "pole2 = "+poles[1].toString()
  zeros = []
  gain = 1.0
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(p,t)
  if not mp:
    w = 2.0*PI*(fpeak+0.04)
    a1r,a2r = -2.0*r*cos(w),r*r
    print "areverse =",[1,a1r,a2r]
    print "wreverse= "+str(w)
    poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
    print "pole1reverse = "+poles[0].toString()
    print "pole2reverse = "+poles[1].toString()
    zeros = []
    gain = 1.0
    rcf = RecursiveCascadeFilter(poles,zeros,gain)
    rcf.applyReverse(t,x)
  else:
    copy(t,x)
  conv(2,0,[1.0,-0.95],len(x),0,copy(x),len(x),0,x) # attenuate DC
  return x

def getWavelet(fpeak,decay,nh,kh,mp=False):
  x = zerofloat(nh)
  x[-kh] = 1.0
  return addWavelet(fpeak,decay,x,mp)

def plotSequence(x,xmax=None,title=None):
  sp = SimplePlot.asPoints(x)
  if xmax==None:
    xmax = max(abs(max(x)),abs(min(x)))
    xmax *= 1.05
  sp.setVLimits(-xmax,xmax)
  if title:
    sp.setTitle(title)

def plotSequences(st,xs,labels=None,title=None):
  nx = len(xs)
  pp = PlotPanel(nx,1)
  for ix,xi in enumerate(xs):
    pv = pp.addPoints(ix,0,st,xi)
    if labels:
      pp.setVLabel(ix,labels[ix])
  pp.setHLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setVisible(True)
  if title:
    pp.setTitle(title)

def plotWavelets(st,hs,fracWidth,fracHeight,aspectRatio,hmax=None,title=None,pngDir=None,
  onecol=None,twocol=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  lw = 1.8,3,1
  nh = len(hs)
  hsmax = 0
  for ih in range(nh):
    if ih==0:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.NONE)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(3.4)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
    elif hs[ih]:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.SOLID)
      pv.setLineWidth(1)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
  if hmax==None:
    hmax = hsmax*1.05
  sp.setVLimits(-hmax,hmax)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Amplitude (normalized)")
  sp.setSize(960,560)
  if title:
    if pngDir==None:
      sp.setTitle(title)
  if pngDir:
    if onecol:
      sp.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"onecol.png"
      sp.paintToPng(720.0,3.08,pngDir)
    if twocol:
      sp.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"twocol.png"
      sp.paintToPng(720.0,3.0,pngDir)

def plotWavelet(st,h,fracWidth,fracHeight,aspectRatio,hmax=None,title=None,pngDir=None,
  onecol=None,twocol=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  lw = 1.8,3,1
  hsmax = 0
  pv = sp.addPoints(st,h)
  pv.setLineStyle(PointsView.Line.SOLID)
  pv.setLineWidth(1)
  hsmax = max(hsmax,abs(max(h)),abs(min(h)))
  if hmax==None:
    hmax = hsmax*1.05
  sp.setVLimits(-hmax,hmax)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Amplitude (normalized)")
  sp.setSize(960,560)
  if title:
    if pngDir==None:
      sp.setTitle(title)
  if pngDir:
    if onecol:
      sp.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"onecol.png"
      sp.paintToPng(720.0,3.08,pngDir)
    if twocol:
      sp.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"twocol.png"
      sp.paintToPng(720.0,3.0,pngDir)

def plot1Trace(st, f, tmin, tmax, 
  amax, fracWidth, fracHeight, aspectRatio, title=None, pngDir=None, onecol=None, twocol=None):
  pv1 = PointsView(st,f)
  pv1.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  
  pp = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.setHLimits(0,-amax,amax)
  pp.setHInterval(0,2.0)
  pp.setVLimits(tmin,tmax)
  pp.setHLabel(0,"Amplitude")
  pp.setVLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if onecol:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      print "zz = "+str(aspectRatio)
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"twocol.png"
      pf.paintToPng(720.0,3.0,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plot2TracesSideBySide(st, f, g, tmin, tmax, 
  amax, fracWidth, fracHeight, aspectRatio, title=None, pngDir=None, onecol=None, twocol=None):
  pv1 = PointsView(st,f)
  pv1.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  pv2 = PointsView(st,g)
  pv2.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  
  pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(0,1,pv2)
  pp.setHLimits(0,-amax,amax)
  pp.setHLimits(1,-amax,amax)
  pp.setHInterval(0,2.0)
  pp.setHInterval(1,2.0)
  pp.setVLimits(tmin,tmax)
  pp.setHLabel(0,"Amplitude")
  pp.setHLabel(1,"Amplitude")
  pp.setVLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if onecol:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      print "zz = "+str(aspectRatio)
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"twocol.png"
      pf.paintToPng(720.0,3.0,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plot3TracesSideBySide(st, f, g, h, tmin, tmax, 
  amax, fracWidth, fracHeight, aspectRatio, title=None, pngDir=None, onecol=None, twocol=None):
  pv1 = PointsView(st,f)
  pv1.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  pv2 = PointsView(st,g)
  pv2.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  pv3 = PointsView(st,h)
  pv3.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  
  pp = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(0,1,pv2)
  pp.addTiledView(0,2,pv3)
  pp.setHLimits(0,-amax,amax)
  pp.setHLimits(1,-amax,amax)
  pp.setHLimits(2,-amax,amax)
  pp.setHInterval(0,2.0)
  pp.setHInterval(1,2.0)
  pp.setHInterval(2,2.0)
  pp.setVLimits(tmin,tmax)
  pp.setHLabel(0,"Amplitude")
  pp.setHLabel(1,"Amplitude")
  pp.setHLabel(2,"Amplitude")
  pp.setVLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if onecol:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      print "zz = "+str(aspectRatio)
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"twocol.png"
      pf.paintToPng(720.0,3.0,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def addNoise(nrms, seed, f):
  n = len(f)
  r = Random(seed)
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n),0.5))
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
  
def convert1Dto2D(f1,g1,u1):
  nt = len(f1)
  nx = 200
  f = zerofloat(nt,nx)
  g = zerofloat(nt,nx)
  for i in range(0,nx):
    f[i] = f1
    g[i] = g1

  u = zerofloat(nt,nx)
  for i in range(0,nx):
    u[i] = u1
  return f,g,u





#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
