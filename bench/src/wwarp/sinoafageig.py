#############################################################################
#Builds PP and PS figures for CWP presentation.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from wwarp import WaveletWarpingAfAgEig

#############################################################################
#pngdir is the location the images will be sent to.
#pngDir = "./png/figures/"
#pngDir = "./png/figuresAbs/"
#pngDir = "./pres14/PPPSSTiShAH/"
#pngDir = "./pres14/PPPSSTiShAH/poster/"
#pngDir = "./png/200iterationsh1/"
pngDir = None

def main(args):
  goSino()

def goSino():
  ####Sampling parameters
  naf,kaf = 20,-10 # sampling for inverse A of wavelet in PS image
  nag,kag = 20,-10 # sampling for inverse A of wavelet in PS image
  nhf,khf = 81,-40 # sampling for wavelet H in PP image
  nhg,khg = 81,-40 # sampling for wavelet H in PP image
  nt,dt,ft = 501,0.004,0.000 # used for plotting only
  nx,dx,fx = 721,0.015,0.000
  saf = Sampling(naf,dt,kaf*dt)
  sag = Sampling(nag,dt,kag*dt)
  shf = Sampling(nhf,dt,khf*dt)
  shg = Sampling(nhg,dt,khg*dt)
  st = Sampling(nt,dt,ft)
  ntg = 852
  stg = Sampling(ntg,dt,ft)
  sx = Sampling(nx,dx,fx)
  itmin,itmax = 100,400 # PP time window
  sfac = 1.000 # stabilization factor
  wha = 0.000 # weight for HA = I terms
  ####

  ####Plotting of initial images
  ## f represents the pp image, g represents the ps image, and ts represents the pre-calculated
  #time shifts from DTW. 
  rawf,rawg,ts = getSinoImagesRaw()
  f,g,u = getSinoImages() # PP image, PS image, and warping u(t,x)
  title = "rawPP"
  plotImageTimeHalf(st,sx,rawf,"PP time (s)",0.45,0.8,16.0/9.0,itmin=100,
  itmax=400,title=title,pngDir=pngDir,twocol=True)
  title = "rawPS"
  plotImageTimeHalf(stg,sx,rawg,"PS time (s)",0.45,0.8,16.0/9.0,itmin=228,
  itmax=700,title=title,pngDir=pngDir,twocol=True)
  title = "PP"
  plotImageTimeHalf(st,sx,f,"PP time (s)",0.45,0.8,16.0/9.0,itmin=100,
  itmax=400,title=title,pngDir=pngDir,twocol=True)
  title = "PS"
  plotImageTimeHalf(stg,sx,g,"PS time (s)",0.45,0.8,16.0/9.0,itmin=228,
  itmax=700,title=title,pngDir=pngDir,twocol=True)
  title="g"
  plotImageTimeWhole(stg,sx,g,"PS time (s)",0.9,0.8,16.0/9.0,itmin=228,itmax=700,
  title=title,pngDir=pngDir,twocol=True)
  title="f"
  plotImageTimeWhole(st,sx,f,"PP time (s)",0.9,0.8,16.0/9.0,itmin=itmin,itmax=itmax,
  title=title,pngDir=pngDir,twocol=True)
  ####

  ####The amount of squeezing applied to the PS image.  
  ntwin = nt
  uprime = zerofloat(ntwin,nx)
  for ix in range(0,nx):
    for it in range(1,ntwin):
      #uprime[ix][it] = u[ix][itmin+it]-u[ix][itmin+it-1]
      uprime[ix][it] = u[ix][it]-u[ix][it-1]
  title = "uprime"
  plotImageTimeU(st,sx,uprime,"PP time (s)",0.9,0.8,16.0/9.0,itmin=itmin,itmax=itmax,
  title=title,pngDir=pngDir,twocol=True)
  ####
  
  ####This is the warping with wavelets section.
  ww = WaveletWarpingAfAgEig()
  ww.setTimeRange(itmin,itmax)
  ww.setStabilityFactor(sfac)
  afag = ww.getInverseAfAg(naf,kaf,nag,kag,u,f,g)
  awf,awg = separateafag(naf,nag,afag)
  #Eigenvalue and eigenvector analysis
  natot = naf+nag
  for i in range(1):
    print "eigenvector "+str(i)
    dump(ww.getEigVector(i))

  print "eigval0/eigval1 = "+str(ww.getEig01Ratio());
  print "eigval0 = "+str(ww.getEigVal(0));
  print "eigval1 = "+str(ww.getEigVal(1));

  hwf = ww.getWaveletH(naf,kaf,awf,nhf,khf)
  hwg = ww.getWaveletH(nag,kag,awg,nhg,khg)

  nhwf = normalizeMAAWOS(hwf)
  nhwg = normalizeMAAWOS(hwg)

  plotWaveletsPpPs(shf,nhwf,nhwg,0.9,0.8,16.0/9.0,title=title,pngDir=pngDir,twocol=True)

def getSinoImages():
  dataDir = "/Users/Chris/data/sinos/"
  n1f,n1g,d1,f1 = 501,852,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1f,n2)
  g = readImage(dataDir+"ps.dat",n1g,n2)
  u = readImage(dataDir+"shifts.dat",n1f,n2)
  u = add(u,rampfloat(0.0,1.0,0.0,n1f,n2))
  gain(100,f)
  gain(100,g)
  return f,g,u

def getSinoImagesRaw():
  dataDir = "/Users/Chris/data/sinos/"
  n1f,n1g,d1,f1 = 501,852,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1f,n2)
  g = readImage(dataDir+"ps.dat",n1g,n2)
  ts = readImage(dataDir+"shifts.dat",n1f,n2)
  return f,g,ts

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def gain(hw,f):
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def normalizeMax(f):
  return mul(1.0/max(abs(f)),f,f)
 
def normalizeMAAWOS(f):
  minf = min(f)
  maxf = max(f)
  absminf = abs(minf) 
  absmaxf = abs(maxf) 
  if absminf<absmaxf:
    return mul(1.0/maxf,f)
  else:
    return mul(1.0/minf,f)

def normalizeRms(f):
  mul(1.0/rms(f),f,f)
def rms(f):
  return sqrt(sum(mul(f,f))/len(f)/len(f[0]))

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
    #rjp = -rj
    #rjq = sin(2*PI*(1+ji)*rj)
    rjp = -rj
    rjq = -rj

    si.accumulate(tp,rjp,nt,1.0,0.0,p)
    si.accumulate(tq,rjq,nt,1.0,0.0,q)
  return p,q

def addWavelet(fpeak,decay,p,mp=False):
  w = 2.0*PI*fpeak
  if not mp:
    decay *= 2.0
    w -= 2.0*PI*0.04
  r = exp(-decay)
  a1,a2 = -2.0*r*cos(w),r*r
  #print "a =",[1,a1,a2]
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(p,t)
  if not mp:
    w = 2.0*PI*(fpeak+0.04)
    poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
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

def plotWavelets(st,hs,hmax=None,title=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  cs = [Color.BLACK,Color.GRAY,Color.LIGHT_GRAY]
  nh = len(hs)
  hsmax = 0
  for ih in range(nh):
    if hs[ih]:
      #pv = sp.addPoints(st,hs[ih])
      #pv.setLineStyle(ls[ih])
      sv = sp.addSequence(st,hs[ih])
      sv.setColor(cs[ih])
      #pv.setLineWidth(2)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
  if hmax==None:
    hmax = hsmax*1.05
  sp.setVLimits(-hmax,hmax)
  sp.setVLabel("Amplitude (normalized)")
  sp.setHLabel("Time (s)")
  if title:
    sp.setTitle(title)

def plotImage(st,sx,f,zoom=False,png=None):
  wpt = 240
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setFontSizeForPrint(8,wpt)
  sp.setSize(350,450)
  pv = sp.addPixels(st,sx,f)
  pv.setClips(-2.0,2.0)
  if zoom:
    sp.setVLimits(0.400,1.600) # zoom consistent with tmin,tmax = 100,400
  sp.setHLabel("Distance (km)")
  sp.setVLabel("PP time (s)")
  if pngDir and png:
    sp.paintToPng(360,3.08,pngDir+png+".png")
  else: 
    sp.setTitle(png)

def plotImageTimeHalf(st,sx,f,vlabel,fracWidth,fracHeight,aspectRatio,itmin=None,itmax=None,title=None,pngDir=None,halfcol=None,onecol=None,twocol=None):
  dt = st.getDelta()
  wpt = 240
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(350,450)
  pv = sp.addPixels(st,sx,f)
  pv.setClips(-2.0,2.0)
  sp.setVLabel(vlabel)
  sp.setHLabel("Distance (km)")
  sp.setVInterval(0.5)
  sp.setHInterval(2.0)
  if itmin:
    sp.setVLimits(itmin*dt,itmax*dt)

  sp.setSize(480,560)
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

def plotImageTimeWhole(st,sx,f,vlabel,fracWidth,fracHeight,aspectRatio,itmin=None,itmax=None,title=None,pngDir=None,
  halfcol=None,onecol=None,twocol=None):
  dt = st.getDelta()
  wpt = 240
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(st,sx,f)
  pv.setClips(-2.0,2.0)
  sp.setVLabel(vlabel)
  sp.setHLabel("Distance (km)")
  sp.setVInterval(0.1)
  sp.setHInterval(2.0)
  if itmin:
    sp.setVLimits(itmin*dt,itmax*dt)

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
      #sp.paintToPng(720.0,3.0,pngDir)
      sp.paintToPng(720.0,14.0,pngDir)

def plotImageTimeU(st,sx,f,vlabel,fracWidth,fracHeight,aspectRatio,itmin=None,itmax=None,title=None,pngDir=None,
  halfcol=None,onecol=None,twocol=None):
  dt = st.getDelta()
  wpt = 240
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(st,sx,f)
  pv.setColorModel(ColorMap.JET)
  pv.setClips(1.3,1.8)
  sp.addColorBar()
  sp.setVLabel(vlabel)
  sp.setHLabel("Distance (km)")
  sp.setVInterval(0.5)
  sp.setHInterval(2.0)
  if itmin:
    sp.setVLimits(itmin*dt,itmax*dt)

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

def plot2ImagesSideBySide(st, sx, f, g, itmin, itmax, 
  ixmin, ixmax, ticintV, ticintH, amax, fracWidth, fracHeight, aspectRatio,title=None, pngDir=None, onecol=None, twocol=None):
  dt = st.getDelta()
  dx = sx.getDelta()
  xmin = ixmin*dx
  xmax = ixmax*dx
  tmin = itmin*dt
  tmax = itmax*dt
  pv1 = PixelsView(st,sx,f)
  pv1.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv1.setClips(-2.0,2.0)
  pv2 = PixelsView(st,sx,g)
  pv2.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv2.setClips(-2.0,2.0)
  
  pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(0,1,pv2)
  pp.setVLimits(tmin,tmax)
  pp.setHLimits(0,xmin,xmax)
  pp.setHLimits(1,xmin,xmax)
  #pp.setHLabel(0,"Amplitude")
  #pp.setHLabel(1,"Amplitude")
  #pp.setHLabel(2,"Amplitude")
  #pp.setHLabel(3,"Amplitude")
  pp.setHLabel(0,"Distance (km)")
  pp.setHLabel(1,"Distance (km)")
  pp.setVLabel("PP time (s)")
  pp.setHInterval(0,ticintH)
  pp.setHInterval(1,ticintH)
  pp.setVInterval(ticintV)
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


def plot3ImagesSideBySide(st, sx, f, sg, hslag, itmin, itmax, 
  ixmin, ixmax, ticintV, ticintH, amax, teaser, title=None, pngDir=None, onecol=None, twocol=None):
  dt = st.getDelta()
  dx = sx.getDelta()
  xmin = ixmin*dx
  xmax = ixmax*dx
  tmin = itmin*dt
  tmax = itmax*dt
  pv1 = PixelsView(st,sx,f)
  pv1.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv1.setClips(-2.0,2.0)
  pv2 = PixelsView(st,sx,sg)
  pv2.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv2.setClips(-2.0,2.0)
  pv3 = PixelsView(st,sx,hslag)
  pv3.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv3.setClips(-2.0,2.0)
  
  pp = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(0,1,pv2)
  pp.addTiledView(0,2,pv3)
  pp.setVLimits(tmin,tmax)
  pp.setHLimits(0,xmin,xmax)
  pp.setHLimits(1,xmin,xmax)
  pp.setHLimits(2,xmin,xmax)
  #pp.setHLabel(0,"Amplitude")
  #pp.setHLabel(1,"Amplitude")
  #pp.setHLabel(2,"Amplitude")
  #pp.setHLabel(3,"Amplitude")
  pp.setHLabel(0,"Distance (km)")
  pp.setHLabel(1,"Distance (km)")
  pp.setHLabel(2,"Distance (km)")
  pp.setVLabel("PP time (s)")
  pp.setHInterval(0,ticintH)
  pp.setHInterval(1,ticintH)
  pp.setHInterval(2,ticintH)
  pp.setVInterval(ticintV)
  pf = PlotFrame(pp)
  if teaser:
    pf.setSize(500,166)
  else:
    pf.setSize(500,300)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if onecol:
      pf.setFontSizeForPrint(8.0,222.0)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      pf.setFontSizeForPrint(8.0,469.0)
      pngDir = pngDir+title+"twocol.png"
      pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plot4ImagesSideBySide(st, sx, f, sg, hslg, hslag, itmin, itmax, 
  ixmin, ixmax, amax, title=None, pngDir=None, onecol=None, twocol=None):
  dt = st.getDelta()
  dx = sx.getDelta()
  xmin = ixmin*dx
  xmax = ixmax*dx
  tmin = itmin*dt
  tmax = itmax*dt
  pv1 = PixelsView(st,sx,f)
  pv1.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv1.setClips(-2.0,2.0)
  pv2 = PixelsView(st,sx,sg)
  pv2.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv2.setClips(-2.0,2.0)
  pv3 = PixelsView(st,sx,hslg)
  pv3.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv3.setClips(-2.0,2.0)
  pv4 = PixelsView(st,sx,hslag)
  pv4.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv4.setClips(-2.0,2.0)
  
  pp = PlotPanel(1,4,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(0,1,pv2)
  pp.addTiledView(0,2,pv3)
  pp.addTiledView(0,3,pv4)
  pp.setVLimits(tmin,tmax)
  pp.setHLimits(0,xmin,xmax)
  pp.setHLimits(1,xmin,xmax)
  pp.setHLimits(2,xmin,xmax)
  pp.setHLimits(3,xmin,xmax)
  #pp.setHLabel(0,"Amplitude")
  #pp.setHLabel(1,"Amplitude")
  #pp.setHLabel(2,"Amplitude")
  #pp.setHLabel(3,"Amplitude")
  pp.setHLabel(0,"Distance (km)")
  pp.setHLabel(1,"Distance (km)")
  pp.setHLabel(2,"Distance (km)")
  pp.setHLabel(3,"Distance (km)")
  pp.setHInterval(0,2)
  pp.setHInterval(1,2)
  pp.setHInterval(2,2)
  pp.setHInterval(3,2)
  pp.setVInterval(.2)
  pp.setVLabel("PP time (s)")
  pf = PlotFrame(pp)
  pf.setSize(2345,458)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if onecol:
      pf.setFontSizeForPrint(8.0,222.0)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      pf.setFontSizeForPrint(8.0,469.0)
      pngDir = pngDir+title+"twocol.png"
      pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plot4ImagesOnTopOfEach(st, sx, f, hslag, hslg, sg, itmin, itmax, 
  ixmin, ixmax, amax, title=None, pngDir=None, onecol=None, twocol=None):
  dt = st.getDelta()
  dx = sx.getDelta()
  xmin = ixmin*dx
  xmax = ixmax*dx
  tmin = itmin*dt
  tmax = itmax*dt
  pv1 = PixelsView(st,sx,f)
  pv1.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv1.setClips(-2.0,2.0)
  pv2 = PixelsView(st,sx,hslag)
  pv2.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv2.setClips(-2.0,2.0)
  pv3 = PixelsView(st,sx,hslg)
  pv3.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv3.setClips(-2.0,2.0)
  pv4 = PixelsView(st,sx,sg)
  pv4.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv4.setClips(-2.0,2.0)
  
  pp = PlotPanel(4,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(1,0,pv2)
  pp.addTiledView(2,0,pv3)
  pp.addTiledView(3,0,pv4)
  pp.setVLimits(0,tmin,tmax)
  pp.setVLimits(1,tmin,tmax)
  pp.setVLimits(2,tmin,tmax)
  pp.setVLimits(3,tmin,tmax)
  pp.setHLimits(xmin,xmax)
  #pp.setHLabel(0,"Amplitude")
  #pp.setHLabel(1,"Amplitude")
  #pp.setHLabel(2,"Amplitude")
  #pp.setHLabel(3,"Amplitude")
  pp.setHLabel("Distance (km)")
  pp.setVLabel(0,"PP time (s)")
  pp.setVLabel(1,"PP time (s)")
  pp.setVLabel(2,"PP time (s)")
  pp.setVLabel(3,"PP time (s)")
  pf = PlotFrame(pp)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if onecol:
      pf.setFontSizeForPrint(8.0,222.0)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      pf.setFontSizeForPrint(8.0,469.0)
      pngDir = pngDir+title+"twocol.png"
      pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plot4Images4Square(st, sx, f, hslag, hslg, sg, tmin, tmax, 
  xmin, xmax, amax, title=None, pngDir=None, onecol=None, twocol=None):
  dt = st.getDelta()
  dx = sx.getDelta()
  xmin = ixmin*dx
  xmax = ixmax*dx
  tmin = itmin*dt
  tmax = itmax*dt
  pv1 = PixelsView(st,sx,f)
  pv1.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv1.setClips(-2.0,2.0)
  pv2 = PixelsView(st,sx,hslag)
  pv2.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv2.setClips(-2.0,2.0)
  pv3 = PixelsView(st,sx,hslg)
  pv3.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv3.setClips(-2.0,2.0)
  pv4 = PixelsView(st,sx,sg)
  pv4.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv4.setClips(-2.0,2.0)
  
  pp = PlotPanel(2,2,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(1,0,pv2)
  pp.addTiledView(0,1,pv3)
  pp.addTiledView(1,1,pv4)
  pp.setVLimits(0,tmin,tmax)
  pp.setVLimits(1,tmin,tmax)
  pp.setHLimits(0,xmin,xmax)
  pp.setHLimits(1,xmin,xmax)
  #pp.setHLabel(0,"Amplitude")
  #pp.setHLabel(1,"Amplitude")
  #pp.setHLabel(2,"Amplitude")
  #pp.setHLabel(3,"Amplitude")
  pp.setHLabel(0,"Distance (km)")
  pp.setHLabel(1,"Distance (km)")
  pp.setVLabel(0,"PP time (s)")
  pp.setVLabel(1,"PP time (s)")
  pf = PlotFrame(pp)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if onecol:
      pf.setFontSizeForPrint(8.0,222.0)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      pf.setFontSizeForPrint(8.0,469.0)
      pngDir = pngDir+title+"twocol.png"
      pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)


def plotWaveletsPpPs(st,h1,h2,fracWidth,fracHeight,aspectRatio,title=None,
  pngDir=None,halfcol=None,onecol=None,twocol=None):
  wpt = 240
  pp = PlotPanel(2,1)
  h1 = mul(h1,1.0/max(abs(h1)))
  h2 = mul(h2,1.0/max(abs(h2)))
  sv1 = pp.addSequence(0,0,st,h1)
  sv2 = pp.addSequence(1,0,st,h2)
  pp.setVLimits(0,-0.65,1.05)
  pp.setVLimits(1,-0.65,1.05)
  if pngDir:
    pp.setVLabel(0,"Amplitude")
    pp.setVLabel(1,"Amplitude")
  else:
    pp.setVLabel(0,"PP wavelet")
    pp.setVLabel(1,"PS wavelet")

  pp.setHLabel("Time (s)")
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
      #pf.paintToPng(720.0,3.0,pngDir)
      pf.paintToPng(720.0,14.0,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plotWaveletsPpPsSide(st,h1,h2,title=None,
  pngDir=None,halfcol=None,onecol=None,twocol=None):
  wpt = 240
  pp = PlotPanel(1,2)
  h1 = mul(h1,1.0/max(abs(h1)))
  h2 = mul(h2,1.0/max(abs(h2)))
  sv1 = pp.addSequence(0,0,st,h1)
  sv2 = pp.addSequence(0,1,st,h2)
  pp.setVLimits(0,-0.65,1.05)
  if pngDir:
    pp.setVLabel(0,"Amplitude")
  else:
    pp.setVLabel(0,"PP wavelet")

  pp.setHLabel(0,"Time (s)")
  pp.setHLabel(1,"Time (s)")
  pf = PlotFrame(pp)
  pf.setSize(350,200)
  pf.setFontSizeForPrint(8,wpt)
  pf.setVisible(True)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if halfcol:
      pf.setFontSizeForPrint(8.0,111.0)
      pngDir = pngDir+title+"halfcol.png"
      pf.paintToPng(720.0,1.54,pngDir)
    if onecol:
      pf.setFontSizeForPrint(8.0,222.0)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      pf.setFontSizeForPrint(8.0,469.0)
      pngDir = pngDir+title+"twocol.png"
      pf.paintToPng(720.0,6.51,pngDir)

  #if pngDir and png:
    #pf.paintToPng(720,wpt/72.0,pngDir+png+".png")
#Plot amplitude spectrum for a single trace
def plotAmplitudeSpectrumT(st, p, itmin, itmax, title, amax=None):
  #Time sampling for the time window specified.
  nt = itmax-itmin
  dt = st.getDelta()
  ft = st.getValue(itmin)

  subp = zerofloat(nt)
  subst = Sampling(nt,dt,ft)
  for it in range(0,nt):
    subp[it] = p[itmin+it]

  #Frequency sampling
  nfft = FftReal.nfftSmall(4*nt)#more time sample, the finer freq. samples
  nf = nfft/2+1
  df = 1.0/(nfft*dt)
  ff = 0.0
  fs = Sampling(nf,df,ff)
  amp = computeAmplitudeSpectrum(subst,fs,nfft,subp)
  plotSpectrum(fs,amp,title,amax=amax)

def computeAmplitudeSpectrum(st, fs, nfft, p):
  nt = st.getCount()
  dt = st.getDelta()
  ft = st.getFirst()
  nf = fs.getCount()
  df = fs.getDelta()
  ff = fs.getFirst()

  # Real-to-complex fast Fourier transform.
  fft = FftReal(nfft)
  cf = zerofloat(2*nf)
  copy(nt,p,cf)
  fft.realToComplex(-1,cf,cf)

  #Adjust phase for possibly non-zero time of first sample.
  wft = rampfloat(0.0,-2.0*FLT_PI*df*ft,nf);
  cf = cmul(cf,cmplx(cos(wft),sin(wft)));

  af = cabs(cf);
  #Amplitude spectrum normalized
  #float amax = max(max(af),FLT_EPSILON);
  #af = mul(1.0f/amax,af);
  return af;

def calcMax(x):
  nt = len(x[0])
  nx = len(x)
  max = 0.0
  for ix in range(nx):
    for it in range(nt):
      if max < x[ix][it]:
        max = x[ix][it]
  return max

def max(x):
  nt = len(x)
  max = 0.0
  for it in range(nt):
    if max < x[it]:
      max = x[it]
  return max

def separateafag(naf,nag,afag):
  na = naf+nag
  af = zerofloat(naf)
  ag = zerofloat(nag)
  for i in range(0,naf):
    af[i] = afag[i]
  ii=0
  for i in range(naf,na):
    ag[ii] = afag[i]
    ii = ii + 1
  return af,ag

def plotSpectrum(sf,spec,title,amax=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setVLabel("Amplitude")
  sp.setHLabel("Frequency (Hz)")
  sp.setSize(750,400)
  sp.addTitle(title)
  if amax:
    sp.setVLimits(0,amax)
  pv = sp.addPoints(sf,spec)

def plotRMSDiff(st,ea,fracWidth,fracHeight,aspectRatio,hmin,hmax,color=None,title=None,pngDir=None,
  onecol=None,twocol=None):
  sp = SimplePlot()
  nh = len(ea)
  pv = sp.addPoints(st,ea)
  pv.setLineStyle(PointsView.Line.SOLID)
  pv.setLineWidth(2)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(20)
  if color:
    pv.setLineColor(color)
    pv.setMarkColor(color)
  sp.setVLimits(hmin,hmax)
  sp.setHLabel("Iteration")
  sp.setVLabel("RMS differences")
  
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

def plotRMSDiff3(st,ea,eb,ec,fracWidth,fracHeight,aspectRatio,hmin,hmax,title=None,pngDir=None,
  onecol=None,twocol=None):
  sp = SimplePlot()
  nh = len(ea)
  pv = sp.addPoints(st,ea)
  color = Color.BLACK
  pv.setLineColor(color)
  pv.setMarkColor(color)
  pv.setLineStyle(PointsView.Line.SOLID)
  pv.setLineWidth(2)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(20)
  pv = sp.addPoints(st,eb)
  color = Color.RED
  pv.setLineColor(color)
  pv.setMarkColor(color)
  pv.setLineStyle(PointsView.Line.SOLID)
  pv.setLineWidth(2)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(20)
  pv = sp.addPoints(st,ec)
  color = Color.BLUE
  pv.setLineColor(color)
  pv.setMarkColor(color)
  pv.setLineStyle(PointsView.Line.SOLID)
  pv.setLineWidth(2)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(20)
  sp.setVLimits(hmin,hmax)
  sp.setHLabel("Iteration")
  sp.setVLabel("RMS differences")
  
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

def plotRMSDiff2(st,ea,eb,fracWidth,fracHeight,aspectRatio,hmin,hmax,color1=None,color2=None,title=None,pngDir=None,
  onecol=None,twocol=None):
  sp = SimplePlot()
  nh = len(ea)
  pv = sp.addPoints(st,ea)
  if color1:
    pv.setLineColor(color1)
    pv.setMarkColor(color1)
  pv.setLineStyle(PointsView.Line.SOLID)
  pv.setLineWidth(2)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(20)
  pv = sp.addPoints(st,eb)
  if color2:
    pv.setLineColor(color2)
    pv.setMarkColor(color2)
  pv.setLineStyle(PointsView.Line.SOLID)
  pv.setLineWidth(2)
  pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pv.setMarkSize(20)
  sp.setVLimits(hmin,hmax)
  sp.setHLabel("Iteration")
  sp.setVLabel("RMS differences")
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



#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
