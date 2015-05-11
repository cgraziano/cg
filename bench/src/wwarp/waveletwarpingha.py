#############################################################################
# Demo wavelet estimation from warping.
# This version minimizes || f - HSLAg ||, using iterative alternating
# updates to H and A, where HA need not equal the identity operator I.
# (1) find h such that f ~ SLAg * h (where h is a shaping filter)
# (2) find a such that f ~ HSLGa (solving non-Toeplitz system)

from imports import *

from edu.mines.jtk.dsp.Conv import *
from wwarp import WaveletWarpingHA

#############################################################################

#pngDir = "./png/figures/"
#pngDir = "./png/figuresAbs/"
#pngDir = "./png/200iterationsh1/"
pngDir = None

def main(args):
  #goSimpleTest()
  goSino()

def goSino():
  na,ka = 81,-40 # sampling for inverse A of wavelet in PS image
  nh,kh = 81,-40 # sampling for wavelet H in PP image
  nt,dt,ft = 501,0.004,0.000 # used for plotting only
  nx,dx,fx = 721,0.015,0.000
  sa = Sampling(na,dt,ka*dt) # sampling of the inverse wavelet in the PS image
  sh = Sampling(nh,dt,kh*dt) # sampling of the wavelet in the PP image
  st = Sampling(nt,dt,ft)  # sampling of the PP image
  ntg = 852
  stg = Sampling(ntg,dt,ft) # sampling of the PS image
  sx = Sampling(nx,dx,fx) # spatial sampling for both the PP and PS images
  itmin,itmax = 100,500 # PP time window
  sfac = 1.000 # stabilization factor
  wha = 0.000 # weight for HA = I terms
  f,g,u = getSinoImages() # PP image, PS image, and warping u(t,x)

  
  ntwin = nt
  uprime = zerofloat(ntwin,nx)
  for ix in range(0,nx):
    for it in range(1,ntwin):
      uprime[ix][it] = u[ix][it]-u[ix][it-1]
  SimplePlot.asPixels(u)
  SimplePlot.asPoints(u[50])
  SimplePlot.asPoints(uprime[50])
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(uprime)
  pv.setClips(1.4,1.72)
  sp.addColorBar()

  plotAmplitudeSpectrumT(st, f[300], 100, 400, "f amp", amax=None)#Amplitude spectrum single trace
  plotAmplitudeSpectrumT(stg, g[300], 228, 700, "g amp", amax=None)

  ww = WaveletWarpingHA()
  ww.setTimeRange(itmin,itmax)# Identify the time range that is used for wavelet estimation.
  ww.setStabilityFactor(sfac)# Stability factor if matrix is not positive definite. (unsure if 
  #best solution)
  slg = ww.applyS(u,ww.applyL(u,g)) # PS warping without wavelets
  plotImageTime(st,sx,f,vlabel="PP time (s)",itmin=100,itmax=400,title="PPsamples",pngDir=pngDir,\
  halfcol=False,onecol=True,twocol=False)
  plotImageTime(stg,sx,g,vlabel="PS time (s)",itmin=228,itmax=700,title="PSsamples",\
  pngDir=pngDir,halfcol=False,onecol=True,twocol=False)
  e1 = ww.rms(sub(f,slg))
  for niter in [0,11]:
    for wha in [0]:#wha controls the influence of making h and a inverses of one another.
      print "niter =",niter," wha =",wha
      suffix = str(niter)+str(int(wha))
      ww.setWeightHA(wha)
      
      ag = zerofloat(na); ag[-ka] = 1.0 # initial inverse a in g
      #ag = mul(ag,10)
      hf = ww.getWaveletH(nh,kh,na,ka,ag,u,f,g) # wavelet h in f
      #hf = mul(hf,10)
      hslag = ww.applyHSLA(na,ka,ag,nh,kh,hf,u,g) # PS warping wavelets
      ew = ww.rms(sub(f,hslag))
      print "ew",ew
      print "ag:"
      dump(ag)
      print "hf:"
      dump(hf)
      #main iteration loop and what is needed to get estimate fo wavelets.
      for jiter in range(niter):
        print "iteration",jiter
        ag = ww.getInverseA(na,ka,nh,kh,hf,u,f,g) # inverse a in g
        hf = ww.getWaveletH(nh,kh,na,ka,ag,u,f,g) # wavelet h in f

        print "ag:"
        dump(ag)
        print "hf:"
        dump(hf)

        hslag = ww.applyHSLA(na,ka,ag,nh,kh,hf,u,g) # PS warping wavelets
        ew = ww.rms(sub(f,hslag))
        print "ew",ew
        hg = ww.getWaveletH(na,ka,ag,nh,kh) # wavelet in g
        title = "wavelets"+str(jiter)+" ew = "+str(ew)
        plotWaveletsPpPs(sh,hf,hg,title=title,pngDir=pngDir,\
        halfcol=True,onecol=False,twocol=False)
      
      aG = ww.applyA(na,ka,ag,g)
      plotImage(stg,sx,aG,zoom=False,png="Ag"+suffix)
      hg = ww.getWaveletH(na,ka,ag,nh,kh) # wavelet in g
      sg = ww.applyS(u,ww.applyL(u,g))
      hslag = ww.applyHSLA(na,ka,ag,nh,kh,hf,u,g) # PS warping with wavelets

      #hslag = mul(hslag,ww.rms(f)/ww.rms(hslag))
      ew = ww.rms(sub(f,hslag))
      print "  e1 =",e1," ew =",ew
      #plotImage(st,sx,hslag,zoom=True,png="psww"+suffix)

      title="wavelets"+suffix
      plotWaveletsPpPs(sh,hf,hg,title=title,pngDir=pngDir,\
      halfcol=True,onecol=False,twocol=False)
      title="wavelets"+suffix+"side"
      plotWaveletsPpPsSide(sh,hf,hg,title=title,pngDir=pngDir,\
      halfcol=False,onecol=True,twocol=False)

      
      #Make images within specified window have rms of 1.
      frms = ww.rms(f)
      sgrms = ww.rms(sg)
      frms1 = div(f,frms)
      sgrms1 = div(sg,sgrms)

      if niter==0:
        hslg0 = ww.applyHSLA(na,ka,ag,nh,kh,hf,u,g) # PS warping with wavelets
        #Make images within specified window have rms of 1.
        hslg0rms = ww.rms(hslg0)
        hslg0 = div(hslg0,hslg0rms)

      if niter!=0:
        hslag11 = ww.applyHSLA(na,ka,ag,nh,kh,hf,u,g) # PS warping with wavelets
        #Make images within specified window have rms of 1.
        hslag11rms = ww.rms(hslag11)
        hslag11 = div(hslag11,hslag11rms)

        ixmin=int(0/dx)
        ixmax=int(10.8/dx)
        itmin=int(.6/dt)
        itmax=int(1.0/dt)
        ticintV = 0.1
        ticintH = 1
        amax=2 
        title = "3Imagesteaser"+suffix#zoomed version of other image
        teaser = True
        plot3ImagesSideBySide(st,sx,frms1,hslag11,sgrms1,\
        itmin,itmax,ixmin,ixmax,ticintV,ticintH,amax,teaser,title=title,\
        pngDir=pngDir, onecol=None, twocol=True)

        title = "3ImagesUnzoomedTeaser"+suffix#zoomed version of other image
        ixmin=0
        ixmax=nx-1
        itmin=100
        itmax=400
        ticintV = 0.2
        ticintH = 2
        teaser = False
        plot3ImagesSideBySide(st,sx,f,sg,hslag11,\
        itmin,itmax,ixmin,ixmax,ticintV,ticintH,amax,teaser,title=title,\
        pngDir=pngDir,onecol=None,twocol=True)


  itmax = 212 
  title="4imageComp"
  plot4ImagesSideBySide(st,sx,frms1,hslag11,hslg0,sgrms1,\
  itmin,itmax,ixmin,ixmax,amax,\
  title=title, pngDir=pngDir, onecol=None, twocol=True)

  
#synthetic tests of iterative wavelet estimation.
def goSimpleTest():
  nt,ni = 481,10 # number of time samples; number of impulses
  freq,decay = 0.08,0.05 # peak frequency and decay for wavelet
  #na,ka = 11,-5 # sampling for inverse wavelet A
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.00
  wha = 10.0 
  st = Sampling(nt,dt,ft)
  for mp in [False]: # True, for minimum-phase; False for other
    hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
    for r in [2.0]: # 0.5 for stretch; 2.0 for squeeze
      aw = zerofloat(na); aw[-ka] = 1.0
      hw = zerofloat(nh); hw[-kh] = 1.0
      p,q = makeImpulses(r,nt,ni)
      f = addWavelet(freq,decay,p,mp)
      g = addWavelet(freq,decay,q,mp)
      u = rampfloat(0.0,r,nt)
      if r<=1.0:
        u = add((1.0-r)*(nt-1),u)
      ww = WaveletWarpingHA()
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)
      ww.setWeightHA(wha)
      ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
      #dump(ak)
      hw = ww.getWaveletH(nh,kh,na,ka,aw,u,f,g) # estimated wavelet
      for iter in range(1):
        aw = ww.getInverseA(na,ka,nh,kh,hw,u,f,g) # estimated inverse
        print "aw:"
        dump(aw)
        hw = ww.getWaveletH(nh,kh,na,ka,aw,u,f,g) # estimated wavelet
        print "hw:"
        dump(hw)
        hslag = ww.applyHSLA(na,ka,aw,nh,kh,hw,u,g) # PS warping wavelets
        ew = ww.rms(sub(f,hslag))
        print "ew",ew

        #dump(aw)
      #hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
      sg = ww.applyS(u,g)
      ag = ww.applyA(na,ka,aw,g)
      lag = ww.applyL(u,ag) # lowpass, if squeezing
      slag = ww.applyS(u,lag)
      hslag = ww.applyH(nh,kh,hw,slag)
      normalizeMax(hw)
      normalizeMax(hk)
      title = "r = "+str(r)
      #plotSequences(st,[f,g],labels=["f","g"],title=title)
      #plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[f,hslag],labels=["f","HSLAg"],title=title)
      plotWavelets(Sampling(nh,dt,kh*dt),[hw,hk],title=title)

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

def plotImageTime(st,sx,f,vlabel,itmin=None,itmax=None,title=None,pngDir=None,
  halfcol=None,onecol=None,twocol=None):
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

  if title:
    if pngDir==None:
      sp.setTitle(title)
  if pngDir:
    if halfcol:
      sp.setFontSizeForPrint(8.0,111.0)
      pngDir = pngDir+title+"halfcol.png"
      sp.paintToPng(720.0,1.54,pngDir)
    if onecol:
      sp.setFontSizeForPrint(8.0,222.0)
      pngDir = pngDir+title+"onecol.png"
      sp.paintToPng(720.0,3.23,pngDir)
    if twocol:
      sp.setFontSizeForPrint(8.0,469.0)
      pngDir = pngDir+title+"twocol.png"
      sp.paintToPng(720.0,6.51,pngDir)

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


def plotWaveletsPpPs(st,h1,h2,title=None,
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
  pf.setSize(175,450)
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

def plotSpectrum(sf,spec,title,amax=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setVLabel("Amplitude")
  sp.setHLabel("Frequency (Hz)")
  sp.setSize(750,400)
  sp.addTitle(title)
  if amax:
    sp.setVLimits(0,amax)
  pv = sp.addPoints(sf,spec)



#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
