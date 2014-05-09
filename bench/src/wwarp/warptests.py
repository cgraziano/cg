from imports import *
from edu.mines.jtk.dsp.Conv import *
from wwarp import Warp, WaveletWarping
from testing import Synthetic
from linalgebra import ToeplitzRecursion

#pngDir = "./png/figures/"
pngDir = None

def main(args):
  enhancedWarping()

def enhancedWarping():
  ##################Parameters##############################
  #Stretch or squeeze coefficient
  r = .5
  #Time window for wavelet estimation
  tmin,tmax = 0,1.92
  #Frequency window for wavelet estimation
  fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
  #Inverse wavelet parameters estimated
  na,ka = 81,-20
  #Wavelet parameters estimated
  nh,kh = 181,-90
  #Peak frequency and decay of wavelet
  fpeak,decay = 0.08,0.05
  #Number of impulses in each trace
  ni = 2
  #Sampling of data
  nt,dt,ft = 481,.004,0.000
  st = Sampling(nt,dt,ft)
  #Time window for display
  dtmin,dtmax = tmin,tmax
  #amax for display
  damax = 3.5
  ##########################################################
  mp = False#Minimum Phase

  hf,hg = makeTraces(st,r,ni,fpeak,decay,mp=mp)
  u = rampfloat(0.0,r,nt)
  print "u before mod"
  dump(u)
  if r<=1.0:
    u = add((1.0-r)*(nt-1),u)


  title = "StartingWarpedTraces"
  plot2TracesSideBySide(st,hf,hg,dtmin,dtmax,damax,title=title,
  pngDir=pngDir,onecol=True)

  dt = st.getDelta()
  odt = 1.0/dt
  warp = Warp()
  ww = WaveletWarping(warp)
  print int(tmin*odt)
  print int(tmax*odt)
  ww.setTimeRange(int(tmin*odt),int(tmax*odt))
  ww.setStabilityFactor(1.00)

  ww.setFrequencyRange(fmin,fmax)
  a = ww.estimateInverseAtimes(na,ka,st,st,u,hf,hg)

  """ 
  bdf,bdg,bsdg,d,bd = True,False,False,False,False
  df0 = ww.getDifferencePlots(0,bdf,bdg,bsdg,d,bd)
  df1 = ww.getDifferencePlots(1,bdf,bdg,bsdg,d,bd)
  df2 = ww.getDifferencePlots(2,bdf,bdg,bsdg,d,bd)
  title = "Delayed f"
  plot3TracesSideBySide(st,df0,df1,df2,dtmin,dtmax,
  damax,title=title,pngDir=pngDir,onecol=True)

  bdf,bdg,bsdg,d,bd = False,False,True,False,False
  sdg0 = ww.getDifferencePlots(0,bdf,bdg,bsdg,d,bd)
  sdg1 = ww.getDifferencePlots(1,bdf,bdg,bsdg,d,bd)
  sdg2 = ww.getDifferencePlots(2,bdf,bdg,bsdg,d,bd)
  title = "Delayed and warped g"
  plot3TracesSideBySide(st,sdg0,sdg1,sdg2,dtmin,dtmax,
  damax,title=title,pngDir=pngDir,onecol=True)
  
  bdf,bdg,bsdg,d,bd = False,False,False,True,False
  d0 = ww.getDifferencePlots(0,bdf,bdg,bsdg,d,bd)
  d1 = ww.getDifferencePlots(1,bdf,bdg,bsdg,d,bd)
  d2 = ww.getDifferencePlots(2,bdf,bdg,bsdg,d,bd)
  title = "DifferenceTraces"
  plot3TracesSideBySide(st,d0,d1,d2,dtmin,dtmax,
  damax,title=title,pngDir=pngDir,onecol=True)
  title = "DifferenceTraceszoomed"
  plot3TracesSideBySide(st,d0,d1,d2,dtmin,dtmax,
  damax,title=title,pngDir=pngDir,onecol=True,twocol=True)
  """

  h = ww.getWaveletH(na,ka,a,nh,kh)
  hsyn = getWavelet(fpeak,decay,nh,kh,mp) # synthetic/known wavelet
  
  h = normalize(h)
  hsyn = normalize(h)

  title = "Estimated Wavelet_"+str(r)+"_time"
  plotWavelets(Sampling(nh,st.getDelta(),kh*st.getDelta()),
    [hsyn,h],title=title,pngDir=pngDir,onecol=True,twocol=True)
  
  #ahg = WaveletWarping.applyFilter(na,ka,a,hg)
  #slaf = WaveletWarping.applySLAg(na,ka,a,r,st,hg)
  #title = "Impulse_Warpt"
  #plot2TracesSideBySide(st,ahg,slaf,dtmin,dtmax,
  #damax,title=title,pngDir=pngDir,onecol=True)
#
#  hslag = WaveletWarping.applyHLSAg(nh,kh,h,na,ka,a,r,st,hg)
#  title = "FinalComparison"
#  plot3TracesSideBySide(st,hf,hslag,sdg0,dtmin,dtmax,
#  damax,title=title,pngDir=pngDir,onecol=True,twocol=True)
#  #ww.plotAmplitudeSpectrum(Sampling(nh,st.getDelta(),kh*st.getDelta()),
#  #h,False,"h Amplitude Spectrum");
#  #ww.plotAmplitudeSpectrum(st,hg,False,"hg");
#  print "r = "+str(r)
#  dump(a)

def enhancedWarpigNoise():
  #gauss distribution
  print "not done"

def enhancedWarpigNonExactWarping():
  print "not done"

def enhancedWarpigDiffWavelets():
  print "not done"

def makeTraces(st, r, ni, fpeak, decay, mp=None):
  nt = st.getCount()
  f,g = makeImpulses(r,nt,ni)
  hf = addWavelet(fpeak,decay,f,mp)
  hg = addWavelet(fpeak,decay,g,mp)
  return hf,hg

def makeReflections(st, rt, ra, r):
  nt = st.getCount()
  dt = st.getDelta()
  ft = st.getFirst()
  nref = len(rt)
  rt1 = rt
  rt2 = mul(r,rt)
  ra1 = ra
  ra2 = ra
  f = zerofloat(nt)
  g = zerofloat(nt)
  si = SincInterp.fromErrorAndFrequency(0.01, 0.45)
  for ir in range(0,nref):
    rti = rt1[ir]
    rai = ra1[ir]
    si.accumulate(rti,rai,nt,dt,ft,f)
    rti = rt2[ir]
    rai = ra2[ir]
    si.accumulate(rti,rai,nt,dt,ft,g)
  return f,g

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


def getArWavelet(st,fpeak,decay,nh,kh,mp=False):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  print "a1 =",a1," a2 =",a2
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  if mp:
    gain *= sqrt(1.0+a1*a1+a2*a2)
  x = zerofloat(nh)
  t = zerofloat(nh)
  h = zerofloat(nh)
  x[-kh] = 1.0
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(x,t)
  if mp:
    rcf.applyReverse(t,h)
  else:
    copy(t,h)
  return h

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


def plot2TracesSideBySide(st, f, g, tmin, tmax, 
  amax, title=None, pngDir=None, onecol=None, twocol=None):
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
  pp.setVLimits(tmin,tmax)
  pp.setHLabel(0,"Amplitude")
  pp.setHLabel(1,"Amplitude")
  pp.setVLabel("Time (s)")
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

def plot3TracesSideBySide(st, f, g, h, tmin, tmax, 
  amax, title=None, pngDir=None, onecol=None, twocol=None):
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
  pp.setVLimits(tmin,tmax)
  pp.setHLabel(0,"Amplitude")
  pp.setHLabel(1,"Amplitude")
  pp.setHLabel(2,"Amplitude")
  pp.setVLabel("Time (s)")
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

def plotTrace(st, p, tmin, tmax, amax, title): 
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  sp.setVLimits(tmin,tmax)
  sp.setHLimits(-amax,amax)
  sp.addTitle(title)
  pv = sp.addPoints(st,p)

def plotWavelets(st,hs,hmax=None,title=None,pngDir=None,
  onecol=None,twocol=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  lw = 1,3,1
  nh = len(hs)
  hsmax = 0
  for ih in range(nh):
    if hs[ih]:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(ls[ih])
      pv.setLineWidth(lw[ih])
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
  if hmax==None:
    hmax = hsmax*1.05
  sp.setVLimits(-hmax,hmax)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Amplitude")
  if title:
    if pngDir==None:
      sp.setTitle(title)
  if pngDir:
    if onecol:
      sp.setFontSizeForPrint(8.0,222.0)
      pngDir = pngDir+title+"onecol.png"
      sp.paintToPng(720.0,3.08,pngDir)
    if twocol:
      sp.setFontSizeForPrint(8.0,469.0)
      pngDir = pngDir+title+"twocol.png"
      sp.paintToPng(720.0,6.51,pngDir)

def balance(sigma,f):
  f = add(max(f)*0.00001,f)
  ff = mul(f,f)
  RecursiveExponentialFilter(sigma).apply1(ff,ff)
  return div(f,sqrt(ff))

def normalize(h):
  return div(h,max(max(h),-min(h)))





#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())


