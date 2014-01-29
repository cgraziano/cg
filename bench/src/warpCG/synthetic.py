from imports import *
from edu.mines.jtk.dsp.Conv import *
from dwarp import Warp, DynamicWarpingW 
from testing import Synthetic
from linalgebra import ToeplitzRecursion

pngDir = "./png/synthetic/"
#pngDir = None

def main(args):
  awarp1 = .3
  tmin,tmax = .1,3
  #for i in range(20):
  #  tmin = ftmin+i*.1
  estimateWaveletNoOSamp(awarp1,tmin,tmax)
  #estimateWaveletNoOSamp(awarp2)


def estimateWaveletNoOSamp(awarp,tmin,tmax):
  na,ka = 3,0
  nh,kh = 800,-125
  fpeak,decay = 20,0.02
  nt,dt,ft = 2000,.002,0
  st = Sampling(nt,dt,ft)
  rt = .2,1.2
  ra = 1,1
  hf,hg = makeTraces(st,rt,ra,awarp,fpeak,decay,zp=False)
  plotTrace(st,hf,0,3,3,"hf")
  plotTrace(st,hg,0,3,3,"hg")

  dt = st.getDelta()
  odt = 1.0/dt
  warp = Warp()
  dww = DynamicWarpingW(warp)
  dww.setTimeRange(int(tmin*odt),int(tmax*odt))
  dww.setStabilityFactor(1.00)
  halfNyq = .5*.5*odt
  if awarp>2:
    maxf = halfNyq*2
  else:
    maxf = halfNyq*awarp
    
  dww.setFrequencyRange(0*dt,maxf*dt)
  a = dww.getInverseAWarpNoOSamp(na,ka,st,st,awarp,hf,hg)
  h = dww.getWaveletH(na,ka,a,nh,kh)
  hsyn = getArWavelet(st,fpeak,decay,nh,kh,zp=False)
  title = "Estimated Wavelet, "+str(awarp)+\
  " time "+str(tmin)+", "+str(tmax)
  plotWavelets(Sampling(nh,st.getDelta(),kh*st.getDelta()),
    [hsyn,h],title=title)
  bdf,bbdf,bdg,bsdg,bbsdg,d = False,False,True,True,False,False
  #dww.plotDifferencePlots(bdf,bbdf,bdg,bsdg,bbsdg,d,0,3,3)
  
  bdg,bsdg = True,True
  #dww.plotDifferenceSpectrums(bdg,bsdg,0)
  #dww.plotAmplitudeSpectrum(Sampling(nh,st.getDelta(),kh*st.getDelta()),
  #hsyn,False,"WaveletSpectrum");
  print "a awarp = "+str(awarp)
  dump(a)
  return a

def estimateWavelet(awarp):
  na,ka = 3,0
  nh,kh = 250,-125
  fpeak,decay = 50,0.15
  nt,dt,ft = 500,.002,0
  st = Sampling(nt,dt,ft)
  rt = .2,.8
  ra = 1,1
  hf,hg = makeTraces(st,rt,ra,awarp,fpeak,decay,zp=None)
  plotTrace(st,hf,0,0.998,2,"hf")
  plotTrace(st,hg,0,0.998,2,"hg")

  dt = st.getDelta()
  odt = 1.0/dt
  #hf = balance(40,hf)
  #hg = balance(40,hg)
  #plotTrace(st,hf,0,0.998,2,"balance hf")
  #plotTrace(st,hg,0,0.998,2,"balance hg")
  warp = Warp()
  dww = DynamicWarpingW(warp)
  dww.setTimeRange(int(.15*odt),int(.7*odt))
  dww.setStabilityFactor(1.00)
  dww.setFrequencyRange(5*dt,85*dt)
  a = dww.getInverseAWarp(na,ka,st,st,awarp,hf,hg)
  h = dww.getWaveletH(na,ka,a,nh,kh)
  hsyn = getArWavelet(st,fpeak,decay,nh,kh)
  SimplePlot.asPoints(hsyn)
  title = "Estimated Wavelet, "+str(awarp)
  plotWavelets(Sampling(nh,st.getDelta(),kh*st.getDelta()),
    [h,hsyn],title=title)
  xcorfunc = zerofloat(nh+kh);
  xcor(nh,kh,hsyn,nh,kh,h,nh+kh,0,xcorfunc)
  bdf,bdg,bdgO,bsdgO,bsdg,d = True,True,True,True,True,True
  #dww.plotDifferencePlots(bdf,bdg,bdgO,bsdgO,bsdg,d,0,0.998,3)
  #dww.plotDifferencePlots(bdf,bdg,bsdg,d,0,0.998,3)
  bdg,bdgO,bsdgO,bsdg = True,True,True,True
  #dww.plotDifferenceSpectrums(bdg,bdgO,bsdgO,bsdg,0)
  #dww.plotDifferenceSpectrums(bdg,bsdg,0)
  print "a awarp = "+str(awarp)
  dump(a)
  return a

  
def makeTraces(st, rt, ra, awarp, fpeak, decay, zp=None):
  f,g = makeReflections(st,rt,ra,awarp)
  plotTrace(st,f,0,3,3,"f")
  plotTrace(st,g,0,3,3,"g")
  hf = addArWavelet(st,f,fpeak,decay,zp=zp)
  hg = addArWavelet(st,g,fpeak,decay,zp=zp)
  return hf,hg

def makeReflections(st, rt, ra, awarp):
  nt = st.getCount()
  dt = st.getDelta()
  ft = st.getFirst()
  nref = len(rt)
  rt1 = rt
  rt2 = mul(awarp,rt)
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

def getArWavelet(st,fpeak,decay,nh,kh,zp=False):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  print "a1 =",a1," a2 =",a2
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  if zp:
    gain *= sqrt(1.0+a1*a1+a2*a2)
  x = zerofloat(nh)
  t = zerofloat(nh)
  h = zerofloat(nh)
  x[-kh] = 1.0
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(x,t)
  if zp:
    rcf.applyReverse(t,h)
  else:
    copy(t,h)
  return h


def addArWavelet(st, p, fpeak, decay, zp=False):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  if zp:
    gain *= sqrt(1.0+a1*a1+a2*a2)
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(p,t)
  if zp:
    rcf.applyReverse(t,x)
  else:
    copy(t,x)
  return x

def plotTrace(st, p, tmin, tmax, amax, title): 
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  sp.setVLimits(tmin,tmax)
  sp.setHLimits(-amax,amax)
  sp.addTitle(title)
  pv = sp.addPoints(st,p)

def plotWavelets(st,hs,hmax=None,title=None,pngDir=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  nh = len(hs)
  hsmax = 0
  for ih in range(nh):
    if hs[ih]:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(ls[ih])
      pv.setLineWidth(2)
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
    sp.paintToPng(1000,5,pngDir+title+".png")

def balance(sigma,f):
  f = add(max(f)*0.00001,f)
  ff = mul(f,f)
  RecursiveExponentialFilter(sigma).apply1(ff,ff)
  return div(f,sqrt(ff))




#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

