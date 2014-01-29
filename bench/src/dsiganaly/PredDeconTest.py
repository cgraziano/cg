from imports import *
from dsiganaly import PredDecon
from dsiganaly import ObjectiveFunction 
from edu.mines.jtk.util.ArrayMath import *
from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.lapack import DMatrix
from edu.mines.jtk.interp.CubicInterpolator import *

from warp import WaveletNmo,WarpedWavelet,ShapingFilter
from data import SUDataGrabber
from edu.mines.jtk.awt import ColorMap;

from random import *
 
def main(args):
  nref = 1 
  na = 5
  ka = -1
  nh = 201
  kh = -50
  #goEstimateWaveletFromSynthetic(na,nref)
  #goEstimateWaveletFromSyntheticPEF(na, ka, nh, kh, nref)
  #goEstimateWaveletFromSynthetic_More3NA_PEF(na, ka, nh, kh, nref)

  #testInverseCalcSynthetic(na,0,nref)
  #makeObjectFunctPlotSynthetic(nref)

  #testInverseCalcReal(na,0)
  makeObjectFunctPlotReal()

def testInverseCalcSynthetic(na,ka,nref):
  #################### Build wavelet and CMP gather with const v##########
  st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
  nref,vnmo = nref,2.0 # number of reflectors and NMO velocity
  freq,decay = 30.0,0.1 # peak frequency and decay for wavelet
  p = makeCmpReflections(vnmo,nref,st,sx) # cmp gather without wavelet
  hp = addArWavelet(freq,decay,st,sx,p) # cmp gather with wavelet
  ak = zerofloat(na) # array for the known inverse wavelet a
  r,w = exp(-decay),2.0*PI*freq*st.delta # radius and frequency of poles
  a1,a2 = -2.0*r*cos(w),r*r # coefficients for inverse wavelet
  ak[0-ka] = 1.0
  ak[1-ka] = a1
  ak[2-ka] = a2
  print "known inverse coefficients"
  dump(ak)
  wn = WaveletNmo(st,sx,vnmo)
  fa = zerofloat(nt,nx)
  plotGather(st,sx,hp,0,2,100)
  for x in range(0,nx):
    conv(na,0,ak,nt,0,hp[x],nt,0,fa[x])
  nmofa = wn.applyNmo(fa)
  print PredDecon.semblance(nmofa)
  test = zerofloat(nt,nx)
  for i in range(0,nx):
    test[i][125] = 2.0
  plotGatherPNG(st,sx,test,"test zoom",.4,.6)
  plotGatherPNG(st,sx,nmofa,"nmofa zoom",.4,.6)
  print PredDecon.semblance(test)
  plotGather(st,sx,fa,0,2,100)
  plotGather(st,sx,nmofa,0,2,100)
  #plotGather(st,sx,fa,0,2,100)


  #################### Enhanced NMO Estimation #########################
  aenha = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  print "enhance nmo inverse coefficients"
  dump(aenha)

  pd = PredDecon(hp,na,ka)
  apredC = pd.getPredErrorCoef()
  apredD = pd.getInverseAPef(na, 0, hp) 
  print "Chris predictive decon inverse coefficients"
  dump(apredC)
  print "predictive decon inverse coefficients"
  dump(apredD)
  fa = zerofloat(nt,nx)
  for x in range(0,nx):
    conv(na,0,apred,nt,0,hp[x],nt,0,fa[x])
  #plotGather(st,sx,fa,0,2,100)
  #plotGather(st,sx,wn.applyNmo(fa),0,2,100)

  ####Use predictive deconvolution to get prediction error coefficients
  #to get the inverse of the wavelet
  #pd = PredDecon(wn.applyNmo(hp),na,ka)
  #pd = PredDecon(hp,na,ka)
  #pd = PredDecon(hp,3)
  #aprednmo = pd.getPredErrorCoef()
  #print "predictive decon nmo inverse coefficients"
  #dump(aprednmo)

def testInverseCalcReal(na,ka):
  #################### Build wavelet and CMP gather with const v##########
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first
  cmp = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=1300.strip"
  
  hp = SUDataGrabber.grab2DFile(cmp, nx, nt)
  hp = tpow(3.0,st,hp)
  tnmo = [0.010,0.646,0.846,1.150,1.511,1.825,2.490,2.985]
  vnmo = [1.510,1.575,1.687,1.817,1.938,1.980,2.446,2.735]
  vnmoInterp = interpolateVel(tnmo,st,vnmo)
  

  #################### Enhanced NMO Estimation #########################
  wn = WaveletNmo(st,sx,vnmoInterp)
  aenha = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  print "enhance nmo inverse coefficients"
  dump(aenha)

  pd = PredDecon(hp,na,ka)
  apred = pd.getPredErrorCoef()
  print "predictive decon inverse coefficients"
  dump(apred)

  ####Use predictive deconvolution to get prediction error coefficients
  #to get the inverse of the wavelet
  #pd = PredDecon(wn.applyNmo(hp),na,ka)
  #pd = PredDecon(hp,na,ka)
  #apred = pd.getPredErrorCoef()
  #print "predictive decon nmo inverse coefficients"
  dump(apred)

#Creates a synthetic trace and plots the objective function
def makeObjectFunctPlotSynthetic(nref):
  st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
  nref,vnmo = nref,2.0 # number of reflectors and NMO velocity
  freq,decay = 30.0,0.1 # peak frequency and decay for wavelet
  p = makeCmpReflections(vnmo,nref,st,sx) # cmp gather without wavelet
  hp = addArWavelet(freq,decay,st,sx,p) # cmp gather with wavelet
  na,ka = 3,0 # sampling for inverse wavelet a
  ak = zerofloat(na) # array for the known inverse wavelet a
  r,w = exp(-decay),2.0*PI*freq*st.delta # radius and frequency of poles
  a1,a2 = -2.0*r*cos(w),r*r # coefficients for inverse wavelet
  ak[0-ka] = 1.0
  ak[1-ka] = a1
  ak[2-ka] = a2
  print "known inverse coefficients"
  dump(ak)
  fa1 = a1-2.5
  fa2 = a2-2.5
  sa1 = Sampling(50,.1,fa1)
  sa2 = Sampling(50,.1,fa2)
  #ojg = PredDecon.constructPEFGatherObjFunction(hp,sa1,sa2,sx,st)
 

  #nmo corrected gather
  wn = WaveletNmo(st,sx,vnmo)
  ojgnmo = ObjectiveFunction.PEFGather(hp,sa1,sa2,sx,st,wn,vnmo)

  #weird "semblance" term
  ojgnmosem = ObjectiveFunction.SemblanceGather(hp,sa1,sa2,sx,st,wn,vnmo)
  fa = zerofloat(nt,nx)
  for x in range(0,nx):
    conv(na,0,ak,nt,0,hp[x],nt,0,fa[x])
  nmofa = wn.applyNmo(fa)
  plotGatherPNG(st,sx,hp,0,2,100,"Synthetic Gather")  
  plotGather(st,sx,fa,0,2,100)
  plotGather(st,sx,nmofa,0,2,100)
  
  print ak
  #plotObjectiveFunct(sa1,sa2,ojg,"Synth Gather objective function",ak)
  plotObjectiveFunctWa(sa1,sa2,ojgnmo,"Synth. Obj. Function PEF Term",ak)
  plotObjectiveFunctWa(sa1,sa2,ojgnmosem,"Synth. Obj. Function Semblance Term",ak)
  
def makeObjectFunctPlotReal():
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first
  cmp = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp=1301.strip"
  
  hp = SUDataGrabber.grab2DFile(cmp, nx, nt)
  hp = tpow(3.0,st,hp)
  
  #tnmo =     [0.0380286,0.475358,0.950715,1.25,1.75882,2.33876,3.54617,4.35428,5.68528]
  #vnmo =     [1.50597  ,1.52462 ,1.79203 ,1.87,2.0,2.28953,2.45122,2.65644,2.69997]
  #ts=0.010,0.646,0.846,1.150,1.511,1.825,2.490,2.985
  #vs=1.510,1.575,1.687,1.817,1.938,1.980,2.446,2.735
  ts=0.010,0.646,0.846,1.150,1.511,1.825,2.490,2.985
  vs=1.510,1.575,1.687,1.817,1.938,2.105,2.446,2.735
  vnmoInterp = interpolateVel(ts,st,vs)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addPoints(st, vnmoInterp)

  wn = WaveletNmo(st,sx,vnmoInterp)
  na = 3
  ka = 0
  pd = PredDecon(hp,na,ka)
  apred = pd.getPredErrorCoef()
  fa1 = apred[1]
  fa2 = apred[2]
  sa1 = Sampling(1,.1,fa1)
  sa2 = Sampling(1,.1,fa2)

  #nmo corrected gather
  wn = WaveletNmo(st,sx,vnmoInterp)
  ojgnmo = ObjectiveFunction.PEFGather(hp,sa1,sa2,sx,st,wn,vnmo)

  #weird "semblance" term
  ojgnmosem = ObjectiveFunction.SemblanceGather(hp,sa1,sa2,sx,st,wn,vnmo)
  
  
  plotObjectiveFunctWa(sa1,sa2,ojgnmo,"VG Obj. Function PEF Term",apred)
  plotObjectiveFunctWa(sa1,sa2,ojgnmosem,"VG Obj. Function Semblance Term",apred)
  fa = zerofloat(nt,nx)
  for x in range(0,nx):
    conv(na,0,apred,nt,0,hp[x],nt,0,fa[x])
  nmofa = wn.applyNmo(fa)
  plotGatherPNG(st,sx,hp,0,6,100,"Viking Graben CDP 1300")  
  plotGather(st,sx,fa,0,6,100)  

  plotGather(st,sx,nmofa,0,6,100)


def goEstimateWaveletFromSynthetic(na,nref):
  #Input data
  st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
  nref,vnmo = nref,2.0 # number of reflectors and NMO velocity
  freq,decay = 30.0,0.1 # peak frequency and decay for wavelet
  p = makeCmpReflections(vnmo,1,st,sx) # cmp gather without wavelet
  hp = addArWavelet(freq,decay,st,sx,p) # cmp gather with wavelet
  plotGather(st,sx,hp,0,2,100,"CMP gather")
  na,ka = na,0 # sampling for inverse wavelet a
  nh,kh = 201,-50 # sampling for wavelet h
  ak = zerofloat(na) # array for the known inverse wavelet a
  r,w = exp(-decay),2.0*PI*freq*st.delta # radius and frequency of poles
  a1,a2 = -2.0*r*cos(w),r*r # coefficients for inverse wavelet
  ak = zerofloat(na) # array for the known inverse wavelet a
  ak[0-ka] = 1.0
  ak[1-ka] = a1
  ak[2-ka] = a2

  wn = WaveletNmo(st,sx,vnmo)
  a = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
  g = wn.applyHNmoA(na,ka,a,nh,kh,h,hp) #convolve inverse, apply NMO, convolve wavelet
  #apply NMO to original data
  e = wn.applyNmo(hp)
  
  hk = wn.getWaveletH(na,ka,ak,nh,kh); # estimate wavelet
  
  tmin,tmax,perc = 0.0,2,100.0
  #plotGather(st,sx,d,tmin=tmin,tmax=tmax,perc=perc,title="Difference")
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  #plotGather(st,sx,hp,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  sa = Sampling(na,st.delta,ka*st.delta)
  sh = Sampling(nh,st.delta,kh*st.delta)
  plot2Sequences(sa,sa,ak,a,title="inverse")
  plot2Sequences(sh,sh,normalize(hk),normalize(h),title="estimated wavelet")
 
def goEstimateWaveletFromSyntheticPEF(na, ka, nh, kh, nref):
  #Input data
  st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
  nref,vnmo = nref,2.0 # number of reflectors and NMO velocity
  freq,decay = 30.0,0.1 # peak frequency and decay for wavelet
  p = makeCmpReflections(vnmo,1,st,sx) # cmp gather without wavelet
  hp = addArWavelet(freq,decay,st,sx,p) # cmp gather with wavelet
  plotGather(st,sx,hp,0,2,100,"CMP gather")
  na,ka = na,ka # sampling for inverse wavelet a
  nh,kh = nh,kh # sampling for wavelet h
  ak = zerofloat(na) # array for the known inverse wavelet a
  r,w = exp(-decay),2.0*PI*freq*st.delta # radius and frequency of poles
  a1,a2 = -2.0*r*cos(w),r*r # coefficients for inverse wavelet
  ak = zerofloat(na) # array for the known inverse wavelet a
  ak[0-ka] = 1.0
  ak[1-ka] = a1
  ak[2-ka] = a2

  wn = WaveletNmo(st,sx,vnmo)
  a = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
  g = wn.applyHNmoA(na,ka,a,nh,kh,h,hp) #convolve inverse, apply NMO, convolve wavelet
  #apply NMO to original data
  e = wn.applyNmo(hp)
  
  hk = wn.getWaveletH(na,ka,ak,nh,kh); # estimate wavelet

  pd = PredDecon(hp,na,ka)
  apred = pd.getPredErrorCoef()
  hpred = wn.getWaveletH(na,ka,apred,nh,kh); # estimate wavelet
  
  
  tmin,tmax,perc = 0.0,2,100.0
  #plotGather(st,sx,d,tmin=tmin,tmax=tmax,perc=perc,title="Difference")
  #plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  #plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  #plotGather(st,sx,hp,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  sa = Sampling(na,st.delta,ka*st.delta)
  sh = Sampling(nh,st.delta,kh*st.delta)
  plot3Sequences(sa,sa,sa,ak,a,apred,title="inverse")
  plot3Sequences(sh,sh,sh,normalize(hk),normalize(h),normalize(hpred),title="estimated wavelet")
  
def goEstimateWaveletFromSynthetic_More3NA_PEF(na, ka, nh, kh, nref):
  #Input data
  st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
  nref,vnmo = nref,2.0 # number of reflectors and NMO velocity
  freq,decay = 30.0,0.1 # peak frequency and decay for wavelet
  p = makeCmpReflections(vnmo,1,st,sx) # cmp gather without wavelet
  hp = addArWavelet(freq,decay,st,sx,p) # cmp gather with wavelet
  plotGather(st,sx,hp,0,2,100,"CMP gather")
  na,ka = na,ka # sampling for inverse wavelet a
  nh,kh = nh,kh # sampling for wavelet h

  wn = WaveletNmo(st,sx,vnmo)
  a = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
  g = wn.applyHNmoA(na,ka,a,nh,kh,h,hp) #convolve inverse, apply NMO, convolve wavelet

  pd = PredDecon(hp,na,ka)
  apred = pd.getPredErrorCoef()
  hpred = wn.getWaveletH(na,ka,apred,nh,kh); # estimate wavelet
  
  
  tmin,tmax,perc = 0.0,2,100.0
  #plotGather(st,sx,d,tmin=tmin,tmax=tmax,perc=perc,title="Difference")
  #plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  #plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  #plotGather(st,sx,hp,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  sa = Sampling(na,st.delta,ka*st.delta)
  sh = Sampling(nh,st.delta,kh*st.delta)
  plot2Sequences(sa,sa,a,apred,title="inverse")
  plot2Sequences(sh,sh,normalize(h),normalize(hpred),title="estimated wavelet")
  

def makeCmpReflections(vel,nref,st,sx):
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  p = zerofloat(nt,nx)
  ts = add(ft,mul((nt-1)*dt,randfloat(nref)))
  ts = [0.5]
  rs = sub(mul(2.0,randfloat(nref)),1.0)
  rs = [2.0]
  si = SincInterp.fromErrorAndFrequency(0.01,0.45)
  for jx in range(nx):
    xj = sx.getValue(jx)
    cj = (xj*xj)/(vel*vel)
    for jr in range(nref):
      tj = ts[jr]
      rj = rs[jr]
      tj = sqrt(tj*tj+cj)
      si.accumulate(tj,rj,nt,dt,ft,p[jx])
  return p

def addArWavelet(fpeak,decay,st,sx,p):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  #print "a1 =",a1," a2 =",a2
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  x = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.apply1Forward(p,x)
  return x

def plotObjectiveFunctWa(sa1,sa2,objf,title,a):
  pv = PixelsView(sa1,sa2,objf)
  pv.setColorModel(ColorMap.JET)#PRISM)#HUE_BLUE_TO_RED)
  pp = PlotPanel()
  pov = PointsView([a[1]],[a[2]])
  pov.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  pov.setMarkSize(10)
  pov.setMarkColor(Color.GREEN)
  pp.addTiledView(pv)
  pp.addTiledView(pov)
  pp.setHLabel("a1")
  pp.setVLabel("a2")
  pp.setTitle(title)
  pp.addColorBar()
  pf = PlotFrame(pp)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pf.setFontSize(18)
  pf.paintToPng(100,5,title+".png")

def tpow(power,st,f):
  """Applies t^power gain."""
  nt,dt,ft = st.count,st.delta,st.first
  nx = len(f)
  tp = pow(rampfloat(ft,dt,nt),power) # sampled times raised to power
  g = zerofloat(nt,nx)
  for ix in range(nx):
    mul(tp,f[ix],g[ix])
  return g

#take tnmo and vnmo and interpolates vnmo^2 and turns it back
#to vnmo
def interpolateVel(tnmo, st, vnmoRaw):
  nt = st.getCount();
  dt = st.getDelta();
  t = []
  for i in range(0,nt):
    t.append(i*dt)
  #Build squared nmo values
  nvnmo = len(vnmoRaw)
  vnmoRawSq = []
  for i in range(0,nvnmo):
    vnmoRawSq.append(vnmoRaw[i]*vnmoRaw[i])

  ci = CubicInterpolator(Method.MONOTONIC, tnmo, vnmoRawSq)
  vnmoInterpSq = ci.interpolate(t)

  #Square Root vnmo squared values
  vnmoInterp = []
  for i in range(0,nt):
    vnmoInterp.append(sqrt(vnmoInterpSq[i]))
  return vnmoInterp


def plotGather(st,sx,p,tmin=None,tmax=None,perc=None,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.setHLabel("Offset (km)")
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  sp.setVLimits(tmin,tmax)
  pv = sp.addPixels(st,sx,p)
  if perc:
    pv.setPercentiles(100-perc,perc)
  sp.getPlotPanel()

def plotSequence(s,x,xmax=None,title=None):
  sp = SimplePlot.asPoints(s,x)
  if xmax:
    sp.setVLimits(-xmax,xmax)
  if title:
    sp.setTitle(title)

def plot2Sequences(sr,sg,r,g,xmax=None,title=None):
  pvr = PointsView(sr,r)
  pvr.setLineColor(Color.RED)
  pvg = PointsView(sg,g)
  pvg.setLineColor(Color.GREEN)

  pp = PlotPanel()
  pp.addTiledView(pvr)
  pp.addTiledView(pvg)
  if xmax:
    pp.setVLimits(-xmax,xmax)
  if title:
    pp.setTitle(title)

  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pf.setVisible(True)

def plot3Sequences(sr,sg,sb,r,g,b,xmax=None,title=None):
  pvr = PointsView(sr,r)
  pvr.setLineColor(Color.RED)
  pvg = PointsView(sg,g)
  pvg.setLineColor(Color.GREEN)
  pvb = PointsView(sb,b)
  pvb.setLineColor(Color.BLUE)

  pp = PlotPanel()
  pp.addTiledView(pvr)
  pp.addTiledView(pvg)
  pp.addTiledView(pvb)
  if xmax:
    pp.setVLimits(-xmax,xmax)
  if title:
    pp.setTitle(title)

  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pf.setVisible(True)

def normalize(h):
  return div(h,max(h))
  
def plotGatherPNG(st,sx,p,tmin=None,tmax=None,perc=None,title=None):
  pv = PixelsView(st,sx,p)
  pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.addTiledView(pv)
  pp.setTitle(title)
  if tmin:
    pp.setVLimits(tmin,tmax)
  pp.setHLabel("Offset (km)")
  pp.setVLabel("Time (s)")
  pp.addColorBar()
  pf = PlotFrame(pp)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pf.setFontSize(18)
  pf.paintToPng(100,5,title+".png")
  pp = PlotPanel()


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
