#############################################################################
# Demo wavelet estimation from NMO stretch
# This script estimates multiple wavelets at different times for one defined CMP 
# gather in synthetic data. 
#############################################################################


from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import NormalMoveout,WaveletNmo

#############################################################################

#pngDir = "./png/nmotime/"
pngDir = None

def main(args):
  timewindow = [.5,1.75]
  #timewindow = [.5,1.75,1.0,1.75]
  #timewindow = [0,.5,.5,1.0,1.0,1.5,1.5,2.0]
  funcert = 1.0
  luncert = 1.0
  goEstimateWaveletFromGather(args[1],timewindow,funcert,luncert)

def goEstimateWaveletFromGather(name,timewindow,funcert,luncert):
  """ Estimates wavelet from a gather sampled in time and offset """
  print name
  if name == "syn1": # Synthetic with one reflector
    st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
    nref,vnmo = 1,2.0 # number of reflectors and NMO velocity
    tran,tbed = False,False
    freq,decay = 20.0,0.05 # peak frequency and decay for wavelet
    fmin,fmax,sfac = 0.0,50.0,1.00
    texp,tbal = 0.00,0
    checkTimeWindow(st, timewindow)
    
    zp = True# zero-phase?
    if zp:
      na,ka = 5,-2
      nh,kh = 251,-125
      decay *= 4.0
    else:
      na,ka = 11,0
      nh,kh = 151,-25
    hsyn = getArWavelet(freq,decay,st,nh,kh,zp)
  elif name == "synr": # Synthetic with random reflectors
    st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
    tran,tbed = True,False
    nref,vnmo = 40,2.0 # number of reflectors and NMO velocity
    freq,decay = 20.0,0.05 # peak frequency and decay for wavelet
    fmin,fmax,sfac = 0.0,50.0,1.00
    texp,tbal = 0.00,0
    checkTimeWindow(st, timewindow)
    zp = True # zero-phase?
    na,ka = 5,-2 # sampling for inverse wavelet a
    nh,kh = 251,-125 # sampling for wavelet h
    hsyn = getArWavelet(freq,decay,st,nh,kh,zp)
  elif name == "synt": # Synthetic with random thin beds
    st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
    sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
    tran,tbed = True,True
    nref,vnmo = 40,2.0 # number of reflectors and NMO velocity
    freq,decay = 20.0,0.05 # peak frequency and decay for wavelet
    fmin,fmax,sfac = 0.0,50.0,1.0001
    texp,tbal = 0.00,0
    checkTimeWindow(st, timewindow)
    zp = True # zero-phase?
    if zp:
      na,ka = 5,-2
      nh,kh = 251,-125
      decay *= 4.0
    else:
      na,ka = 11,0
      nh,kh = 151,-25
    hsyn = getArWavelet(freq,decay,st,nh,kh,zp)

  f = zerofloat(nt,nx)
  p = makeCmpReflections(vnmo,nref,st,sx,random=tran,thinBeds=tbed)
  print "zero phase = ",zp
  f = addArWavelet(freq,decay,st,sx,p,zp)
  f = tpow(texp,st,f)
  if tbal>0:
    f = balance(tbal,f)

  ntw = len(timewindow)
  #for whole gather
  #tmin = 0.00
  #tmax = 2.00
  #itmin = int(tmin/dt)
  #itmax = int(tmax/dt)
  #nw = itmax-itmin+1
  #w = rampWeights(funcert,luncert,nw)
  #inverseAndWaveletCalc(sx,st,f,vnmo,fmin,fmax,sfac,
  #    na,ka,nh,kh,hsyn,tmin,tmax,w=w)
  hpefs = zerofloat(nh,ntw/2)
  hnmos = zerofloat(nh,ntw/2)
  t = 0
  while t<ntw:
    tmin = timewindow[t]
    tmax = timewindow[t+1]
    itmin = int(tmin/dt)
    itmax = int(tmax/dt)
    nw = itmax-itmin+1
    w = rampWeights(funcert,luncert,nw)
    hpef,hnmo = inverseAndWaveletCalc(sx,st,f,vnmo,fmin,fmax,sfac,
    na,ka,nh,kh,hsyn,tmin,tmax,w=w)
    #hpef,hnmo = inverseAndWaveletCalc(sx,st,f,vnmo,fmin,fmax,sfac,
    #na,ka,nh,kh,hsyn,tmin,tmax)
    hpefs[t/2] = hpef
    hnmos[t/2] = hnmo
    t+=2
  title = "Time Window Wavelets "
  t=0
  while t<ntw:
    title = title+" "+str(timewindow[t])+"-"+str(timewindow[t+1])+","
    t += 2
  titlePEF = "PEF "+title
  titleNMO = "NMO "+title

  hmax = max(hpefs)
  if hmax < max(hnmos):
    hmax = max(hnmos)

  plotTimeWindowWavelets(Sampling(nh,st.delta,kh*st.delta),hpefs,hsyn,hmax,
               title=titlePEF,pngDir=pngDir)
  plotTimeWindowWavelets(Sampling(nh,st.delta,kh*st.delta),hnmos,hsyn,hmax,
               title=titleNMO,pngDir=pngDir)

def stackError(g):
  nt,nx = len(g[0]),len(g)
  s = copy(g[0])
  for ix in range(1,nx):
    s = add(s,g[ix])
  s = div(s,nx)
  g = copy(g)
  for ix in range(nx):
    g[ix] = sub(g[ix],s)
  return g

def rms(g):
  nt,nx = len(g[0]),len(g)
  return sqrt(sum(mul(g,g))/nt/nx)

def tpow(power,st,f):
  """Applies t^power gain."""
  nt,dt,ft = st.count,st.delta,st.first
  nx = len(f)
  tp = pow(rampfloat(ft,dt,nt),power) # sampled times raised to power
  g = zerofloat(nt,nx)
  for ix in range(nx):
    mul(tp,f[ix],g[ix])
  return g

def balance(sigma,f):
  f = add(max(f)*0.00001,f)
  ff = mul(f,f)
  RecursiveExponentialFilter(sigma).apply1(ff,ff)
  return div(f,sqrt(ff))

def normalize(h):
  return div(h,max(max(h),-min(h)))

def makeCmpReflections(vel,nref,st,sx,random=False,thinBeds=False):
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  p = zerofloat(nt,nx)
  if random:
    rand = Random(21) #14, 3141, 6, 17
    ts = add(ft,mul((nt-1)*dt,randfloat(rand,nref)))
    rs = sub(mul(2.0,randfloat(rand,nref)),1.0)
  else:
    ts = rampfloat(nt*dt/(nref+1),nt*dt/(nref+1),nref)
    rs = fillfloat(1.0,nref)
  if thinBeds:
    tsc = copy(ts)
    rsc = copy(rs)
    nref *= 2
    ts = zerofloat(nref)
    rs = zerofloat(nref)
    for iref in range(0,nref,2):
      ts[iref] = tsc[iref/2]
      ts[iref+1] = ts[iref]+5.0*dt
      rs[iref] = rsc[iref/2]
      rs[iref+1] = -rs[iref]
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

def getArWavelet(fpeak,decay,st,nh,kh,zp=False):
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

def addArWavelet(fpeak,decay,st,sx,p,zp=False):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  print "a1 =",a1," a2 =",a2
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  if zp:
    gain *= sqrt(1.0+a1*a1+a2*a2)
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.apply1Forward(p,t)
  if zp:
    rcf.apply1Reverse(t,x)
  else:
    copy(t,x)
  return x

def plotGather(st,sx,p,tmin=None,tmax=None,perc=None,
    title=None,pngDir=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    if pngDir==None:
      sp.setTitle(title)
  sp.setHLabel("Offset (km)")
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  if tmin!=None and tmax!=None:
    sp.setVLimits(tmin,tmax)
  pv = sp.addPixels(st,sx,p)
  if perc:
    pv.setPercentiles(100-perc,perc)
  if pngDir:
    sp.paintToPng(1000,5,pngDir+title+".png")

def plotSequence(s,x,xmax=None,title=None,pngDir=None):
  sp = SimplePlot.asPoints(s,x)
  if xmax==None:
    xmax = max(abs(max(x)),abs(min(x)))
    xmax *= 1.05
  sp.setVLimits(-xmax,xmax)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Amplitude")
  if title:
    if pngDir==None:
      sp.setTitle(title)
  if pngDir:
    sp.paintToPng(1000,5,pngDir+title+".png")

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


def plotTimeWindowWavelets(st,hs,hsyn,hmax,title=None,pngDir=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  lc = [Color.RED,Color.GREEN,Color.BLUE]
  nh = len(hs)
  hsmax = 0
  c = 0
  s = 0
  pv = sp.addPoints(st,hsyn)
  for ih in range(nh):
    if hs[ih]:
      pv = sp.addPoints(st,hs[ih])
      if s > 2:
        c += 1
        s = 0
      pv.setLineStyle(ls[s])
      pv.setLineColor(lc[c])
      pv.setLineWidth(2)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
      s += 1
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

def checkTimeWindow(st, timewindow):
  if acceptableTimeWindow(st, timewindow) == False:
      print "TIME WINDOWS NOT ACCEPTABLE!"
  else:
    print "Time Windows:"
    ntw = len(timewindow)
    t = 0
    while t<ntw:
      print str(timewindow[t]) + " " + str(timewindow[t+1])
      t += 2

def acceptableTimeWindow(st, timewindow):
  st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  ntw = len(timewindow)
  if ntw==0: 
    print "Need time window(s)"
    return False
  elif ntw%2==1:
    print "Need pair of times for time window"
    return False
  elif ft > timewindow[0]:
    print "1st time point before gather starts"
    return False
  elif nt*dt < timewindow[ntw-1]:
    print "last time point after gather ends"
    return False
  else:
    return True

def inverseAndWaveletCalc(sx,st,f,vnmo,fmin,fmax,sfac,
    na,ka,nh,kh,hsyn,tmin,tmax,w=None):
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  vnmo = fillfloat(vnmo,nt)
  nmo = NormalMoveout()
  #nmo.setStretchMax(9.0)
  if w:
    wn = WaveletNmo(nmo,w)
  else:
    wn = WaveletNmo(nmo)
  wn.setFrequencyRange(fmin*dt,fmax*dt)
  wn.setTimeRange(int(tmin/dt),int(tmax/dt))
  wn.setStabilityFactor(sfac)
  e = nmo.apply(st,sx,vnmo,f)
  apef = wn.getInverseAPef(na,ka,f)
  anmo = wn.getInverseANmo(na,ka,st,sx,vnmo,f)
  hpef = wn.getWaveletH(na,ka,apef,nh,kh)
  hnmo = wn.getWaveletH(na,ka,anmo,nh,kh)
  title = "Estimated Wavelets "+str(tmin)+" s to "+str(tmax)+" s"
  plotWavelets(Sampling(nh,st.delta,kh*st.delta),[hnmo,hpef,hsyn],
               title=title,pngDir=pngDir)
  title = "Input Gather "+str(tmin)+" s to "+str(tmax)+" s"
  plotGather(st,sx,f,tmin=tmin,tmax=tmax,title=title,pngDir=pngDir)
  title = "Conventional NMO "+str(tmin)+" s to "+str(tmax)+" s"
  plotGather(st,sx,e,title=title,pngDir=pngDir)
  alist = [apef,anmo]
  hlist = [hpef,hnmo]
  tlist = ["PEF","NMO"]
  for ia in range(0,len(alist)):
    a = alist[ia]
    h = hlist[ia]
    t = tlist[ia]
    nah = na+nh
    kah = ka+kh
    ah = zerofloat(nah)
    conv(na,ka,a,nh,kh,h,nah,kah,ah)
    g = wn.applyHNmoA(na,ka,a,nh,kh,h,st,sx,vnmo,f)
    epef = wn.getVariancePef(na,ka,a,f)
    enmo = wn.getVarianceNmo(na,ka,a,st,sx,vnmo,f)
    enor = wn.getNormalizedVarianceNmo(na,ka,a,st,sx,vnmo,f)
    print t+": epef =",epef," enmo =",enmo," enor =",enor
    print " a ="; dump(a)
    title = t+" improved NMO "+str(tmin)+" s to "+str(tmax)+" s"
    plotGather(st,sx,g,title=title,pngDir=pngDir)
  return hpef,hnmo

def rampWeights(fs, ls, n1):
  m = (ls-fs)/(n1-1)
  s = rampfloat(fs,m,n1)
  ss = mul(s,s)
  w = div(1.0,ss)
  return w


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

