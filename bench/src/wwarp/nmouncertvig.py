#############################################################################
# Demo wavelet estimation from NMO stretch
# This script give more or less uncertainty for different parts of
# the Viking Graben data.
#############################################################################



from imports import *

from edu.mines.jtk.dsp.Conv import *
from warp import NormalMoveout,WaveletNmo

#############################################################################

pngDir = "./png/nmouncertvig/"
#pngDir = None

def main(args):
  tmin = 0
  tmax = 3
  funcert = .5
  luncert = .01
  goEstimateWaveletFromGather(tmin,tmax,funcert,luncert)

def goEstimateWaveletFromGather(tmin,tmax,funcert,luncert):
  """ Estimates wavelet from a gather sampled in time and offset """
  st,sx,f = getVikingGrabenCDP(800)
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  ts=0.010,0.646,0.846,1.150,1.511,1.825,2.490,2.985
  vs=1.510,1.575,1.687,1.817,1.938,1.980,2.446,2.735
  vnmoP = CubicInterpolator(ts,vs).interpolate(rampfloat(ft,dt,nt))
  vnmoM = fillfloat(1.500,nt) # water velocity, for multiples
  vnmo = vnmoP
  na,ka = 16,0 # sampling for inverse wavelet a
  nh,kh = 151,-25 # sampling for wavelet h
  fmin,fmax,sfac = 5.0,75.0,1.00
  texp,tbal,smax = 2.0,100,9.00
  perc = 98.0
  #f = tpow(texp,st,f)
  if tbal>0:
    f = balance(tbal,f)
  #for whole gather
  #tmin = 0.00
  #tmax = 2.00
  #inverseAndWaveletCalc(sx,st,f,vnmo,fmin,fmax,sfac,
  #    na,ka,nh,kh,hsyn,tmin,tmax)
  itmin = int(tmin/dt)
  itmax = int(tmax/dt)
  nw = itmax-itmin+1
  w = zerofloat(nw,3)
  w[0] = fillfloat(1.0,nw)
  w[1] = rampWeights(funcert,luncert,nw)
  w[2] = rampWeights(luncert,funcert,nw)
  uncerts = ["Constant Uncertainty","Ramp Up Uncertainty","Ramp Down Uncertainty"]
  hpefs = zerofloat(nh,3)
  hnmos = zerofloat(nh,3)
  for u in range(0,3):
    hpef,hnmo = inverseAndWaveletCalc(sx,st,f,vnmo,smax,fmin,fmax,sfac,
    na,ka,nh,kh,tmin,tmax,perc,uncerts[u],w=w[u])
    #hpef,hnmo = inverseAndWaveletCalc(sx,st,f,vnmo,fmin,fmax,sfac,
    #na,ka,nh,kh,tmin,tmax)
    hnmos[u] = hnmo
    hpefs[u] = hpef
  title = "Constant Uncert., Uncert. "+str(funcert)+" "+str(luncert)+", Uncert. "+str(luncert)+" "+str(funcert)
  titlePEF = "PEF "+title
  titleNMO = "NMO "+title

  hmax = max(hpefs)
  if hmax < max(hnmos):
    hmax = max(hnmos)

  plotTimeWindowWavelets(Sampling(nh,st.delta,kh*st.delta),hpefs,hmax,
               title=titlePEF,pngDir=pngDir)
  plotTimeWindowWavelets(Sampling(nh,st.delta,kh*st.delta),hnmos,hmax,
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
def plotTimeWindowWavelets(st,hs,hmax,title=None,pngDir=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  lc = [Color.RED,Color.GREEN,Color.BLUE]
  nh = len(hs)
  hsmax = 0
  c = 0
  s = 0
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

def inverseAndWaveletCalc(sx,st,f,vnmo,smax,fmin,fmax,sfac,
    na,ka,nh,kh,tmin,tmax,perc,uncert,w=None):
  tmino = 0
  tmaxo = 3
  hmax = 5  
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  nmo = NormalMoveout()
  nmo.setStretchMax(smax)
  if w:
    wn = WaveletNmo(nmo,w)
  else:
    wn = WaveletNmo(nmo)
  wn.setTimeRange(int(tmin/dt),int(tmax/dt))
  wn.setFrequencyRange(fmin*dt,fmax*dt)
  wn.setStabilityFactor(sfac)
  e = nmo.apply(st,sx,vnmo,f)
  apef = wn.getInverseAPef(na,ka,f)
  anmo = wn.getInverseANmo(na,ka,st,sx,vnmo,f)
  hpef = wn.getWaveletH(na,ka,apef,nh,kh);
  hnmo = wn.getWaveletH(na,ka,anmo,nh,kh);
  title = "Estimated Wavelets "+str(uncert)+" "+str(tmin)+" s to "+str(tmax)+" s"
  plotWavelets(Sampling(nh,st.delta,kh*st.delta),[hnmo,hpef],hmax,
               title=title,pngDir=pngDir)
  title = "Input Gather "+str(uncert)+" "+str(tmin)+" s to "+str(tmax)+" s"
  plotGather(st,sx,f,tmin=tmin,tmax=tmax,title=title,perc=perc,pngDir=pngDir)
  title = "Conventional NMO "+str(uncert)+" "+str(tmin)+" s to "+str(tmax)+" s"
  plotGather(st,sx,e,tmin=tmino,tmax=tmaxo,perc=perc,title=title,pngDir=pngDir)
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
    print tlist[ia]+" epef =",epef," enmo =",enmo," enor =",enor
    print " a ="; dump(a)
    title = t+" improved NMO "+str(uncert)+" "+str(tmin)+" s to "+str(tmax)+" s"
    plotGather(st,sx,g,tmin=tmino,tmax=tmaxo,perc=perc,
        title=title,pngDir=pngDir)
    #plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,
    #  title=t+": stack error")
    #plotSequence(Sampling(na,st.delta,ka*st.delta),normalize(a),
    #             title="inverse wavelet")
    #plotSequence(Sampling(nah,st.delta,kah*st.delta),normalize(ah),
    #             title="unit impulse")
  #d = wn.getDifferenceGathers(na,ka,st,sx,vnmo,f)
  #for ia in range(na):
  #  plotGather(st,sx,d[ia],
  #             tmin=tmin,tmax=tmax,perc=perc,title="lag="+str(ka+ia))

  return hpef,hnmo

def rampWeights(fs, ls, n1):
  m = (ls-fs)/(n1-1)
  s = rampfloat(fs,m,n1)
  ss = mul(s,s)
  w = div(1.0,ss)
  return w

def getVikingGrabenCDP(icdp):
  if 145 <= icdp | icdp <= 147:
    fileName = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp145_147.dat"
    return getCDP(icdp, 145, fileName)
  elif 400 <= icdp | icdp <= 800:
    fileName = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp400_800.dat"
    return getCDP(icdp, 400, fileName)
  elif 1100 <= icdp | icdp <= 1500:
    fileName = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp1100_1500.dat"
    return getCDP(icdp, 1100, fileName)
  elif 1700 <= icdp | icdp <= 2000:
    fileName = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/cdp1700_2000.dat"
    return getCDP(icdp, 1700, fileName)
  else:
    print "Requested CDP not available"
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
#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())


