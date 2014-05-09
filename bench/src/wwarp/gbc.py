from imports import *
from edu.mines.jtk.dsp.Conv import *
from dwarp import DynamicWarpingW 

datadir = "C:/Users/Chris/Documents/CWP/Research/research/gbc/dat"
pngDir = "./png/balancetest/"
#pngDir = None
#Tests of the warping with wavelets algorithm with another dataset.
###########################Read#########################
#The script that should be used to try warping with wavelets on your own 
#2D data set is waveletwarpingha.py
#######################################################

def main(args):
  #Available Options: pp, ps1, ps2, shiftl, shiftc
  pp = "pp"
  ps = "ps1"
  shift = "shiftl"
  iline = 70 
  ka = 0 
  na = 50
  kh = 0
  nh = 251

  pp = getData(pp,iline)
  ps = getData(ps,iline)
  shift = getData(shift,iline)
  stpp = Sampling(2000,.002,0.0)
  stps = Sampling(2000,.002,0.0)
  sx = Sampling(145,0.033529,0.0)
  ###################################d#
  #balanceCheckS(iline, sx, stpp, stps, pp, ps, shift, 80, 80)
  testInverse(ka, na, kh, nh, stpp, stps, sx, pp, ps, shift)
  

def testInverse(ka, na, kh, nh, stpp, stps, sx, pp, ps, shift):
  fmin = 5
  fmax = 85
  tppmin = .5
  tppmax = 2.5
  tpsmin = .7
  tpsmax = 3.5
  amax = 6
  dt = stpp.delta
  dww = DynamicWarpingW()
  dww.setStabilityFactor(1.01)
  dww.setFrequencyRange(5*dt, 85*dt)
  dww.setTimeRange(300,1150)

  pp = balance(80,pp)
  ps = balance(80,ps)
  e = dww.applyShifts(stps,ps,shift)
  adww = dww.getInverseAWarp(na,ka,stpp,stps,sx,shift,pp,ps)
  apef = dww.getInverseAPef(na,ka,ps)

  plotDifferenceGathers(dww,na,ka,stpp,stps,sx,shift,pp,ps)

  #plotSlice(stps,sx,ps,tmin=0.0,tmax=3.998,title="ps")
  hdww = dww.getWaveletH(na, ka, adww, nh, kh)
  hpef = dww.getWaveletH(na, ka, apef, nh, kh)
  title = "Estimated Wavelets"
  plotWavelets(Sampling(nh,stpp.delta,kh*stpp.delta),[hdww,hpef],
               title=title)
  title = "PS Data"
  plotSlice(stps,sx,ps,tmin=tpsmin,tmax=tpsmax,amax=amax,title=title)
  title = "Conventional Warping"
  plotSlice(stpp,sx,e,tmin=tppmin,tmax=tppmax,amax=amax,title=title)
  #alist = [adww]
  #hlist = [hdww]
  #tlist = ["DWW"]
  alist = [apef,adww]
  hlist = [hpef,hdww]
  tlist = ["PEF","DWW"]
  for ia in range(0,len(alist)):
    a = alist[ia]
    h = hlist[ia]
    t = tlist[ia]
    nah = na+nh
    kah = ka+kh
    ah = zerofloat(nah)
    conv(na,ka,a,nh,kh,h,nah,kah,ah)
    g = dww.applyHSA(na,ka,a,nh,kh,h,stpp,stps,sx,shift,pp,ps)
    epef = dww.getVariancePef(na,ka,a,ps)
    edww = dww.getVarianceDww(na,ka,a,stpp,stps,sx,shift,pp,ps)
    #enor = dww.getNormalizedVarianceNmo(na,ka,a,st,sx,vnmo,f)
    print t+": epef =",epef," edww =",edww#," enor =",enor
    print " a ="; dump(a)
    title = t+" improved Warping"
    plotSlice(stpp,sx,g,tmin=tppmin,tmax=tppmax,amax=amax,title=title)


def balanceCheckS(iline, sx, stpp, stps, pp, ps, shift, sigpp, sigwarp):
  print "Available Options: pp, ps1, ps2, shiftl, shiftc"
  amax = 6
  dww = DynamicWarpingW()
  print "sigpp = "+str(sigpp)+" sigps = "+str(sigwarp)

  print "rms pp = "+str(dww.rms(stpp,sx,0.5,2.4,pp))
  plotSlice(stpp,sx,pp,tmin=0.5,tmax=2.4,amax=amax,title="pp")
  if sigpp!=0:
    pp = balance(sigpp,pp)

  print "balanced pp rms = "+str(dww.rms(stpp,sx,0.5,2.4,pp))
  plotSlice(stpp,sx,pp,tmin=0.5,tmax=2.4,amax=amax,title="balanced pp")

  print "ps rms = "+str(dww.rms(stps,sx,0.5,3.4,ps))
  plotSlice(stps,sx,ps,tmin=0.5,tmax=3.4,amax=amax,title="ps")

  #balance then warp
  if sigwarp!=0:
    ps = balance(sigwarp,ps)

  print "balanced ps rms = "+str(dww.rms(stps,sx,0.5,3.4,ps))
  plotSlice(stps,sx,ps,tmin=0.5,tmax=3.4,amax=amax,title="balanced ps")

  warp = dww.applyShifts(stps,ps,shift)
  print "warp balanced ps rms = "+str(dww.rms(stpp,sx,0.5,2.4,warp))
  plotSlice(stpp,sx,warp,tmin=0.5,tmax=2.4,amax=amax,\
      title="warp balanced ps")
  
  #warp then balance
  """
  warp = dww.applyShifts(stps,ps,shift)
  print "warp ps rms = "+str(dww.rms(stpp,sx,0.5,2.4,warp))
  plotSlice(stps,sx,warp,tmin=0.5,tmax=2.4,amax=amax,title="warp ps")

  if sigwarp!=0:
    warp = balance(sigwarp,warp)

  print "balanced warp ps rms = "+str(dww.rms(stps,sx,0.5,3.4,warp))
  plotSlice(stps,sx,warp,tmin=0.5,tmax=2.4,amax=amax,\
      title="balanced warp ps")
  """
  d = sub(warp,pp)
  plotSlice(stpp,sx,d,tmin=0.5,tmax=2.4,amax=amax,title="d")
  print "rms d = "+str(dww.rms(stpp,sx,0.5,2.4,d))

def plotDifferenceGathers(dww, na, ka, stf, stg, sx, shifts, f, g):
  d = dww.getDifferenceGathers(na,ka,stf,stg,sx,shifts,f,g)
  lag = ka
  for ia in range(0,na):
    title = "d lag = "+str(lag)
    plotSlice(stf,sx,d[ia],tmin=0.5,tmax=2.4,amax=6,title=title)
    lag += 1




def getData(data,il=None):
  if data == "pp":
    f = readImage(datadir+"/pp.dat",2000,150,145)
  elif data == "ps1":
    f = readImage(datadir+"/ps1.dat",2000,150,145)
  elif data == "ps2":
    f = readImage(datadir+"/ps2.dat",2000,150,145)
  elif data == "shiftl":
    f = readImage(datadir+"/u_sag_linear.dat",2000,150,145)
  elif data == "shiftc":
    f = readImage(datadir+"/u_sag_monotonic.dat",2000,150,145)

  if il:
    inline = zerofloat(2000,145)
    for i3 in range(0,145):
      for i1 in range(0,2000):
        inline[i3][i1] = f[il][i3][i1]
    #sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    #sp.addPixels(inline)
    
    
    return inline 
  else:
    #sf = SimpleFrame.asImagePanels(f)
    return f

def readImage(fileName,n1,n2,n3=1):
  if n3==1:
    x = zerofloat(n1,n2)
  else:
    x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def plotSlice(st,sx,p,tmin=None,tmax=None,amin=None,amax=None,
    title=None,pngDir=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    if pngDir==None:
      sp.setTitle(title)
  sp.setHLabel("Offset (km)")
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  sp.addColorBar()
  if tmin!=None and tmax!=None:
    sp.setVLimits(tmin,tmax)
  pv = sp.addPixels(st,sx,p)
  #pv.setColorModel(ColorMap.getRedWhiteBlue())
  if amax:
    if amin:
      pv.setClips(amin,amax)
    else:
      pv.setClips(-amax,amax)
  if pngDir:
    sp.paintToPng(1000,5,pngDir+title+".png")

def balance(sigma,f):
  f = add(max(f)*0.00001,f)
  ff = mul(f,f)
  RecursiveExponentialFilter(sigma).apply1(ff,ff)
  return div(f,sqrt(ff))

def sexp(x):
  return mul(sgn(x),log(add(abs(x),1.0)))

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



#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

