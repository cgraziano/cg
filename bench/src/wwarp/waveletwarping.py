#############################################################################
# Demo wavelet estimation from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from wwarp import WaveletWarping
from wwarp import ShapingFilter
from java.util import Random

############################################################################

#pngDir = "./png/figures/"
pngDir = "./png/figuresAbs/"
#pngDir = None

def main(args):
  goSimpleTest()

  #goSFacTest()

  fnrms,gnrms = .0001,.1
  #goNoiseTest(fnrms,gnrms,fNoise=True,gNoise=True)

  fdecay,gdecay = 0.02,0.05
  ffreq,gfreq = 0.10,0.08
  #goDiffWaveletsTest(fdecay,gdecay,ffreq,gfreq)

  rcor = [0.5,2]
  rincor = [0.51,2.01]#.51 gets error,.501 does not
  #goIncorrectWarp(rcor,rincor)

def goSimpleTest():
  nt,ni = 481,2 # number of time samples; number of impulses
  freq,decay = 0.08,0.05 # peak frequency and decay for wavelet
  na,ka = 10,-2#81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  amax = 5
  st = Sampling(nt,dt,ft)
  for mp in [False]: # True, for minimum-phase; False for other
    hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
    for r in [2]: # for stretch and squeeze, ...
      fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
      p,q = makeImpulses(r,nt,ni)
      title = "impulsespsqq"
      f = addWavelet(freq,decay,p,mp)
      g = addWavelet(freq,decay,q,mp)
      #f = addNoise(0.001,42,f)
      #g = addNoise(0.001,42,g)
      u = rampfloat(0.0,r,nt)#PS time as a function of PP time
      if r<=1.0:
        u = add((1.0-r)*(nt-1),u)
        #u is used in warp, which is an interpolator. In the case of
        #r<=1, if this 
        #if statement block was not here, the values to interpolate at (u)
        #woudl not include the wavelets.
      
      #main piece of wavelet estimation code
      ww = WaveletWarping()
      ww.setFrequencyRange(fmin,fmax)
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)
      lq = ww.applyL(u,q)
      sq = ww.applyS(u,lq)
      #plot3TracesSideBySide(st,p,q,sq,tmin*dt,1.6,amax,title=title,pngDir=pngDir,onecol=True)
      ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
      aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
      #SimplePlot.asPoints(aw)
      hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet

      dump(ak)
      dump(aw)
      dump(hw)
      sg = ww.applyS(u,g)
      lg = ww.applyL(u,g)
      slg = ww.applyS(u,lg)
      af = ww.applyA(na,ka,aw,f)
      ag = ww.applyA(na,ka,aw,g)
      lag = ww.applyL(u,ag) # lowpass, if squeezing
      slag = ww.applyS(u,lag)
      baf = ww.applyB(af)
      bslag = ww.applyB(slag)
      hslag = ww.applyH(nh,kh,hw,slag)
      nhw = normalize(hw)
      nhk = normalize(hk)
      title = "r = "+str(r)
      plotSequences(st,[p,q],labels=["p","q"],title=title)
      plotSequences(st,[p,sq],labels=["p","Sq"],title=title)
      plotSequences(st,[p,sq,sub(p,sq)],labels=["p","Sq","diff"],title=title)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[sg,slg],labels=["sg","slg"],title=title)
      plotSequences(st,[af,ag],labels=["Af","Ag"],title=title)
      plotSequences(st,[af,lag],labels=["Af","LAg"],title=title)
      plotSequences(st,[af,slag],labels=["Af","SLAg"],title=title)
      plotSequences(st,[baf,bslag],labels=["BAf","BSLAg"],title=title)
      plotSequences(st,[f,hslag],labels=["f","HSLAg"],title=title)
      plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],title=title)
      #Impulse Analysis
      #########################################
      #sq = ww.applyS(u,q)
      #bsq = ww.applyB(sq)
      #bp = ww.applyB(p)
      #title = "r = "+str(r)
      #plotSequences(st,[p,q],labels=["p","q"],title=title)
      #plotSequences(st,[p,sq],labels=["p","Sq"],title=title)
      #plotSequences(st,[bp,bsq],labels=["Bp","BSq"],title=title)

      

      title = "initialTraces r = "+str(r)+" mp = "+str(mp)

      amax = 5 
      #plot3TracesSideBySide(st,f,g,sg,tmin*dt,1.6,amax,title=title,
      #pngDir=pngDir,onecol=True)
      title = "warpedComparison r = "+str(r)+" mp = "+str(mp)

      #plot3TracesSideBySide(st,f,hslag,sg,tmin*dt,1.6,
      #amax,title=title,pngDir=pngDir,onecol=True)
      
      title = "estimatedWavelet r = "+str(r)+" mp = "+str(mp)
      #plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],title=title,
      #pngDir=pngDir,onecol=True)
      amax = 5
 
      SimplePlot.asPoints(f)
def goSFacTest():
  nt,ni = 481,2 # number of time samples; number of impulses
  freq,decay = 0.08,0.05 # peak frequency and decay for wavelet
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  st = Sampling(nt,dt,ft)
  amax = 2
  #for mp in [True,False]: # True, for minimum-phase; False for other
  for mp in [False]: # True, for minimum-phase; False for other
    hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
    #for r in [0.5,2.0]: # for stretch and squeeze, ...
    for r in [.5]: # for stretch and squeeze, ...
      for sfac in [1.00,1.00001,1.0001,1.001,1.01,1.1,1.5]:
        fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
        p,q = makeImpulses(r,nt,ni)
        f = addWavelet(freq,decay,p,mp)
        g = addWavelet(freq,decay,q,mp)

        #diagnostic
        st = Sampling(nt,dt,0.00)
        title = "starting traces"
        #plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,2,title=title)
        
        u = rampfloat(0.0,r,nt)
        print "u before mod"
        #dump(u)
        if r<=1.0:
          u = add((1.0-r)*(nt-1),u)
        print "u after mod"
        #dump(u)
        ww = WaveletWarping()
        ww.setFrequencyRange(fmin,fmax)
        ww.setTimeRange(tmin,tmax)
        ww.setStabilityFactor(sfac)
        ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
        aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
        hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
        #dump(ak)
        #dump(aw)

        sg = ww.applyS(u,g)
        af = ww.applyA(na,ka,aw,f)
        ag = ww.applyA(na,ka,aw,g)
        lag = ww.applyL(u,ag) # lowpass, if squeezing
        slag = ww.applyS(u,lag)
        baf = ww.applyB(af)
        bslag = ww.applyB(slag)
        hslag = ww.applyH(nh,kh,hw,slag)
        nhw = normalize(hw)
        nhk = normalize(hk)

        #title = " r = "+str(r)+" mp = "+str(mp)
        #plotSequences(st,[f,g],labels=["f","g"],title=title)
        #plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
        #plotSequences(st,[f,g],labels=["f","g"],title=title)
        #plotSequences(st,[af,ag],labels=["Af","Ag"],title=title)
        #plotSequences(st,[af,lag],labels=["Af","LAg"],title=title)
        #plotSequences(st,[af,slag],labels=["Af","SLAg"],title=title)
        #plotSequences(st,[baf,bslag],labels=["BAf","BSLAg"],title=title)
        #plotSequences(st,[f,hslag],labels=["f","HSLAg"],title=title)
        #title = "initialTraces"
        #plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,amax,title=title,
        #pngDir=pngDir,onecol=True)
        #title = "warpedComparison"
        #plot3TracesSideBySide(st,f,sg,hslag,tmin*dt,tmax*dt,
        #amax,title=title,pngDir=pngDir,onecol=True)
        
        title = "estimatedWavelet sfac = "+str(sfac)+\
        " r = "+str(r)+" mp = "+str(mp)
        plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],title=title,
        pngDir=pngDir,onecol=True)

def goNoiseTest(fnrms,gnrms,fNoise=None,gNoise=None):
  nt,ni = 481,2 # number of time samples; number of impulses
  freq,decay = 0.08,0.05 # peak frequency and decay for wavelet
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  st = Sampling(nt,dt,ft)
  amax = 2
  for mp in [True,False]: # True, for minimum-phase; False for other
  #for mp in [False]: # True, for minimum-phase; False for other
    hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
    for r in [0.5,1.5]: # for stretch and squeeze, ...
    #for r in [.5]: # for stretch and squeeze, ...
      fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
      p,q = makeImpulses(r,nt,ni)
      f = addWavelet(freq,decay,p,mp)
      g = addWavelet(freq,decay,q,mp)
      if fNoise:
        f = addNoise(fnrms,42,f)
      if gNoise:
        g = addNoise(gnrms,2,g)

      #diagnostic
      st = Sampling(nt,dt,0.00)
      title = "starting traces"
      plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,2,title=title)
      
      u = rampfloat(0.0,r,nt)
      print "u before mod"
      #dump(u)
      if r<=1.0:
        u = add((1.0-r)*(nt-1),u)
      print "u after mod"
      #dump(u)
      ww = WaveletWarping()
      ww.setFrequencyRange(fmin,fmax)
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)
      ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
      aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
      hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
      #dump(ak)
      #dump(aw)
      sg = ww.applyS(u,g)
      af = ww.applyA(na,ka,aw,f)
      ag = ww.applyA(na,ka,aw,g)
      lag = ww.applyL(u,ag) # lowpass, if squeezing
      slag = ww.applyS(u,lag)
      baf = ww.applyB(af)
      bslag = ww.applyB(slag)
      hslag = ww.applyH(nh,kh,hw,slag)
      nhw = normalize(hw)
      nhk = normalize(hk)
      title = "fnrms = "+str(fnrms)+" gnrms = "+str(gnrms)+" r = "+str(r)+\
      " mp = "+str(mp)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[af,ag],labels=["Af","Ag"],title=title)
      plotSequences(st,[af,lag],labels=["Af","LAg"],title=title)
      plotSequences(st,[af,slag],labels=["Af","SLAg"],title=title)
      plotSequences(st,[baf,bslag],labels=["BAf","BSLAg"],title=title)
      plotSequences(st,[f,hslag],labels=["f","HSLAg"],title=title)
      title = "initialTraces"
      #plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,amax,title=title,
      #pngDir=pngDir,onecol=True)
      #title = "warpedComparison"
      #plot3TracesSideBySide(st,f,sg,hslag,tmin*dt,tmax*dt,
      #amax,title=title,pngDir=pngDir,onecol=True)
      

      title = "estWave fnrms = "+str(fnrms)+" gnrms = "+str(gnrms)+" r = "+\
      str(r)+" mp = "+str(mp)
      plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],title=title,
      pngDir=pngDir,onecol=True)
 

def goDiffWaveletsTest(fdecay,gdecay,ffreq,gfreq):
  nt,ni = 481,2 # number of time samples; number of impulses
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  st = Sampling(nt,dt,ft)
  amax = 2
  #for mp in [True,False]: # True, for minimum-phase; False for other
  for mp in [False]: # True, for minimum-phase; False for other
    hk = getWavelet(gfreq,gdecay,nh,kh,mp) # known wavelet
    #for r in [0.5,2.0]: # for stretch and squeeze, ...
    for r in [.5]: # for stretch and squeeze, ...
      fmin,fmax = 0.0,min(0.5,0.5*r) # bandpass (lowpass), if stretching
      p,q = makeImpulses(r,nt,ni)
      f = addWavelet(ffreq,fdecay,p,mp)
      g = addWavelet(gfreq,gdecay,q,mp)

      #diagnostic
      st = Sampling(nt,dt,0.00)
      title = "starting traces"
      plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,2,title=title)
      
      u = rampfloat(0.0,r,nt)
      print "u before mod"
      #dump(u)
      if r<=1.0:
        u = add((1.0-r)*(nt-1),u)
      print "u after mod"
      #dump(u)
      ww = WaveletWarping()
      ww.setFrequencyRange(fmin,fmax)
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)
      ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
      aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
      hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
      #dump(ak)
      #dump(aw)
      sg = ww.applyS(u,g)
      af = ww.applyA(na,ka,aw,f)
      ag = ww.applyA(na,ka,aw,g)
      lag = ww.applyL(u,ag) # lowpass, if squeezing
      slag = ww.applyS(u,lag)
      baf = ww.applyB(af)
      bslag = ww.applyB(slag)
      hslag = ww.applyH(nh,kh,hw,slag)
      nhw = normalize(hw)
      nhk = normalize(hk)
      title = "r = "+str(r)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[af,ag],labels=["Af","Ag"],title=title)
      plotSequences(st,[af,lag],labels=["Af","LAg"],title=title)
      plotSequences(st,[af,slag],labels=["Af","SLAg"],title=title)
      plotSequences(st,[baf,bslag],labels=["BAf","BSLAg"],title=title)
      plotSequences(st,[f,hslag],labels=["f","HSLAg"],title=title)
      title = "initialTraces"
      #plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,amax,title=title,
      #pngDir=pngDir,onecol=True)
      #title = "warpedComparison"
      #plot3TracesSideBySide(st,f,sg,hslag,tmin*dt,tmax*dt,
      #amax,title=title,pngDir=pngDir,onecol=True)
      
      title = "estimatedWavelet"
      plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],title=title,
      pngDir=pngDir,onecol=True)
 
def goIncorrectWarp(rcor,rincor):
  nt,ni = 481,2 # number of time samples; number of impulses
  freq,decay = 0.08,0.05 # peak frequency and decay for wavelet
  na,ka = 81,-20 # sampling for inverse wavelet A
  nh,kh = 181,-90 # sampling for wavelet H
  dt,ft = 0.004,0.000 # used for plotting only
  tmin,tmax = 0,nt-1
  sfac = 1.000
  st = Sampling(nt,dt,ft)
  amax = 2
  #for mp in [True,False]: # True, for minimum-phase; False for other
  nr = len(rcor)
  for mp in [False]: # True, for minimum-phase; False for other
    hk = getWavelet(freq,decay,nh,kh,mp) # known wavelet
    for ir in range(0,nr): # for stretch and squeeze, ...
      fmin,fmax = 0.0,min(0.5,0.5*rincor[ir]) # bandpass (lowpass), if stretching
      p,q = makeImpulses(rcor[ir],nt,ni)
      f = addWavelet(freq,decay,p,mp)
      g = addWavelet(freq,decay,q,mp)

      #diagnostic
      st = Sampling(nt,dt,0.00)
      title = "starting traces"
      plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,2,title=title)
      
      u = rampfloat(0.0,rincor[ir],nt)
      print "u before mod"
      #dump(u)
      if rincor[ir]<=1.0:
        u = add((1.0-rincor[ir])*(nt-1),u)
      print "u after mod"
      #dump(u)
      ww = WaveletWarping()
      ww.setFrequencyRange(fmin,fmax)
      ww.setTimeRange(tmin,tmax)
      ww.setStabilityFactor(sfac)
      ak = ww.getWaveletH(nh,kh,hk,na,ka) # known inverse wavelet
      aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
      hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
      #dump(ak)
      #dump(aw)
      sg = ww.applyS(u,g)
      af = ww.applyA(na,ka,aw,f)
      ag = ww.applyA(na,ka,aw,g)
      lag = ww.applyL(u,ag) # lowpass, if squeezing
      slag = ww.applyS(u,lag)
      baf = ww.applyB(af)
      bslag = ww.applyB(slag)
      hslag = ww.applyH(nh,kh,hw,slag)
      nhw = normalize(hw)
      nhk = normalize(hk)
      title = " rcor = "+str(rcor[ir])+" rincor = "+str(rincor[ir])+\
      " mp = "+str(mp)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[f,sg],labels=["f","Sg"],title=title)
      plotSequences(st,[f,g],labels=["f","g"],title=title)
      plotSequences(st,[af,ag],labels=["Af","Ag"],title=title)
      plotSequences(st,[af,lag],labels=["Af","LAg"],title=title)
      plotSequences(st,[af,slag],labels=["Af","SLAg"],title=title)
      plotSequences(st,[baf,bslag],labels=["BAf","BSLAg"],title=title)
      plotSequences(st,[f,hslag],labels=["f","HSLAg"],title=title)
      #title = "initialTraces"
      #plot2TracesSideBySide(st,f,g,tmin*dt,tmax*dt,amax,title=title,
      #pngDir=pngDir,onecol=True)
      #title = "warpedComparison"
      #plot3TracesSideBySide(st,f,sg,hslag,tmin*dt,tmax*dt,
      #amax,title=title,pngDir=pngDir,onecol=True)
      
      title = "estWave rcor = "+str(rcor[ir])+" rincor "+str(rincor[ir])+\
      " mp = "+str(mp)
      plotWavelets(Sampling(nh,dt,kh*dt),[nhw,nhk],title=title,
      pngDir=pngDir,onecol=True)

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
  if not mp:#mp is minimum phase
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

def plotWavelets(st,hs,hmax=None,title=None,pngDir=None,
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
  sp.setSize(720,400)
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
  pp.setHInterval(0,2.0)
  pp.setHInterval(1,2.0)
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
  pp.setHInterval(0,2.0)
  pp.setHInterval(1,2.0)
  pp.setHInterval(2,2.0)
  pp.setVLimits(tmin,tmax)
  pp.setHLabel(0,"Amplitude")
  pp.setHLabel(1,"Amplitude")
  pp.setHLabel(2,"Amplitude")
  pp.setVLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setSize(720,400)
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
