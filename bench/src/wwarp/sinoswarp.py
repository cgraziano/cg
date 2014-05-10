#########################################################
#Early script to test sinopec data.
#########################################################

from imports import *
from edu.mines.jtk.dsp.Conv import *
from wwarp import Warp, WaveletWarping, ShapingFilter
from testing import Synthetic
from linalgebra import ToeplitzRecursion
from java.util import Random


#pngDir = "./png/figures/"
pngDir = None

def main(args):
  #na = [9,9,9,8,7,13,13,30]
  #ka = [5,4,2,2,2,3,2,5]
  #nna = len(na)
  #for na in range(20,40):
  #  ka = int(na/2)+1
  #  enhancedWarping(na,-ka)
  #for ina in range(0,nna):
  #  enhancedWarping(na[ina],-ka[ina])
  enhancedWarping()

  #enhancedWarpingSyntheticShifts()

def enhancedWarping():
  ##################Parameters##############################
  #Sampling of data
    #PP,PSwarped,shifts
  ntf,dtf,ftf = 501,0.004,0.0 
  nx,dx,fx = 721,0.0150,0.000
  sft = Sampling(ntf,dtf,ftf)
  sx = Sampling(nx,dx,fx)
    #PS
  ntg,dtg,ftg = 852,0.004,0.0 
  sgt = Sampling(ntg,dtg,ftg)
  #Time window for wavelet estimation
  itminf,itmaxf = 124,250#140,313
  itming,itmaxg = 194,375#290,560
  #Space window for wavelet estimation
  ixmin,ixmax = 475,500
  #Frequency window for wavelet estimation
  fmin,fmax = 10*dtf,50*dtg
  #Inverse wavelet parameters estimated
  na,ka = 25,-12
  #na,ka = 81,-40
  #Wavelet parameters estimated
  nh,kh = 51,-25
  #Stability
  sfac = 1.00
  #amax for display
  amax = 5
  #max amount of squeezing seen in the data
  r = 1.7
  ##########################################################

  f = getSubPP()#ntf = 501
  g = getSubPS()#ntg = 852
  rightsg = getSubPSWarped()#501
  shifts = getShiftImages()#501

  ww = WaveletWarping()
  #ww.setFrequencyRange(fmin,fmax)
  ww.setTimeRange(itminf,itmaxf)
  ww.setSpaceRange(ixmin,ixmax)
  ww.setStabilityFactor(sfac)
  u = ww.getWarpedSamples(sgt,shifts)
  nu = len(u[0])

  ix = 350
  rwarp = zerofloat(nu)
  for iu in range(1,nu):
    rwarp[iu] = u[ix][iu] - u[ix][iu-1]
  SimplePlot.asPoints(rwarp)
  ssgt = Sampling(nu,0.004,0)

##############StartTests######################
  #starting traces
  plotTrace(sft,f[ix],"f",itmin=itminf,itmax=itmaxf,amax=amax)
  plotTrace(sgt,g[ix],"g",itmin=itming,itmax=itmaxg,amax=amax)

  #balance
  #f = balance(60,f)
  #g = balance(60,g)
  plotImage(sft,sx,f,itmin=0,itmax=500,\
  ixmin=0,ixmax=720,title="balanced f")
  plotImage(sgt,sx,g,itmin=0,itmax=851,\
  ixmin=0,ixmax=720,title="balanced g")
  plotImage(sft,sx,f,itmin=itminf,itmax=itmaxf,\
  ixmin=ixmin,ixmax=ixmax,title="windowed balanced f")
  plotImage(sgt,sx,g,itmin=itming,itmax=itmaxg,\
  ixmin=ixmin,ixmax=ixmax,title="windowed balanced g")
  plotTrace(sft,f[ix],"balanced f",itmin=itminf,itmax=itmaxf,amax=amax)
  plotTrace(sgt,g[ix],"balanced g",itmin=itming,itmax=itmaxg,amax=amax)

  #apply L
  ww.setR(r)
  lg = ww.applyL(u,g)
  plotTrace(sgt,lg[ix],"Lg",itmin=itming,itmax=itmaxg)

  #Test correct warping of g
  slg = ww.applyS(u,lg)
  plotAmplitudeSpectrumT(sgt,lg[ix],itming,itmaxg,"Lg",amax=60)
  plotAmplitudeSpectrumT(sft,slg[ix],itminf,itmaxf,"SLg",amax=60)
  plotAmplitudeSpectrumT(sft,f[ix],itminf,itmaxf,"f",amax=60)
  plotTrace(sft,slg[ix],"SLg",itmin=itminf,itmax=itmaxf,amax=amax)
  plotTrace(sft,f[ix],"f",itmin=itminf,itmax=itmaxf,amax=amax)


  #apply B to f and slg
  #bslg = ww.applyB(slg)
  #bf = ww.applyB(f)
  #plotTrace(sft,bslg[ix],"BSLg",itmin=itminf,itmax=itmaxf)
  #plotTrace(sft,bf[ix],"Bf",itmin=itminf,itmax=itmaxf)
  #Differences
  #d = sub(bslg,bf)
  #plotTrace(sft,d[ix],"d",itmin=itminf,itmax=itmaxf)

  #Apply B to differences, not to f and slg seperately
  d = sub(slg,f)
  bd = ww.applyB(d)
  plotTrace(sft,d[ix],"d",itmin=itminf,itmax=itmaxf)
  plotTrace(sft,bd[ix],"bd",itmin=itminf,itmax=itmaxf)


  #Estimated Wavelet
  #a = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
  niter = 10
  #a = ww.getInverseAIter(na,ka,nh,kh,niter,u,f,g) # estimated inverse wavelet
  a = ww.getInverseAIter2(na,ka,nh,kh,niter,u,f,g) # estimated inverse wavelet
  ga = ww.applyA(na,ka,a,g)
  lga = ww.applyL(u,ga)
  slga = ww.applyS(u,lga)
  h = ShapingFilter.design(nh,kh,nu,0,slga,nu,0,f)

  #h = ww.getWaveletH(na,ka,a,nh,kh) # estimated wavelet
  title = "Estimated Wavelet na = "+str(na)+" ka = "+str(ka)
  plotWavelets(Sampling(nh,sft.getDelta(),kh*sft.getDelta()),
    [h],title=title,pngDir=pngDir,onecol=True,twocol=True)
  title = "Estimated Inverse Wavelet na = "+str(na)+" ka = "+str(ka)
  plotWavelets(Sampling(na,sft.getDelta(),ka*sft.getDelta()),
    [a],title=title,pngDir=pngDir,onecol=True,twocol=True)

  #Test
  hag = ww.applyH(nh,kh,h,ww.applyA(na,ka,a,g))
  plotImage(sgt,sx,hag,itmin=itming,itmax=itmaxg,\
  ixmin=ixmin,ixmax=ixmax,amax=amax,title="hag")

  hslag = ww.applyHSLA(na,ka,a,nh,kh,h,u,g)
  plotImage(sft,sx,hslag,itmin=0,itmax=501,\
  ixmin=0,ixmax=720,amax=amax,title="HSLAg")
  plotImage(sft,sx,hslag,itmin=itminf,itmax=itmaxf,\
  ixmin=ixmin,ixmax=ixmax,amax=.1,title="HSLAg")
  plotImage(sft,sx,f,itmin=itminf,itmax=itmaxf,\
  ixmin=ixmin,ixmax=ixmax,amax=amax,title="f")
  plotImage(sft,sx,slg,itmin=itminf,itmax=itmaxf,\
  ixmin=ixmin,ixmax=ixmax,amax=amax,title="SLg")

  #apply BSLGa
  bslga = ww.applyBSLA(na,ka,a,u,g)
  plotImage(sft,sx,bslga,itmin=itminf,itmax=itmaxf,\
  ixmin=ixmin,ixmax=ixmax,amax=.5,title="BSLGa")
  
  #apply BFa
  bfa = ww.applyA(na,ka,a,f)
  plotImage(sft,sx,bfa,itmin=itminf,itmax=itmaxf,\
  ixmin=ixmin,ixmax=ixmax,amax=.5,title="BFa")


  diffEnh = sub(f,hslag)
  rmsDE = rms(diffEnh,ixmin,ixmax,itminf,itmaxf)
  print "Enhanced RMS Difference",rmsDE

  #apply Bf
  bf = ww.applyB(f)
  plotImage(sft,sx,bf,itmin=itminf,itmax=itmaxf,\
  ixmin=ixmin,ixmax=ixmax,title="Bf")

  #apply BSLg
  bslg = ww.applyBSL(na,ka,a,u,g)
  plotImage(sft,sx,bslg,itmin=itminf,itmax=itmaxf,\
  ixmin=ixmin,ixmax=ixmax,title="BSLg")

  sg = ww.applyS(u,g)
  diffReg = sub(f,sg)
  rmsDR = rms(diffReg,ixmin,ixmax,itminf,itmaxf)
  print "Regular RMS Difference",rmsDR

  #convolve two wavelets together
  nt = (na+ka-1)+(nh+kh-1)+-(ka+kh)
  hg = zerofloat(nt)
  conv(nh,kh,h,na,ka,a,nt,(ka+kh),hg)
  title = "h*a na = "+str(na)+" ka = "+str(ka)
  plotTrace(Sampling(nt,dtf,(kh+ka)*dtf),hg,title)
  #amax = 4
  #ia = 40
  #bf,bg,bsdg,bd,bbd = True,True,True,True,True
  #ww.plotDiffTraces(sx,sft,sgt,itminf,itmaxf,itming,itmaxg,ix,amax,\
  #na,ka,ia,bf,bg,bsdg,bd,bbd)
  bf,bg,bldg,bsldg,bd,bbd = False,False,False,False,False,True
  #ww.plotDiffGathers(sx,sft,sgt,itminf,itmaxf,itming,itmaxg,ixmin,ixmax,\
  #amax,na,ka,bf,bg,bldg,bsldg,bd,bbd)

  #print "Improved variance: "+str(ww.getVarianceEnhancedWarp(na,ka,a,\
  #sft,sx,u,g,f))
  #print "Regular variance: "+str(ww.getVarianceRegularWarp(na,ka,a,\
  #sft,sx,u,g,f))

  #df,dg,sdg,d,bd = True,False,True,True,True
  #ia = 39 
  #lg = ww.applyL(u,g)


  #plotAmplitudeSpectrumT(sgt,g[ix],itming,itmaxg,"g")
  #plotAmplitudeSpectrumT(sgt,lg[ix],itming,itmaxg,"lg")
  #sg = ww.applyS(u,lg)
  #plotAmplitudeSpectrumT(sft,sg[ix],itminf,itmaxf,"sg")
  #plotTrace(sft,f[ix],"f",itmin=itminf,itmax=itmaxf)
  #plotTrace(sgt,g[ix],"g",itmin=itming,itmax=itmaxg)

  #plotImage(sft,sx,f,amax=damax,itmin=itminf,itmax=itmaxf,\
  #ixmin=ixmin,ixmax=ixmax,title="f")
  #plotImage(sgt,sx,g,amax=damax,itmin=itming,itmax=itmaxg,\
  #ixmin=ixmin,ixmax=ixmax,title="g")

  #h = ww.getWaveletH(na,ka,a,nh,kh)
  #hsyn = getWavelet(fpeak,decay,nh,kh,mp) # synthetic/known wavelet

  #title = "Estimated Wavelet_"+str(r)+"_time"
  #plotWavelets(Sampling(nh,st.getDelta(),kh*st.getDelta()),
  #  [hsyn,h],title=title,pngDir=pngDir,onecol=True,twocol=True)
  
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

def enhancedWarpingP(na,ka):
  ##################Parameters##############################
  #Sampling of data
    #PP,PSwarped,shifts
  ntf,dtf,ftf = 501,0.004,0.0 
  nx,dx,fx = 721,0.0150,0.000
  sft = Sampling(ntf,dtf,ftf)
  sx = Sampling(nx,dx,fx)
    #PS
  ntg,dtg,ftg = 852,0.004,0.0 
  sgt = Sampling(ntg,dtg,ftg)
  #Time window for wavelet estimation
  itminf,itmaxf = 140,313#int(1.0/dt),int(1.5/dt)
  itming,itmaxg = 290,560
  #Space window for wavelet estimation
  ixmin,ixmax = 290,410
  #Frequency window for wavelet estimation
  fmin,fmax = 5*dtf,85*dtg
  #Inverse wavelet parameters estimated
  #na,ka = 81,-40
  #na,ka = 10,-4
  na,ka = na,ka
  #Wavelet parameters estimated
  nh,kh = 181,-90
  #Stability
  sfac = 1.00
  #amax for display
  amax = 5
  #max amount of squeezing seen in the data
  r = 1.5
  ##########################################################

  f = getSubPP()#ntf = 501
  g = getSubPS()#ntg = 852
  rightsg = getSubPSWarped()#501
  shifts = getShiftImages()#501

  ww = WaveletWarping()
  ww.setFrequencyRange(fmin,fmax)
  ww.setTimeRange(itminf,itmaxf)
  ww.setSpaceRange(ixmin,ixmax)
  ww.setStabilityFactor(sfac)
  u = ww.getWarpedSamples(sgt,shifts)
  nu = len(u[0])

  ix = 350
  #for iu in range(1,nu):
  #  print u[ix][iu]
  ssgt = Sampling(nu,0.004,0)

##############StartTests######################
  #starting traces
  #plotTrace(sft,f[ix],"f",itmin=itminf,itmax=itmaxf)
  #plotTrace(sgt,g[ix],"g",itmin=itming,itmax=itmaxg)

  #balance
  f = balance(30,f)
  g = balance(30,g)
  #plotTrace(sft,f[ix],"balanced f",itmin=itminf,itmax=itmaxf,amax=amax)
  #plotTrace(sgt,g[ix],"balanced g",itmin=itming,itmax=itmaxg,amax=amax)

  #apply L
  ww.setR(r)
  lg = ww.applyL(u,g)
  #plotTrace(sgt,lg[ix],"Lg",itmin=itming,itmax=itmaxg)

  #Test correct warping of g
  slg = ww.applyS(u,lg)
  #plotAmplitudeSpectrumT(sgt,lg[ix],itming,itmaxg,"Lg",amax=60)
  #plotAmplitudeSpectrumT(sft,slg[ix],itminf,itmaxf,"SLg",amax=60)
  #plotTrace(sft,slg[ix],"SLg",itmin=itminf,itmax=itmaxf,amax=amax)
  #plotTrace(sft,f[ix],"f",itmin=itminf,itmax=itmaxf,amax=amax)


  #apply B to f and slg
  #bslg = ww.applyB(slg)
  #bf = ww.applyB(f)
  #plotTrace(sft,bslg[ix],"BSLg",itmin=itminf,itmax=itmaxf)
  #plotTrace(sft,bf[ix],"Bf",itmin=itminf,itmax=itmaxf)
  #Differences
  #d = sub(bslg,bf)
  #plotTrace(sft,d[ix],"d",itmin=itminf,itmax=itmaxf)

  #Apply B to differences, not to f and slg seperately
  d = sub(slg,f)
  bd = ww.applyB(d)
  #plotTrace(sft,bd[ix],"bd",itmin=itminf,itmax=itmaxf)


  #Estimated Wavelet
  a = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
  h = ww.getWaveletH(na,ka,a,nh,kh) # estimated wavelet
  title = "Estimated Wavelet na = "+str(na)+" ka = "+str(ka)
  plotWavelets(Sampling(nh,sft.getDelta(),kh*sft.getDelta()),
    [h],title=title,pngDir=pngDir,onecol=True,twocol=True)

  #hslag = ww.applyHSLA(na,ka,a,nh,kh,h,u,g)
  #plotImage(sft,sx,hslag,itmin=itminf,itmax=itmaxf,\
  #ixmin=ixmin,ixmax=ixmax,title="HSLAg")
  nt = len(h)
  g = zerofloat(nt)
  conv(nh,kh,h,na,ka,a,nt,0,g)
  plotTrace(Sampling(nt,dtf,(kh+ka)*dtf),g,title)
  #amax = 4
  #ia = 40
  #bf,bg,bsdg,bd,bbd = True,True,True,True,True
  #ww.plotDiffTraces(sx,sft,sgt,itminf,itmaxf,itming,itmaxg,ix,amax,\
  #na,ka,ia,bf,bg,bsdg,bd,bbd)

  #df,dg,sdg,d,bd = True,False,True,True,True
  #ia = 39 
  #lg = ww.applyL(u,g)


  #plotAmplitudeSpectrumT(sgt,g[ix],itming,itmaxg,"g")
  #plotAmplitudeSpectrumT(sgt,lg[ix],itming,itmaxg,"lg")
  #sg = ww.applyS(u,lg)
  #plotAmplitudeSpectrumT(sft,sg[ix],itminf,itmaxf,"sg")
  #plotTrace(sft,f[ix],"f",itmin=itminf,itmax=itmaxf)
  #plotTrace(sgt,g[ix],"g",itmin=itming,itmax=itmaxg)

  #plotImage(sft,sx,f,amax=damax,itmin=itminf,itmax=itmaxf,\
  #ixmin=ixmin,ixmax=ixmax,title="f")
  #plotImage(sgt,sx,g,amax=damax,itmin=itming,itmax=itmaxg,\
  #ixmin=ixmin,ixmax=ixmax,title="g")

  #h = ww.getWaveletH(na,ka,a,nh,kh)
  #hsyn = getWavelet(fpeak,decay,nh,kh,mp) # synthetic/known wavelet

  #title = "Estimated Wavelet_"+str(r)+"_time"
  #plotWavelets(Sampling(nh,st.getDelta(),kh*st.getDelta()),
  #  [hsyn,h],title=title,pngDir=pngDir,onecol=True,twocol=True)
  
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


def enhancedWarpingSyntheticShifts():
  ##################Parameters##############################
  #Sampling of data
    #PP,PSwarped,shifts
  ntf,dtf,ftf = 501,0.004,0.0 
  nx,dx,fx = 721,0.0150,0.000
  sft = Sampling(ntf,dtf,ftf)
  sx = Sampling(nx,dx,fx)
    #PS
  ntg,dtg,ftg = 852,0.004,0.0 
  sgt = Sampling(ntg,dtg,ftg)
  #Time window for wavelet estimation
  itminf,itmaxf = 140,313#int(1.0/dt),int(1.5/dt)
  itming,itmaxg = 290,560
  #Space window for wavelet estimation
  ixmin,ixmax = 290,410
  #Frequency window for wavelet estimation
  fmin,fmax = 5*dtf,85*dtg
  #Inverse wavelet parameters estimated
  na,ka = 81,-40
  #Wavelet parameters estimated
  nh,kh = 181,-90
  #Stability
  sfac = 1000000
  #Balance
  tbal = 100
  #amax for display
  damax = 3.5
  #max amount of squeezing seen in the data
  r =  2
  #position
  ix = 350
  ##########################################################

  f = getSubPP()#ntf = 501
  g = getSubPS()#ntg = 852
  rightsg = getSubPSWarped()#501
  shifts = getShiftImages()#501

  ww = WaveletWarping()
  ww.setFrequencyRange(fmin,fmax)
  ww.setTimeRange(itminf,itmaxf)
  ww.setSpaceRange(ixmin,ixmax)
  ww.setStabilityFactor(sfac)
  u = ww.getWarpedSamples(sgt,shifts)
  #aw = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
  #hw = ww.getWaveletH(na,ka,aw,nh,kh) # estimated wavelet
  nu = len(u[0])
  print nu
  print len(g[0])
  ix = 350
  u = rampfloat(0,r,0,nu,nx)#constant factor of 2
  for iu in range(itminf,itmaxf):
    print u[ix][iu]-u[ix][iu-1]
  ssgt = Sampling(nu,0.004,0)


##############StartTests######################
  #starting trace
  plotTrace(sgt,g[ix],"g",itmin=itming,itmax=itmaxg)

  #apply L
  ww.setR(r)
  lg = ww.applyL(u,g)
  plotTrace(sgt,lg[ix],"Lg",itmin=itming,itmax=itmaxg)

  #Test correct warping of g
  slg = ww.applyS(u,lg)
  plotAmplitudeSpectrumT(sgt,lg[ix],itming,itmaxg,"Lg",amax=60)
  plotAmplitudeSpectrumT(ssgt,slg[ix],int(itming/r),int(itmaxg/r),\
  "SLg",amax=60)
  plotTrace(ssgt,slg[ix],"SLg",itmin=int(itming/r),itmax=int(itmaxg/r))
  plotTrace(sgt,lg[ix],"Lg",itmin=itming,itmax=itmaxg)

  #apply B to f and slg
  bslg = ww.applyB(slg)
  bg = ww.applyB(g)
  plotTrace(ssgt,bslg[ix],"BSLg",itmin=int(itming/r),itmax=int(itmaxg/r))

  #Differences
  #d = sub(bslg,bg)
  #plotTrace(sgt,d[ix],"d",itmin=itminf,itmax=itmaxf)

  #Apply B to differences, not to f and slg seperately
  #d = sub(slg,bg)
  #bd = ww.applyB(d)
  #plotTrace(sgt,bd[ix],"bd",itmin=itminf,itmax=itmaxf)


  #Estimated Wavelet
  #a = ww.getInverseA(na,ka,u,f,g) # estimated inverse wavelet
  #h = ww.getWaveletH(na,ka,a,nh,kh) # estimated wavelet
  #amax = 4
  #ia = 40
  #bf,bg,bsdg,bd,bbd = True,True,True,True,True
  #ww.plotDiffTraces(sx,sft,sgt,itminf,itmaxf,itming,itmaxg,ix,amax,\
  #na,ka,ia,bf,bg,bsdg,bd,bbd)



def plotImage(st, sx, f, itmin=None, itmax=None, ixmin=None, ixmax=None,
  amin=None, amax=None, title=None, pngDir=None, onecol=None, twocol=None):
  pv = PixelsView(st,sx,f)
  pv.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  if amax and amin:
    pv.setClips(amin,amax)
  elif amax:
    pv.setClips(-amax,amax)
  
  pp = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv)
  pp.setHLabel(0,"Position (km)")
  pp.setVLabel("Time (s)")
  if itmin:
    dt = st.getDelta()
    tmin = itmin*dt
    tmax = itmax*dt
    pp.setVLimits(tmin,tmax)
  if ixmin:
    dx = sx.getDelta()
    xmin = ixmin*dx
    xmax = ixmax*dx
    pp.setHLimits(xmin,xmax)
  
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



def plot2TracesSideBySide(st, f, g, itmin, itmax, 
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
  dt = st.getDelta()
  tmin = itmin*dt
  tmax = itmax*dt
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

def plot3TracesSideBySide(st, f, g, h, itmin, itmax, 
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
  dt = st.getDelta()
  tmin = itmin*dt
  tmax = itmax*dt
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

def plotTrace(st, p, title, itmin=None, itmax=None , amax=None): 
  dt = st.getDelta()
  ft = st.getFirst()
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  if itmin:
    tmin = ft+itmin*dt
    tmax = ft+itmax*dt
    tmin = itmin*dt
    tmax = itmax*dt
    sp.setVLimits(tmin,tmax)
  if amax:
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

def balance(sigma,f):
  f = add(max(f)*0.00001,f)
  ff = mul(f,f)
  RecursiveExponentialFilter(sigma).apply1(ff,ff)
  return div(f,sqrt(ff))

def getSinoImages():
  dataDir = "/Users/Chris/data/sinos/"
  n1,d1,f1 = 2001,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"z260.dat",n1,n2)
  g = readImage(dataDir+"x260.dat",n1,n2)
  gain(100,f)
  gain(100,g)
  return f,g

def getSubPP():
  dataDir = "/Users/Chris/data/sinos/"
  n1,d1,f1 = 501,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1,n2)
  gain(100,f)
  return f

def getSubPS():
  dataDir = "/Users/Chris/data/sinos/"
  n1,d1,f1 = 852,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"ps.dat",n1,n2)
  gain(100,f)
  return f

def getSubPSWarped():
  dataDir = "/Users/Chris/data/sinos/"
  n1,d1,f1 = 501,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  g = readImage(dataDir+"pswarped.dat",n1,n2)
  gain(100,g)
  return g

def getShiftImages():
  dataDir = "/Users/Chris/data/sinos/"
  n1,d1,f1 = 501,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  shifts = readImage(dataDir+"shifts.dat",n1,n2)
  return shifts 

def gain(hw,f):
  """ normalize RMS amplitude within overlapping windows, half-width hw """
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def normalize(h):
  return div(h,max(max(h),-min(h)))

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

#Plot amplitude spectrum for a gather trace
def plotAmplitudeSpectrumG(st, p, ixmin, ixmax, itmin, itmax, nnfft, title):
  #Time sampling for the time window specified.
  nt = ithi-itlo
  dt = st.getDelta()
  ft = st.getValue(itmin)

  nx = ixmax-ixmin

  subp = zerofloat(nt,nx)
  subst = Sampling(nt,dt,ft)
  for ix in range(0,nx):
    for it in range(0,nt):
      subp[ix][it] = p[ix][itmin+it]

  #Frequency sampling
  nfft = FftReal.nfftSmall(2*nt)#more time sample, the finer freq. samples
  nf = nfft/2+1
  df = 1.0/(nfft*dt)
  ff = 0.0
  fs = Sampling(nf,df,ff)
  ampSum = zerofloat(nf)
  for ix in range(0,nx):
    amp = computeAmplitudeSpectrum(subst,fs,nfft,p)
    ampSum = add(ampSum,amp)
  ampMean = div(ampSum,nx)
  plotSpectrum(fs,ampMean,title)

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

def rms(f,i2min,i2max,i1min,i1max):
  n2 = i2max-i2min
  n1 = i1max-i1min
  sum = 0.0
  for i2 in range(i2min,i2max):
    for i1 in range(i1min,i1max):
      sum = sum + f[i2][i1]*f[i2][i1]
  mean = sum/(n2*n1)
  return sqrt(mean)






#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())



