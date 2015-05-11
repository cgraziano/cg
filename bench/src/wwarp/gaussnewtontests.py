#############################################################################
# Demo of 2 wavelet estimations from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.awt.ColorMap import *
from edu.mines.jtk.lapack import *
from wwarp import WaveletWarpingCBGN, Warper, Start
from wwarp import ShapingFilter 
import synthetic
import plotting
from java.util import Random

############################################################################

def main(args):
  testGNOnSynthetic1D()
  #testGNOnSynthetic2D()
  #testCommuteSynthetic()
  #testCommuteSino()
  #chooseAlphaTestsLessSqueezing()
  #chooseAlphaTestsMoreSqueezing()
  #chooseAlphaTestsLessSqueezingNoise()
  #chooseAlphaTestsMoreSqueezingNoise()
  #chooseAlphaTestsSino()


#A separate test to figure out if the Gauss-Newton method is working.
def testGNOnSynthetic1D():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  maxpc = 0.025#or0.001
  niter = 300
  ww = WaveletWarpingCBGN()
  ww.setMaxPercentChange(maxpc)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setLineSearchMinScale(0.0000)
  
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)

  #Get Gauss-Newton Information
  dataResRmsInitS = ww.getDataResRmsInitialS()
  dataResRmsFinaS = ww.getDataResRmsFinalS()
  dataRes2NormSqInitS = ww.getDataRes2NormSqInitialS()
  dataRes2NormSqFinaS = ww.getDataRes2NormSqFinalS()
  bPenaltyRes2NormSqInitS = ww.getBPenaltyRes2NormSqInitialS()
  bPenaltyRes2NormSqFinaS = ww.getBPenaltyRes2NormSqFinalS()
  cPenaltyRes2NormSqInitS = ww.getCPenaltyRes2NormSqInitialS()
  cPenaltyRes2NormSqFinaS = ww.getCPenaltyRes2NormSqFinalS()
  allResRmsInitS = ww.getAllResRmsInitialS()
  allResRmsFinaS = ww.getAllResRmsFinalS()
  allResRmsAllS = ww.getAllResRmsAllS()
  allRes2NormSqInitS = ww.getAllRes2NormSqInitialS()
  allRes2NormSqFinaS = ww.getAllRes2NormSqFinalS()
  stepLengthS = ww.getStepLengthS()
  conditionNumberS = ww.getConditionNumberS()
  b2NormS = ww.getB2NormS()
  c2NormS = ww.getC2NormS()
  deltaB2NormS = ww.getDeltaB2NormS()
  deltaC2NormS = ww.getDeltaC2NormS()
  gradient2NormS = ww.getGradient2NormS()
  alphaBS = ww.getAlphaB()
  alphaCS = ww.getAlphaC()
  rmsPercentChange = ww.getRmsPercentChangeS()
  lastIter = ww.getLastIter()

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  ak = ww.getWaveletC(nc,kc,ck,na,ka)
  bk = ww.getWaveletC(nd,kd,dk,nb,kb)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sg = warp.applyS(u,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))
  plotAmplitudeSpectrumT(Sampling(nb,0.004,0.004*kb), bguess, 0, nb, "b Spectrum")
  plotAmplitudeSpectrumT(Sampling(nb,0.004,0.004*kb), bw, 0, nb, "b Spectrum")
  plotAmplitudeSpectrumT(Sampling(nc,0.004,0.004*kc), dw, 0, nc, "d Spectrum")
  plotAmplitudeSpectrumT(Sampling(nc,0.004,0.004*kc), cw, 0, nc, "c Spectrum")
  plotAmplitudeSpectrumT(Sampling(nt,0.004,0.00), g, tmin, tmax, "g spectrum")
  plotAmplitudeSpectrumT(Sampling(nt,0.004,0.00), sub(g,noiseg), tmin, tmax, "No noise g spectrum")
  plotAmplitudeSpectrumT(Sampling(nt,0.004,0.00), bg, tmin, tmax, "Bg spectrum")
  plotAmplitudeSpectrumT(Sampling(nt,0.004,0.00), sbg, tmin, tmax, "SBg spectrum")
  plotAmplitudeSpectrumT(Sampling(nt,0.004,0.00), f, tmin, tmax, "f spectrum")

  #Print
  dt = 0.004
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)
  print "Number of time samples: "+str(nt)
  print "Number of reflection coefficients in f between tmin and tmax: "+str(ni)
  print "Random reflection coefficients: "+str(randomi)
  print "Reflection coefficients exist past "+\
  "tmax and corresponding u(tmax) in f and g, respecitively: "+str(moreps)
  print "Initial amount of squeezing: "+str(r0)
  print "Final amount of squeezing: "+str(r1)
  print "Shift between f and g: "+str(v)
  print "Wavelet c's peak frequency: "+str(freqc)
  print "Wavelet d's peak frequency: "+str(freqd)
  print "Wavelet c's decay: "+str(decayc)
  print "Wavelet d's decay: "+str(decayd)
  print "Wavelet c is minimum phase: "+str(mpc)
  print "Wavelet d is minimum phase: "+str(mpd)
  print "Noise to signal ratio in f and g: "+str(nrmsf)

  #Plotting
    #Data

  pngDir = None
  title= "[f,csbg] [f,hsg] [f,sg]"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-8.0
  hmax1,hmax2 = 1.0,8.0
  hlabel = ["[f,CSBg]","[f,HSg]","[f,Sg]"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,[[f,csbg],[f,hsg],[f,sg]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = None
  title= "f csbg p sbg bg q g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-8.0
  hmax1,hmax2 = 1.0,8.0
  hlabel = ["f","HSg","CSBg","p","SBg","Bg","q","g"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin1,hmax1],\
  [hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin2,hmax2]]
  hint = [None,None,None,None,None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,csbg,p,sbg,bg,q,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

    #GN Measurements
  title = "Total Initial RMS of Residuals (Black) Final RMS of Residuals (Blue)"
  for i in range(lastIter,niter):
    allResRmsInitS[i] = allResRmsInitS[lastIter]
    allResRmsFinaS[i] = allResRmsFinaS[lastIter]
  maxrmsri = max(allResRmsInitS)
  minrmsrf = min(allResRmsFinaS)
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [allResRmsInitS,allResRmsFinaS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)

  title = "Total RMS of Residuals (Black) "
  maxrmsri = max(allResRmsAllS)
  minrmsrf = min(allResRmsAllS)
  siter = Sampling(niter+1,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)


  title = "Total Initial 2NormSq of Residuals (Black) Final 2NormSq of Residuals (Blue)"
  maxrmsri = max(allRes2NormSqInitS)
  minrmsrf = min(allRes2NormSqFinaS)
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "2NormSq of residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [allRes2NormSqInitS,allRes2NormSqFinaS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)


  title = "Data Initial RMS of Residuals (Black) Final RMS of Residuals (Blue)"
  maxrmsri = max(dataResRmsInitS)
  minrmsrf = min(dataResRmsFinaS)
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [dataResRmsInitS,dataResRmsFinaS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)

  title = "Initial 2NormSq of individual and all residuals"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,
  [allRes2NormSqInitS,dataRes2NormSqInitS,bPenaltyRes2NormSqInitS,cPenaltyRes2NormSqInitS,alphaBS,alphaCS],\
  color=None,\
  vlabel=["allRes","dataRes","bPenRes","cPenRes","alphaB","alphaC"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  title = "Final 2NormSq of individual and all residuals"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,
  [allRes2NormSqFinaS,dataRes2NormSqFinaS,bPenaltyRes2NormSqFinaS,cPenaltyRes2NormSqFinaS,alphaBS,alphaCS],\
  color=None,\
  vlabel=["allRes","dataRes","bPenRes","cPenRes","alphaB","alphaC"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  title = "Step Length and Two Norm of Gradient"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,[stepLengthS,gradient2NormS],\
  color=None,\
  vlabel=["StepLength","2NormGrad"],\
  vminmax=[None,None],\
  vint=[None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  title = "2NormDeltaB,2NormDeltaC,2NormB,2NormC,ConditionNumber,RMS Percent Change"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,[deltaB2NormS,deltaC2NormS,\
  b2NormS,c2NormS,conditionNumberS,rmsPercentChange],\
  color=None,\
  vlabel=["2NormDelB","2NormDelC","2NormB","2NormC","CondNum","RMSPercChange"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  hsize=960,vsize=760,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)


    #Normalize
  ncguess = normalizeM(cguess)
  nbguess = normalizeM(bguess)
  ncw = normalizeM(cw)
  ndw = normalizeM(dw)
  nbw = normalizeM(bw)
  nck = normalizeM(ck)
  ndk = normalizeM(dk)
  nak = normalizeM(ak)
  nbk = normalizeM(bk)

    #Wavelets
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known c"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known d"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ndw,ndk],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known b"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbw,nbk],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

#A separate test to figure out if the Gauss-Newton method is working.
def testGNOnSynthetic2D():
  #Synthetic parameters
  nt,nx,ni,randomi = 1081,5,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  freqd,decayd = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn2D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nx,nt,ni,randomi,moreps)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  niter = 500
  ww = WaveletWarpingCBGN()
  ww.setMaxPercentChange(0.01)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setLineSearchMinScale(0.0000)
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  #Get Gauss-Newton Information
  dataResRmsInitS = ww.getDataResRmsInitialS()
  dataResRmsFinaS = ww.getDataResRmsFinalS()
  dataRes2NormSqInitS = ww.getDataRes2NormSqInitialS()
  dataRes2NormSqFinaS = ww.getDataRes2NormSqFinalS()
  bPenaltyRes2NormSqInitS = ww.getBPenaltyRes2NormSqInitialS()
  bPenaltyRes2NormSqFinaS = ww.getBPenaltyRes2NormSqFinalS()
  cPenaltyRes2NormSqInitS = ww.getCPenaltyRes2NormSqInitialS()
  cPenaltyRes2NormSqFinaS = ww.getCPenaltyRes2NormSqFinalS()
  allResRmsInitS = ww.getAllResRmsInitialS()
  allResRmsFinaS = ww.getAllResRmsFinalS()
  allRes2NormSqInitS = ww.getAllRes2NormSqInitialS()
  allRes2NormSqFinaS = ww.getAllRes2NormSqFinalS()
  stepLengthS = ww.getStepLengthS()
  conditionNumberS = ww.getConditionNumberS()
  b2NormS = ww.getB2NormS()
  c2NormS = ww.getC2NormS()
  deltaB2NormS = ww.getDeltaB2NormS()
  deltaC2NormS = ww.getDeltaC2NormS()
  gradient2NormS = ww.getGradient2NormS()
  alphaBS = ww.getAlphaB()
  alphaCS = ww.getAlphaC()
  rmsPercentChange = ww.getRmsPercentChangeS()
  lastIter = ww.getLastIter()

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  ak = ww.getWaveletC(nc,kc,ck,na,ka)
  bk = ww.getWaveletC(nd,kd,dk,nb,kb)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  #Print
  dt = 0.004
  dx = 0.015
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)
  print "Number of time samples: "+str(nt)
  print "Number of reflection coefficients in f between tmin and tmax: "+str(ni)
  print "Random reflection coefficients: "+str(randomi)
  print "Reflection coefficients exist past "+\
  "tmax and corresponding u(tmax) in f and g, respecitively: "+str(moreps)
  print "Initial amount of squeezing: "+str(r0)
  print "Final amount of squeezing: "+str(r1)
  print "Shift between f and g: "+str(v)
  print "Wavelet c's peak frequency: "+str(freqc)
  print "Wavelet d's peak frequency: "+str(freqd)
  print "Wavelet c's decay: "+str(decayc)
  print "Wavelet d's decay: "+str(decayd)
  print "Wavelet c is minimum phase: "+str(mpc)
  print "Wavelet d is minimum phase: "+str(mpd)
  print "Noise to signal ratio in f and g: "+str(nrmsf)

  #Plotting
    #Data
  pngDir = None
  title= "f csbg p sbg bg q g"
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hlabel = ["f(km)","HSg(km)","CSBg(km)","p(km)","SBg(km)","Bg(km)","q(km)","g(km)"]
  hminmax = None
  hint = None
  tilespacing = 5
  clip = 8.0
  plotting.plotImagesSideBySide(st,sx,[f,hsg,csbg,p,sbg,bg,q,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

    #GN Measurements
  title = "Total Initial RMS of Residuals (Black) Final RMS of Residuals (Blue)"
  maxrmsri = max(allResRmsInitS)
  minrmsrf = min(allResRmsFinaS)
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [allResRmsInitS,allResRmsFinaS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)

  title = "Total Initial 2NormSq of Residuals (Black) Final 2NormSq of Residuals (Blue)"
  maxrmsri = max(allRes2NormSqInitS)
  minrmsrf = min(allRes2NormSqFinaS)
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "2NormSq of residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [allRes2NormSqInitS,allRes2NormSqFinaS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)


  title = "Data Initial RMS of Residuals (Black) Final RMS of Residuals (Blue)"
  maxrmsri = max(dataResRmsInitS)
  minrmsrf = min(dataResRmsFinaS)
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [dataResRmsInitS,dataResRmsFinaS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)

  title = "Initial 2NormSq of individual and all residuals"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,
  [allRes2NormSqInitS,dataRes2NormSqInitS,bPenaltyRes2NormSqInitS,cPenaltyRes2NormSqInitS,alphaBS,alphaCS],\
  color=None,\
  vlabel=["allRes","dataRes","bPenRes","cPenRes","alphaB","alphaC"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  title = "Final 2NormSq of individual and all residuals"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,
  [allRes2NormSqFinaS,dataRes2NormSqFinaS,bPenaltyRes2NormSqFinaS,cPenaltyRes2NormSqFinaS,alphaBS,alphaCS],\
  color=None,\
  vlabel=["allRes","dataRes","bPenRes","cPenRes","alphaB","alphaC"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  title = "Step Length and Two Norm of Gradient"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,[stepLengthS,gradient2NormS],\
  color=None,\
  vlabel=["StepLength","2NormGrad"],\
  vminmax=[None,None],\
  vint=[None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  title = "2NormDeltaB,2NormDeltaC,2NormB,2NormC,ConditionNumber,RMS Percent Change"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,[deltaB2NormS,deltaC2NormS,\
  b2NormS,c2NormS,conditionNumberS,rmsPercentChange],\
  color=None,\
  vlabel=["2NormDelB","2NormDelC","2NormB","2NormC","CondNum","RMSPercChange"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  hsize=960,vsize=760,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)


    #Normalize
  ncguess = normalizeM(cguess)
  nbguess = normalizeM(bguess)
  ncw = normalizeM(cw)
  nbw = normalizeM(bw)
  nck = normalizeMAAWOS(ck)
  ndk = normalizeMAAWOS(dk)
  nak = normalizeMAAWOS(ak)
  nbk = normalizeMAAWOS(bk)

    #Wavelets
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known c"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known b"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbw,nbk],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  #A test to see if C and B commute with S with synthetics.
def testCommuteSynthetic():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  ##r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  r0,r1 = 1.6,1.4#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  freqd,decayd = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  alpha = 0.1
  penb = False
  penc = False
  pendeltab = False
  pendeltac = False
  niter = 5000
  ww = WaveletWarpingCBGN()
  ww.setMaxPercentChange(0.01)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setLineSearchMinScale(0.0000)
  ww.setStabilityFactor(0.0)
  ww.penalize(alpha,penb,penc,pendeltab,pendeltac)
  
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  #Get Gauss-Newton Information
  rmsri = ww.getRMSRi()
  rmsrf = ww.getRMSRf()
  twonormgrad = ww.getTwoNormGrad()
  steplength = ww.getStepLength()
  twonormdeltab = ww.getTwoNormDeltaB()
  twonormdeltac = ww.getTwoNormDeltaC()
  condnum = ww.getCondNum()
  twonormb = ww.getTwoNormB()
  twonormc = ww.getTwoNormC()
  rmspercentchange = ww.getRMSPercentChange()
  lastIter = ww.getLastIter()

  print "rms I love you Mom,Dad,Frank,Demi = "+str(ww.rms(ww.computeDataResidual(nc,kc,hw,nb,kb,bone,u,f,g)))
  

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  ak = ww.getWaveletC(nc,kc,ck,na,ka)
  bk = ww.getWaveletC(nd,kd,dk,nb,kb)

  #Processing
  warp = Warper()
  sg = warp.applyS(u,g)
  bsg = ww.applyC(nb,kb,bw,sg)
  cbsg = ww.applyC(nc,kc,cw,bsg)
  bg = ww.applyC(nb,kb,bw,g)
  cbg = ww.applyC(nc,kc,cw,bg)
  scbg = warp.applyS(u,cbg)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  #Print
  dt = 0.004
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)
  print "Number of time samples: "+str(nt)
  print "Number of reflection coefficients in f between tmin and tmax: "+str(ni)
  print "Random reflection coefficients: "+str(randomi)
  print "Reflection coefficients exist past "+\
  "tmax and corresponding u(tmax) in f and g, respecitively: "+str(moreps)
  print "Initial amount of squeezing: "+str(r0)
  print "Final amount of squeezing: "+str(r1)
  print "Shift between f and g: "+str(v)
  print "Wavelet c's peak frequency: "+str(freqc)
  print "Wavelet d's peak frequency: "+str(freqd)
  print "Wavelet c's decay: "+str(decayc)
  print "Wavelet d's decay: "+str(decayd)
  print "Wavelet c is minimum phase: "+str(mpc)
  print "Wavelet d is minimum phase: "+str(mpd)
  print "Noise to signal ratio in f and g: "+str(nrmsf)

  #Plotting
    #Data
  pngDir = None
  title= "f csbg p sbg bg q g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hmin1,hmin2 = -8.0,-8.0
  hmax1,hmax2 = 8.0,8.0
  hlabel = ["f","HSg","CSBg","CBSg","SCBg","g"]
  hminmax = [[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],\
  [hmin2,hmax2],[hmin2,hmax2],[hmin1,hmax1]]
  hint = [None,None,None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,csbg,cbsg,scbg,g],
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

    #Normalize
  ncguess = normalizeM(cguess)
  nbguess = normalizeM(bguess)
  ncw = normalizeM(cw)
  nbw = normalizeM(bw)
  nck = normalizeM(ck)
  ndk = normalizeM(dk)
  nak = normalizeM(ak)
  nbk = normalizeM(bk)

    #Wavelets
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known c"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known b"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbw,nbk],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  #A test to see if C and B commute with S with synthetics.
def testCommuteSino():
  x0 = 300
  halfwidth = 3
  f,g,u = getSinoTrace(x0)
  ref = RecursiveExponentialFilter(halfwidth)
  ref.apply1(g,g)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  alpha = 0.1
  penb = False
  penc = False
  pendeltab = False
  pendeltac = False
  niter = 5000
  ww = WaveletWarpingCBGN()
  ww.setMaxPercentChange(0.01)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setLineSearchMinScale(0.0000)
  ww.setStabilityFactor(0.0)
  ww.penalize(alpha,penb,penc,pendeltab,pendeltac)
  
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  #Get Gauss-Newton Information
  rmsri = ww.getRMSRi()
  rmsrf = ww.getRMSRf()
  twonormgrad = ww.getTwoNormGrad()
  steplength = ww.getStepLength()
  twonormdeltab = ww.getTwoNormDeltaB()
  twonormdeltac = ww.getTwoNormDeltaC()
  condnum = ww.getCondNum()
  twonormb = ww.getTwoNormB()
  twonormc = ww.getTwoNormC()
  rmspercentchange = ww.getRMSPercentChange()
  lastIter = ww.getLastIter()

  print "rms I love you Mom,Dad,Frank,Demi = "+str(ww.rms(ww.computeDataResidual(nc,kc,hw,nb,kb,bone,u,f,g)))
  

  #Processing
  warp = Warper()
  sg = warp.applyS(u,g)
  bsg = ww.applyC(nb,kb,bw,sg)
  cbsg = ww.applyC(nc,kc,cw,bsg)
  bg = ww.applyC(nb,kb,bw,g)
  cbg = ww.applyC(nc,kc,cw,bg)
  scbg = warp.applyS(u,cbg)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  #Print
  dt = 0.004
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)

  #Plotting
    #Data
  nt = len(f)
  pngDir = None
  title= "f hsg csbg cbsg scbg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[tmin*dt,tmax*dt],None
  hmin1,hmin2 = -8.0,-8.0
  hmax1,hmax2 = 8.0,8.0
  hlabel = ["f","HSg","CSBg","CBSg","SCBg"]
  hminmax = [[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],\
  [hmin2,hmax2],[hmin2,hmax2],[hmin1,hmax1]]
  hint = [None,None,None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,csbg,cbsg,scbg],
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

    #Normalize
  ncguess = normalizeM(cguess)
  nbguess = normalizeM(bguess)
  ncw = normalizeM(cw)
  nbw = normalizeM(bw)

    #Wavelets
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known c"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated and Known b"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

def chooseAlphaTestsLessSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 1081,20,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 1.6,1.4#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  freqd,decayd = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  #alphab = [1.0]
  #alphac = [0.0]
  alphab = [0.0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
  alphac = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  #alphab = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,0.0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  #alphab = [0.0,1e-6]
  #alphac = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,0.0,1e-5,1e-4,1e-3,1e-2,1e-1,1]
  #alphac = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  #alphac = [0.0,0.0]
  nalpha = len(alphab)
  condnums = zerofloat(nalpha)
  rmsrfs = zerofloat(nalpha)

  maxpc = 0.01
  niter = 1000
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)

  for i in range(0,nalpha):
    ww = WaveletWarpingCBGN()
    ww.setMaxPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
    ww.setLineSearchMinScale(0.0000)
    print "alphab = "+str(alphab[i])
    print "alphac = "+str(alphac[i])
    ww.setPenalize(alphab[i],alphac[i])
    
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    #Get Gauss-Newton Information
    rmsri = ww.getRMSRiTotal()
    rmsrf = ww.getRMSRfTotal()
    twonormgrad = ww.getTwoNormGrad()
    steplength = ww.getStepLength()
    twonormdeltab = ww.getTwoNormDeltaB()
    twonormdeltac = ww.getTwoNormDeltaC()
    condnum = ww.getCondNum()
    twonormb = ww.getTwoNormB()
    twonormc = ww.getTwoNormC()
    rmspercentchange = ww.getRMSPercentChange()
    lastIter = ww.getLastIter()
    #mf = MedianFinder(niter)
    #condnums[i] = mf.findMedian(condnum)
    #condnums[i] = sum(condnum)/lastIter
    #condnums[i] = condnum[lastIter]
    condnums[i] = max(condnum)
    rmsrfs[i] = rmsrf[lastIter]
    print "Condition Numbers:"
    dump(condnums)
    print "Final RMS of Residuals:"
    dump(rmsrfs)

    """
    #Individual plots
      #GN Measurements
    pngDir = None
    title = "alphab = "+str(alphab[i])+" RMS ri (Black) RMS rf (Blue)"
    maxrmsri = max(rmsri)
    minrmsrf = min(rmsrf)
    siter = Sampling(niter,1.0,0.0)
    color=[Color.BLACK,Color.BLUE]
    vlabel,vminmax,vint = "RMS of residuals",[minrmsrf,maxrmsri],None
    hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
    plotting.plotMeasInSamePlot(siter, [rmsri,rmsrf],\
    color=color,\
    vlabel=vlabel, vminmax=vminmax, vint=vint,\
    hlabel=hlabel, hminmax=hminmax, hint=hint,\
    title=title, pngDir=pngDir,\
    paper=True, onecol=True, twocol=None)

    title = "alphab = "+str(alphab[i])+" Step Length and Two Norm of Gradient"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,[steplength,twonormgrad],\
    color=None,\
    vlabel=["StepLength","2NormGrad"],\
    vminmax=[None,None],\
    vint=[None,None],\
    hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)

    title = "alphab = "+str(alphab[i])+" 2NormDeltaB,2NormDeltaC,2NormB,2NormC,ConditionNumber,RMS Percent Change"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,[twonormdeltab,twonormdeltac,\
    twonormb,twonormc,condnum,rmspercentchange],\
    color=None,\
    vlabel=["DeltMagb","DeltaMagc","2Normb","2Normc","CondNum","RMSPercChange"],\
    vminmax=[None,None,None,None,None,None],\
    vint=[None,None,None,None,None,None],\
    hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
    hsize=960,vsize=760,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)
    """

  #plotting
  print "condition numbers"
  dump(condnums)
  "rmsrfs"
  dump(rmsrfs)
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,max(rmsrfs))
  colors = [Color.GRAY,Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED,\
  Color.GRAY,Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED,\
  Color.GRAY,Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  mark = [PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,\
  PointsView.Mark.FILLED_SQUARE,PointsView.Mark.FILLED_SQUARE,PointsView.Mark.FILLED_SQUARE,PointsView.Mark.FILLED_SQUARE,PointsView.Mark.FILLED_SQUARE,PointsView.Mark.FILLED_SQUARE,PointsView.Mark.FILLED_SQUARE,PointsView.Mark.FILLED_SQUARE,\
  PointsView.Mark.HOLLOW_SQUARE,PointsView.Mark.HOLLOW_SQUARE,PointsView.Mark.HOLLOW_SQUARE,PointsView.Mark.HOLLOW_SQUARE,PointsView.Mark.HOLLOW_SQUARE,PointsView.Mark.HOLLOW_SQUARE,PointsView.Mark.HOLLOW_SQUARE,PointsView.Mark.HOLLOW_SQUARE,PointsView.Mark.HOLLOW_SQUARE]


  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(mark[i])

  #normalize condnums and rmsrfs to be able to calculate Euclidean distance correctly.
  """
  condnums = div(condnums,max(condnums))
  rmsrfs = div(rmsrfs,max(rmsrfs))
  dump(condnums)
  dump(rmsrfs)
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,1.0)
  colors = [Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  """

def chooseAlphaTestsMoreSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 1081,20,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  freqd,decayd = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  alphab = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
  alphac = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  nalpha = len(alphab)
  condnums = zerofloat(nalpha)
  rmsrfs = zerofloat(nalpha)

  maxpc = 0.00
  niter = 1000
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)

  for i in range(0,nalpha):
    ww = WaveletWarpingCBGN()
    ww.setMaxPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
    ww.setLineSearchMinScale(0.0000)
    ww.setPenalize(alphab[i],alphac[i])
    
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    #Get Gauss-Newton Information
    rmsri = ww.getRMSRiTotal()
    rmsrf = ww.getRMSRfTotal()
    condnum = ww.getCondNum()
    lastIter = ww.getLastIter()
    condnums[i] = condnum[lastIter]
    rmsrfs[i] = rmsrf[lastIter]
  #plotting
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,max(rmsrfs))
  colors = [Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  #normalize condnums and rmsrfs to be able to calculate Euclidean distance correctly.
  """
  condnums = div(condnums,max(condnums))
  rmsrfs = div(rmsrfs,max(rmsrfs))
  dump(condnums)
  dump(rmsrfs)
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,1.0)
  colors = [Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  """

def chooseAlphaTestsLessSqueezingNoise():
  #Synthetic parameters
  nt,ni,randomi = 1081,20,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 1.6,1.4#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  freqd,decayd = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  alphab = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
  alphac = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  nalpha = len(alphab)
  condnums = zerofloat(nalpha)
  rmsrfs = zerofloat(nalpha)

  maxpc = 0.00
  niter = 1000
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)

  for i in range(0,nalpha):
    ww = WaveletWarpingCBGN()
    ww.setMaxPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
    ww.setLineSearchMinScale(0.0000)
    ww.setPenalize(alphab[i],alphac[i])
    
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    #Get Gauss-Newton Information
    rmsri = ww.getRMSRiTotal()
    rmsrf = ww.getRMSRfTotal()
    condnum = ww.getCondNum()
    lastIter = ww.getLastIter()
    condnums[i] = condnum[lastIter]
    rmsrfs[i] = rmsrf[lastIter]
  #plotting
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,max(rmsrfs))
  colors = [Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  #normalize condnums and rmsrfs to be able to calculate Euclidean distance correctly.
  """
  condnums = div(condnums,max(condnums))
  rmsrfs = div(rmsrfs,max(rmsrfs))
  dump(condnums)
  dump(rmsrfs)
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,1.0)
  colors = [Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  """

def chooseAlphaTestsMoreSqueezingNoise():
  #Synthetic parameters
  nt,ni,randomi = 1081,20,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  freqd,decayd = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  alphab = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
  alphac = [0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  nalpha = len(alphab)
  condnums = zerofloat(nalpha)
  rmsrfs = zerofloat(nalpha)

  maxpc = 0.01
  niter = 1000
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)

  for i in range(0,nalpha):
    ww = WaveletWarpingCBGN()
    ww.setMaxPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
    ww.setLineSearchMinScale(0.0000)
    ww.setPenalize(alphab[i],alphac[i])
    
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    #Get Gauss-Newton Information
    rmsri = ww.getRMSRiTotal()
    rmsrf = ww.getRMSRfTotal()
    condnum = ww.getCondNum()
    lastIter = ww.getLastIter()
    condnums[i] = condnum[lastIter]
    rmsrfs[i] = rmsrf[lastIter]
  #plotting
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,max(rmsrfs))
  colors = [Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)

  #normalize condnums and rmsrfs to be able to calculate Euclidean distance correctly.
  """
  condnums = div(condnums,max(condnums))
  rmsrfs = div(rmsrfs,max(rmsrfs))
  dump(condnums)
  dump(rmsrfs)
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,1.0)
  colors = [Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
  """


 
def chooseAlphaTestsSino():
  #get sino images
  x0 = 250
  f,g,u = getSinoTrace(x0)
  st = Sampling(len(f),0.004,0.0)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 21,-10# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  alphab = [0.0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
  #alphac = [0.0,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
  alphac = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
  nalpha = len(alphab)
  condnums = zerofloat(nalpha)
  rmsrfs = zerofloat(nalpha)

  maxpc = 0.001
  niter = 1000
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)

  for i in range(0,nalpha):
    ww = WaveletWarpingCBGN()
    ww.setMaxPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
    ww.setLineSearchMinScale(0.0000)
    ww.setPenalize(alphab[i],alphac[i])
    
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    #Get Gauss-Newton Information
    #Get Gauss-Newton Information
    rmsri = ww.getRMSRiTotal()
    rmsrf = ww.getRMSRfTotal()
    twonormgrad = ww.getTwoNormGrad()
    steplength = ww.getStepLength()
    twonormdeltab = ww.getTwoNormDeltaB()
    twonormdeltac = ww.getTwoNormDeltaC()
    condnum = ww.getCondNum()
    twonormb = ww.getTwoNormB()
    twonormc = ww.getTwoNormC()
    rmspercentchange = ww.getRMSPercentChange()
    lastIter = ww.getLastIter()
    condnums[i] = condnum[lastIter]
    #condnums[i] = max(condnum)
    rmsrfs[i] = rmsrf[lastIter]

    #Individual plots
      #GN Measurements
    pngDir = None
    title = "alphab = "+str(alphab[i])+" RMS ri (Black) RMS rf (Blue)"
    maxrmsri = max(rmsri)
    minrmsrf = min(rmsrf)
    siter = Sampling(niter,1.0,0.0)
    color=[Color.BLACK,Color.BLUE]
    vlabel,vminmax,vint = "RMS of residuals",[minrmsrf,maxrmsri],None
    hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
    plotting.plotMeasInSamePlot(siter, [rmsri,rmsrf],\
    color=color,\
    vlabel=vlabel, vminmax=vminmax, vint=vint,\
    hlabel=hlabel, hminmax=hminmax, hint=hint,\
    title=title, pngDir=pngDir,\
    paper=True, onecol=True, twocol=None)

    title = "alphab = "+str(alphab[i])+" Step Length and Two Norm of Gradient"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,[steplength,twonormgrad],\
    color=None,\
    vlabel=["StepLength","2NormGrad"],\
    vminmax=[None,None],\
    vint=[None,None],\
    hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)

    title = "alphab = "+str(alphab[i])+" 2NormDeltaB,2NormDeltaC,2NormB,2NormC,ConditionNumber,RMS Percent Change"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,[twonormdeltab,twonormdeltac,\
    twonormb,twonormc,condnum,rmspercentchange],\
    color=None,\
    vlabel=["DeltMagb","DeltaMagc","2Normb","2Normc","CondNum","RMSPercChange"],\
    vminmax=[None,None,None,None,None,None],\
    vint=[None,None,None,None,None,None],\
    hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
    hsize=960,vsize=760,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)

  #overall plots
  print "condition numbers"
  dump(condnums)
  "rmsrfs"
  dump(rmsrfs)
  sp = SimplePlot()
  sp.setHLimits(0.0,max(condnums))
  sp.setVLimits(0.0,rmsri[0])
  colors = [Color.GRAY,Color.BLACK,Color.BLUE,Color.CYAN,Color.GREEN,Color.YELLOW,Color.ORANGE,Color.RED]
  for i in range(0,nalpha):
    pv = sp.addPoints([condnums[i]],[rmsrfs[i]])
    pv.setMarkColor(colors[i])
    pv.setLineStyle(PointsView.Line.NONE)
    pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)





def sinc(x):
  l = 1e-10
  if (-l<x and x<l):
    return 1.0
  else:
    return sin(x)/(x)
  
def getMinPhaseInverseWavelet(na,fpeak,decay):
  w = 2.0*PI*fpeak
  #print "w = "+str(w)
  r = exp(-decay)
  #print "r= "+str(r)
  a1,a2 = -2.0*r*cos(w),r*r
  a = zerofloat(na) 
  a[0] = 1
  a[1] = a1
  a[2] = a2
  for i in range(3,na):
    a[i] = 0.000000000000001
  return a



def getSimpleInverseWavelet(fpeak,decay):
  w = 2.0*PI*fpeak
  #print "w = "+str(w)
  r = exp(-decay)
  #print "r= "+str(r)
  a1,a2 = -2.0*r*cos(w),r*r
  a = [1,a1,a2]
  return a

def getInverseWavelet(fpeak,decay,nc,kc,na,ka,mp):
  if mp:
    print "na,ka = 3,0"
    print "nc,kc = "+str(nc)+","+str(kc)
    print "mp = "+str(mp)
    w = 2.0*PI*fpeak
    #print "w = "+str(w)
    r = exp(-decay)
    #print "r= "+str(r)
    a1,a2 = -2.0*r*cos(w),r*r
    a = [1,a1,a2]
    return a
  else:
    print "na,ka = "+str(na)+","+str(ka)
    print "nc,kc = "+str(nc)+","+str(kc)
    print "mp = "+str(mp)
    h = getWavelet(fpeak,decay,nc,kc,mp) # known wavelet
    a = ww.getWaveletC(nc,kc,hk,na,ka) # known inverse wavelet
    return a
    
def getWavelet(fpeak,decay,nc,kc,mp=False):
  x = zerofloat(nc)
  print "nc = "+str(nc)
  print "kc = "+str(kc)
  x[-kc] = 1.0
  return synthetic.addWavelet(fpeak,decay,x,mp)

def calcUPrime(dt, u):
  nt = len(u)
  up = zerofloat(nt)
  for i in range(1,nt):
    up[i] = (u[i]-u[i-1])/dt
  return up


def plotSequence(x,xmax=None,title=None):
  sp = SimplePlot.asPoints(x)
  if xmax==None:
    xmax = max(abs(max(x)),abs(min(x)))
    xmax *= 1.05
  sp.setVLimits(-xmax,xmax)
  if title:
    sp.setTitle(title)

def plotSequences(st,xs,amax=None,tmark=None,pngDir=None,labels=None,title=None):
  nx = len(xs)
  pp = PlotPanel(nx,1)
  for ix,xi in enumerate(xs):
    pv = pp.addPoints(ix,0,st,xi)
    if labels:
      pp.setVLabel(ix,labels[ix])
    if amax:
      pp.setVLimits(ix,-amax[ix]-1e-13,amax[ix]+1e-13)
    if tmark:
      pp.setVInterval(ix,tmark[ix])
  pp.setHLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setVisible(True)
  pf.setSize(800,900)
  if pngDir:
     pf.setFontSizeForPrint(8.0,222.0)
     pngDir = pngDir+title+"1wavelet"+"onecol.png"
     pf.paintToPng(720.0,3.08,pngDir)
  else:
    pp.setTitle(title)

def plotImages(st,sx,xs,tmin=None,tmax=None,amax=None,tmark=None,pngDir=None,labels=None,title=None):
  nx = len(xs)
  pp = PlotPanel(nx,1)
  for ix,xi in enumerate(xs):
    pv = pp.addPixels(0,ix,st,xs,xi)
    if amax:
      pv.setVLimits(ix,-amax[ix]-1e-13,amax[ix]+1e-13)
  pp.setHLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setVisible(True)
  pf.setSize(800,900)
  if tmin:
    pv.setHLimits(0,tmin,tmax)
  if pngDir:
     pf.setFontSizeForPrint(8.0,222.0)
     pngDir = pngDir+title+"1wavelet"+"onecol.png"
     pf.paintToPng(720.0,3.08,pngDir)
  else:
    pp.setTitle(title)


def plotSequenceMeas(st,xs,amax=None,amin=None,tmark=None,hlabel=None,pngDir=None,labels=None,title=None):
  nx = len(xs)
  pp = PlotPanel(1,1)
  color = [Color.BLACK,Color.BLUE]
  for ix,xi in enumerate(xs):
    pv = pp.addPoints(0,0,st,xi)
    pv.setLineColor(color[ix])
  if labels:
    pp.setVLabel(0,labels)
  if amax:
    pp.setVLimits(0,0,amax)
  if amin:
    pp.setVLimits(0,amin,amax)
  if tmark:
    pp.setVInterval(0,tmark)
  pp.setHLabel(hlabel)
  pf = PlotFrame(pp)
  pf.setVisible(True)
  pf.setSize(800,900)
  if pngDir:
     pf.setFontSizeForPrint(8.0,222.0)
     pngDir = pngDir+title+"Meas"+"onecol.png"
     pf.paintToPng(720.0,3.08,pngDir)
  else:
    pp.setTitle(title)


def plotWavelets(st,hs,hmax=None,title=None,pngDir=None,
  onecol=None,twocol=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  lw = 1.8,3,1
  nc = len(hs)
  hsmax = 1.0
  for ih in range(nc):
    if ih==0:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.SOLID)
      #pv.setLineStyle(PointsView.Line.NONE)
      pv.setLineColor(Color.BLUE)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(3.4)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
    if ih==1:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.SOLID)
      pv.setLineColor(Color.RED)
      pv.setLineWidth(1)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
    if ih==2:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.SOLID)
      #pv.setLineStyle(PointsView.Line.NONE)
      pv.setLineColor(Color.GREEN)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(3.4)
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

def plotUnknownWavelets(st,hs,hmax=None,title=None,pngDir=None,
  onecol=None,twocol=None):
  sp = SimplePlot()
  ls = [PointsView.Line.SOLID,PointsView.Line.DASH,PointsView.Line.DOT]
  lw = 1.8,3,1
  nc = len(hs)
  hsmax = 1.0
  for ih in range(nc):
    if ih==0:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.SOLID)
      #pv.setLineStyle(PointsView.Line.NONE)
      pv.setLineColor(Color.BLUE)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(3.4)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
    if ih==1:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.SOLID)
      #pv.setLineStyle(PointsView.Line.NONE)
      pv.setLineColor(Color.GREEN)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(3.4)
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

def getTitle(shift, s, constwarp, r, varywarp, r1, r2):
  if shift:
    title = "shift s = "+str(s)
  elif constwarp:
    title = "constant warp r = "+str(r)
  elif varywarp:
    title = "varying warp r2 = "+str(r2)+" r1 = "+str(r1)
  print title
  return title


#Returns the signal added with noise and the noise that 
#was added.
def addNoise(nrms, seed, f):
  n = len(f)
  r = Random(seed)
  g = mul(2.0,sub(randfloat(r,n),0.5))
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  #print "nrms = "+str(nrms)
  #print "rmsnoise/rmssignal = "+str(rms(g)/rms(f))
  return add(f,g),g

#Returns the signal added with noise and the noise that 
#was added.
def addNoise2D(nrms, seed, f):
  nx = len(f)
  nt = len(f[0])
  noise = zerofloat(nt,nx)
  fnoise = zerofloat(nt,nx)
  r = RandomFloat(50)
  for ix in range(0,nx):
    #nrmsx = r.uniform()*nrms#random noise to signal ratio between 0 and nrms
    nrmsx = nrms
    fnoise[ix],noise[ix] = addNoise(nrmsx,seed,f[ix])
    #print "nrmsx = "+str(nrmsx)
  return fnoise,noise


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

def normalizeMAAWOS(f):
  minf = min(f)
  print "minf = "+str(minf)
  maxf = max(f)
  absminf = abs(minf) 
  absmaxf = abs(maxf) 
  if absminf<absmaxf:
    return mul(1.0/maxf,f)
  else:
    return mul(1.0/minf,f)

def normalizeM(f):
  nf = len(f)
  sums = 0.0
  maxv = max(f)
  print "maxv = "+str(maxv)
  f = div(f,maxv)
  return f


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
  amp = div(amp,max(amp))
  plotSpectrum(fs,amp,title,amax=amax)

#Plot phase spectrum for a single trace
def plotPhaseSpectrumT(st, p, itmin, itmax, title, amax=None):
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
  phase = computePhaseSpectrum(subst,fs,nfft,subp)
  plotPhaseSpectrum(fs,phase,title,amax=amax)


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
  wft = rampfloat(0.0,-2.0*FLT_PI*df*ft,nf)
  cf = cmul(cf,cmplx(cos(wft),sin(wft)))

  af = cabs(cf)
  #Amplitude spectrum normalized
  #float amax = max(max(af),FLT_EPSILON)
  #af = mul(1.0f/amax,af)
  return af

def computePhaseSpectrum(st, fs, nfft, p):
  nt = st.getCount()
  dt = st.getDelta()
  ft = st.getFirst()
  nf = fs.getCount()
  df = fs.getDelta()
  ff = fs.getFirst()

  #Real-to-complex fast Fourier transform.
  fft = FftReal(nfft)
  cf = zerofloat(2*nf)
  copy(nt,p,cf)
  fft.realToComplex(-1,cf,cf)

  #Adjust phase for possibly non-zero time of first sample.
  wft = rampfloat(0.0,-2.0*PI*df*ft,nf)
  cf = cmul(cf,cmplx(cos(wft),sin(wft)))

  #Phase response, in cycles.
  pf = carg(cf)
  pf = mul(0.5/PI,pf)
  return pf

#Computes the backward difference of a sequence.
def computeBackDiff(u):
  nu = len(u)
  bd = zerofloat(nu)
  for i in range(1,nu):
    bd[i] = (u[i] - u[i-1])/1.0
  return bd 
def computeBackDiff2D(u):
  nu2 = len(u)
  nu1 = len(u[0])
  bd = zerofloat(nu1,nu2)
  for i in range(0,nu2):
    bd[i] = computeBackDiff(u[i])
  return bd 

def addZeros(x,nt):
  nx = len(x)
  y = zerofloat(nt)
  if nx==nt:
    return x
  else:
    for i in range(0,nx):
      y[i] = x[i]
    for i in range(nx,nt):
      y[i] = 0.0
    return y
 
def dot(f,g):
  nf = len(f)
  sum = 0.0
  for i in range(nf):
    sum = sum + f[i]*g[i]
  return sum


def searchFor(x,n):
  i = 0
  while x[i]<=n:
    i = i + 1
  return i

def plotMatrix(x,title=None):
  s = SimplePlot()
  pv = s.addPixels(x)
  pv.setColorModel(ColorMap.BLUE_WHITE_RED)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setClips(-max(x),max(x))
  nx = len(x)
  print "corner = "+str(x[0][nx-2])
  s.addColorBar()
  if title:
    s.addTitle(title)

def plotImage(st,sx,f,tmin=None,tmax=None,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(350,450)
  pv = sp.addPixels(st,sx,f)
  if tmin:
    sp.setVLimits(tmin,tmax)
  #pv.setClips(-2.0,2.0)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("PP time (s)")
  if title:
    sp.setTitle(title)
 
def plot4ImagesSideBySide(st, sx, f, g, h, i, tmin, tmax, clip, vlabel, title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):
  pv1 = PixelsView(st,sx,f)
  pv1.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv1.setClips(-clip,clip)
  pv2 = PixelsView(st,sx,g)
  pv2.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv2.setClips(-clip,clip)
  pv3 = PixelsView(st,sx,h)
  pv3.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv3.setClips(-clip,clip)
  pv4 = PixelsView(st,sx,i)
  pv4.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
  pv4.setClips(-clip,clip)
  
  pp = PlotPanel(1,4,PlotPanel.Orientation.X1DOWN_X2RIGHT,
  PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(0,1,pv2)
  pp.addTiledView(0,2,pv3)
  pp.addTiledView(0,3,pv4)
  pp.setVLimits(tmin,tmax)
  pp.setVLabel(vlabel)
  pf = PlotFrame(pp)
  pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      print "zz = "+str(aspectRatio)
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(8.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)

  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)


def plotSpectrum(sf,spec,title,amax=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setVLabel("Amplitude")
  sp.setHLabel("Frequency (Hz)")
  sp.setSize(750,400)
  sp.addTitle(title)
  if amax:
    sp.setVLimits(0,amax)
  pv = sp.addPoints(sf,spec)
  #sp.setVLimits(0,0.054)
  #sp.setHLimits(19,21)

def plotPhaseSpectrum(sf,spec,title,amax=None):
  sp = SimplePlot(SimplePlot.Origin.LOWER_LEFT)
  sp.setVLabel("Cycles")
  sp.setHLabel("Frequency (Hz)")
  sp.setSize(750,400)
  sp.addTitle(title)
  if amax:
    sp.setVLimits(0,amax)
  pv = sp.addPoints(sf,spec)
  #sp.setVLimits(0,0.054)
  #sp.setHLimits(19,21)

def getSinoTrace(x0):
  dataDir = "/Users/Chris/data/sinos/"
  n1f,n1g,d1,f1 = 501,852,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1f,n2)
  g = readImage(dataDir+"ps.dat",n1g,n2)
  u = readImage(dataDir+"shifts.dat",n1f,n2)
  u = add(u,rampfloat(0.0,1.0,0.0,n1f,n2))
  fr = zerofloat(n1f)
  gr = zerofloat(n1g)
  ur = zerofloat(n1f)
  fr = f[x0]
  gr = g[x0]
  ur = u[x0]
  gain(100,fr)
  gain(100,gr)
  return fr,gr,ur

def getSinoImage(x0,nx):
  dataDir = "/Users/Chris/data/sinos/"
  n1f,n1g,d1,f1 = 501,852,0.004,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1f,n2)
  g = readImage(dataDir+"ps.dat",n1g,n2)
  u = readImage(dataDir+"shifts.dat",n1f,n2)
  u = add(u,rampfloat(0.0,1.0,0.0,n1f,n2))
  fr = zerofloat(n1f,nx)
  gr = zerofloat(n1g,nx)
  ur = zerofloat(n1f,nx)
  for ix in range(0,nx):
    fr[ix] = f[x0+ix]
    gr[ix] = g[x0+ix]
    ur[ix] = u[x0+ix]
    print "x0+ix = "+str(x0+ix)
  gain(100,fr)
  gain(100,gr)
  return fr,gr,ur

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

#Creates a subset of an array that contains the sum of square differences for an
#iterative algorithm.
#This subset will contain elements of the original array that are within a certain
#range of the last sum of square differences value. This method walks back from the last
#element of the array and determines if the next element is within this range. 
#If the next element
#is not within this range and all the following elements are within the range, then the next and
#following elements will not exist in the array subset.
def zoomConverge(ssd,width):
  nssd = len(ssd)
  lastv = ssd[nssd-1]
  count = 0
  for i in range(nssd):
    if ((lastv+width/2.0>ssd[i]) and (lastv-width/2.0<ssd[i])):
      count = count + 1
  nsubssd = count
  subssd = zerofloat(nsubssd)
  fiter = 0
  for i in range(nsubssd):
    fiter = nssd-1-i
    subssd[nsubssd-1-i] = ssd[fiter]
  return subssd,fiter

def getMedian(f):
  nf = len(f)
  i = rampint(0,1,nf)
  quickIndexSort(f,i)
  if (nf%2==0):
    median = (f[nf/2-1]+f[nf/2])/2
  else:
    median = f[nf/2]
  return median

  


  
  
  

def getWarpingPartsSpectrums(r,itmin,itmax,f,g,u):
  ww = WaveletWarpingCBGN()
  nu = len(u)
  nuu = int(r*(nu-1)+1)
  ng = len(g)
  nug = int(r*(ng-1)+1)
  nf = len(f)
  dtgu = 1.0/r
  w = 0.20/r
  warp = Warper()
  itminu = int(r*(itmin-1)+1)
  itmaxu = int(r*(itmax-1)+1)

  uu = warp.upSampleLinear(nu,1.0,0.0,u,nuu,dtgu,0.0)
  ug = warp.upSample(ng,1.0,0.0,g,nug,dtgu,0.0)
  wug  = warp.applyW(dtgu,uu,ug)
  lwug = warp.applyL(r,uu,wug)
  dlwug = warp.subSample(r,nu,lwug)
  stf = Sampling(nf,1.0,0.0)
  stg = Sampling(ng,1.0,0.0)
  stug = Sampling(nug,dtgu,0.0)
  stuu = Sampling(nuu,dtgu,0.0)
  stsug = Sampling(nuu,dtgu,0.0)
  maxf = max(f)
  maxfd2 = maxf/2.0
  amax = [maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  print len(f)
  print nf
  #plotSequences(Sampling(len(f),1.0,0.0),[f,dlwug],amax=amax,tmark=tmark,\
  #labels=["f","dlwug"],title="fdlwug")
  #plotSequences(stg,[g],amax=amax,tmark=tmark,\
  #labels=["g"],title="g")
  #amax = [maxf]
  #tmark = [maxfd2]
  #plotSequences(stug,[ug],amax=amax,tmark=tmark,\
  #labels=["ug"],title="ug")
  #plotSequences(stuu,[wug],amax=amax,tmark=tmark,\
  #labels=["sug"],title="sug")
  plotAmplitudeSpectrumT(stf, dlwug, itmin, itmax, "amp DLWUg", amax=None)
  #plotAmplitudeSpectrumT(stsug, lwug, itminu, itmaxu, "amp LWUg", amax=None)
  #plotAmplitudeSpectrumT(stsug, wug, itminu, itmaxu, "amp WUg", amax=None)
  plotAmplitudeSpectrumT(stug, ug, 0, nug, "amp Ug", amax=None)
  plotAmplitudeSpectrumT(stg, g, 0, ng, "amp g", amax=None)
  plotAmplitudeSpectrumT(stf, f, itmin, itmax, "amp f", amax=None)

def usePreviousbw(nb,bw):
  newbw = zerofloat(nb)
  for i in range(0,nb-1):
    newbw[i] = bw[i]
  return newbw

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

