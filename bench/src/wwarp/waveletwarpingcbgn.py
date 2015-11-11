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
  #displaySynthetics()
  #showPowerOfLargeShapingFilterWithConstantSqueezing()
  #showPowerOfLargeShapingFilterWithVaryingSqueezing()
  #showPowerOfReasonableShapingFilterWithConstantSqueezing()
  #showPowerOfReasonableShapingFilterWithVaryingSqueezing()
  #displaySimpleCaseOfWarpingWithWavelets()
  #addingAZeroToBAndCDoesNotMakeADifference()
  #showScaleFactorAmbiguity() 
  #showShiftAmibiguity()
  #showUniqueness()
  #showCaseOfNotFindingWaveletsWithConstantSqueezing()
  #whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedB1pt6To1pt4()
  #whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedB3pt15To1pt55()
  #whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedB1pt6To1pt4Noise()
  #whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedB3pt15To1pt55Noise()
  #whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedBSino()
  #whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedBSinoNoiseReduction()
  #choosePenalizationSynthetic3pt15To1pt55()
  #choosePenalizationSynthetic1pt6To1pt4()
  #syntheticNoiseReduction()
  #sinoNoiseReduction()
  #goSino1000ShapingFilter()
  #goTestNcNbTradeOff()
  #goTestNcNbLessSqueezeNoiseTradeOff()
  #goTestTradeOffMinLessSqueezeNoise()
  #goTestTradeOffMinNoNoise()
  #goTestTradeOffMinNoise()
  #goTestNcNb2DTradeOff()
  #goTestNcNb2DNoiseReductionTradeOff()
  #goTestNcNb2DTradeOffMin()
  #goTestSinopecNoiseReduction()
  #goTestNbNhIncrease()
  #goTest1SinoNbNhIncrease()
  #goTestLessSqueezingNbNhIncrease()
  #goTestNbNh2DIncrease()
  #goTestNbNhAll2DIncrease()
  #goSinopec1()
  #goSinopec()
  #goAccumulateVsAntiAliasingDeltaFunctions()
  #testWarpingandBandPassFilter()

#Simply displays the synthetic of your choosing for quick viewing purposes.
def displaySynthetics():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
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
  pngDir = None
  title= "Synthetic traces p, q, f, and g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1,hint2 = 0.5,5.0
  hmin1,hmin2 = -1.0,-8.0
  hmax1,hmax2 = 1.0,8.0
  hlabel = ["p","q","f","g"]
  hminmax = [[hmin1,hmax1],[hmin1,hmax1],[hmin2,hmax2],[hmin2,hmax2]]
  #hint = [hint1,hint1,hint2,hint2]
  hint = [None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[p,q,f,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

#Shows that a large shaping filter (as large as the time window) can 
#account for time-invariant distortion.
def showPowerOfLargeShapingFilterWithConstantSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 2.0,2.0#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
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

  #Shaping filter estimation parameters
  nh,kh = 500,-250
  hstabfact = 0.000001
  nb,kb = 1,0
  b = [1]

  #Processing
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  warper = Warper()
  hw = ww.getWaveletC(nh,kh,nb,kb,b,hstabfact,u,f,g)
  sg = warper.applyS(u,g)
  hsg = ww.applyH(nh,kh,hw,sg)

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

  #Normalization
  nhw = normalizeM(hw)

  #Plotting
  pngDir = None
  title= "Shaping filter and constant squeezing"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1 = 5.0
  hmin1 = -8.0
  hmax1 = 8.0
  hlabel = ["f","HSg","Sg","g"]
  hminmax = [[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1]]
  hint = [None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,sg,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  title = "h"
  st = Sampling(nh,dt,kh*dt)
  hint = 0.5
  plotting.plotWavelets(st,[nhw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

#Shows that a large shaping filter (as large as the time window) can 
#account for time-variant distortion.
def showPowerOfLargeShapingFilterWithVaryingSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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

  #Shaping filter estimation parameters
  nh,kh = 500,-250
  hstabfact = 0.000001
  nb,kb = 1,0
  b = [1]

  #Processing
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  warper = Warper()
  hw = ww.getWaveletC(nh,kh,nb,kb,b,hstabfact,u,f,g)
  sg = warper.applyS(u,g)
  hsg = ww.applyH(nh,kh,hw,sg)

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

  #Normalization
  nhw = normalizeM(hw)

  #Plotting
  pngDir = None
  title= "Shaping filter and constant squeezing"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1 = 5.0
  hmin1 = -8.0
  hmax1 = 8.0
  hlabel = ["f","HSg","Sg","g"]
  hminmax = [[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1]]
  hint = [None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,sg,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  title = "h"
  hint = 0.5
  st = Sampling(nh,dt,kh*dt)
  plotting.plotWavelets(st,[nhw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

#Shows that a reasonable shaping filter (as large as twice Sg's dominant period (10), so twice
#is 20) can account for most of the time-invariant distortion. Dominant period calculation is at .549-.508s 
#in Sg.
def showPowerOfReasonableShapingFilterWithConstantSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 2.0,2.0#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
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

  #Shaping filter estimation parameters
  nh,kh = 81,-40
  hstabfact = 0.0
  nb,kb = 1,0
  b = [1]

  #Processing
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  warper = Warper()
  hw = ww.getWaveletC(nh,kh,nb,kb,b,hstabfact,u,f,g)
  sg = warper.applyS(u,g)
  hsg = ww.applyH(nh,kh,hw,sg)

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

  #Normalization
  nhw = normalizeM(hw)

  #Plotting
  pngDir = None
  title= "Shaping filter and constant squeezing"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1 = 5.0
  hmin1 = -8.0
  hmax1 = 8.0
  hlabel = ["f","HSg","Sg","g"]
  hminmax = [[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1]]
  hint = [None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,sg,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  title = "h"
  st = Sampling(nh,dt,kh*dt)
  hint = 0.5
  plotting.plotWavelets(st,[nhw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

#Shows that a reasonable shaping filter (as large as twice Sg's dominant period (10), so twice
#is 20) can account for only some of the time-invariant distortion. 
#Dominant period calculation is at .549-.508s in Sg.
def showPowerOfReasonableShapingFilterWithVaryingSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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

  #Shaping filter estimation parameters
  nh,kh = 81,-40
  hstabfact = 0.0
  nb,kb = 1,0
  b = [1]

  #Processing
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  warper = Warper()
  hw = ww.getWaveletC(nh,kh,nb,kb,b,hstabfact,u,f,g)
  sg = warper.applyS(u,g)
  hsg = ww.applyH(nh,kh,hw,sg)

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

  #Normalization
  nhw = normalizeM(hw)

  #Plotting
  pngDir = None
  title= "Shaping filter and constant squeezing"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1 = 5.0
  hmin1 = -8.0
  hmax1 = 8.0
  hlabel = ["f","HSg","Sg","g"]
  hminmax = [[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1]]
  hint = [None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,sg,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  title = "h"
  hint = 0.5
  st = Sampling(nh,dt,kh*dt)
  plotting.plotWavelets(st,[nhw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

#Displays the initial reflectivity sequence, the wavelet used, the synthetic seismogram,
#wavelet distortion, and squeezed impulses.
def displaySimpleCaseOfWarpingWithWavelets():
#Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
  r0,r1 = 2.0,2.0#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  freqd,decayd = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1DSimple(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,randomi,moreps)

  #Get know wavelets
  ww = WaveletWarpingCBGN()
  nb,kb = 5,-2
  nc,kc = 181,-90
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Processing
  warp = Warper()
  warp.setScale(False)
  snoscaleg = warp.applyS(u,g)

  warp = Warper()
  sq = warp.applyS(u,q)
  bg = ww.applyC(nb,kb,bk,g)
  sg = warp.applyS(u,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,ck,sbg)
  cbsg = ww.applyC(nc,kc,ck,ww.applyC(nb,kb,bk,sg))
  scbg = warp.applyS(u,ww.applyC(nc,kc,ck,bg))
  SimplePlot.asPoints(cbsg)
  SimplePlot.asPoints(scbg)
  q0 = copy(sbg)
  q1 = ww.delay(1,q0)
  q2 = ww.delay(2,q0)
  p0 = ww.applyC(nc,kc,ck,warp.applyS(u,ww.delay(0,g)))
  p1 = ww.applyC(nc,kc,ck,warp.applyS(u,ww.delay(1,g)))
  p2 = ww.applyC(nc,kc,ck,warp.applyS(u,ww.delay(2,g)))



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
  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces p, q"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[p,q],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces p, Sq"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[p,sq],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  """"
  #Wavelet interpolation
  dt = 0.004
  scale = 8
  nc2 = scale*nc
  dt2 = 0.004/scale
  si = SincInterpolator.fromErrorAndFrequency(0.00001,0.45)
  cki = zerofloat(nc2)
  si.interpolate(nc,dt,-kc*dt,ck,nc2,dt2,-kc*dt,cki)
  """

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title = "Known Wavelet (c)"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  linecolor=[Color.BLACK,Color.BLACK,Color.BLACK]
  markstyle=[PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  plotting.plotWavelets(st,[ck],hint=hint,linecolor=linecolor,markstyle=markstyle,\
  title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
 
  """
  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title = "Interpolated Known Wavelet (c)"
  hint = None
  st = Sampling(nc2,dt2,kc*dt)
  linecolor=[Color.BLACK,Color.BLACK,Color.BLACK]
  markstyle=[PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  plotting.plotWavelets(st,[cki],hint=hint,linecolor=linecolor,markstyle=markstyle,\
  title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  """

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces f, g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces f, Sg No scale"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,snoscaleg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces f, Sg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,sg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces f, Bg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,bg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces f, SBg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,sbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces f, CSBg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces f, SCBg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,scbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title = "Synthetic traces f, CBSg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 5
  hmin = -12
  hmax = 12
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,cbsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces q0,q1,q2"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 1.0
  hmin = -1.5
  hmax = 1.5
  hlabel = ["Amplitude","Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[q0,q1,q2],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  
  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces q0"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,250*dt
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint = 1.0
  hmin = -1.5
  hmax = 1.5
  hlabel = ["Amplitude"]
  hminmax = [[hmin,hmax]]
  hint = [hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[q0],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces q1"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,250*dt
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint = 1.0
  hmin = -1.5
  hmax = 1.5
  hlabel = ["Amplitude"]
  hminmax = [[hmin,hmax]]
  hint = [hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[q1],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces q2"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,250*dt
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint = 1.0
  hmin = -1.5
  hmax = 1.5
  hlabel = ["Amplitude"]
  hminmax = [[hmin,hmax]]
  hint = [hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[q2],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)



  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces p0,p1,p2"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint = 150
  hmin = -250
  hmax = 250
  hlabel = ["Amplitude","Amplitude","Amplitude"]
  hminmax = [[hmin,hmax],[hmin,hmax],[hmin,hmax]]
  hint = [hint,hint,hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[p0,p1,p2],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces p0"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,250*dt
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint = 150
  hmin = -250
  hmax = 250
  hlabel = ["Amplitude"]
  hminmax = [[hmin,hmax]]
  hint = [hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[p0],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces p1"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,250*dt
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint = 150
  hmin = -250
  hmax = 250
  hlabel = ["Amplitude"]
  hminmax = [[hmin,hmax]]
  hint = [hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[p1],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = "./slides/simpleWaveletsandWarping/"
  #pngDir = None
  title= "Synthetic traces p2"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,250*dt
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint = 150
  hmin = -250
  hmax = 250
  hlabel = ["Amplitude"]
  hminmax = [[hmin,hmax]]
  hint = [hint]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[p2],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)













##Learned that you can add zeros to the end of b or c and not change the overall
#RMS of the residuals if you do not change kb or kc. 
def addingAZeroToBAndCDoesNotMakeADifference():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,False# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
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

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax
  
  nc,kc = 81,-40
  nb,kb = 4,-2

  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)

  #Processing
  bg = ww.applyH(nb,kb,b,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,c,sbg)

  bgzero = ww.applyH(nbz,kbz,bzero,g)
  sbgzero = warp.applyS(u,bgzero)
  csbgzero = ww.applyH(ncz,kcz,czero,sbgzero)


  rmsrcb = ww.rms(ww.computeDataResidual(nc,kc,c,nb,kb,b,u,f,g))
  rmsrcbzero = ww.rms(ww.computeDataResidual(ncz,kcz,czero,nbz,kbz,bzero,u,f,g))
  print "rmsrcb = "+str(rmsrcb)
  print "rmsrcbzero = "+str(rmsrcbzero)
  
  """
  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  #pngDir = "./tests/"
  pngDir = None
  title = "Increaseby1"+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max(rmsrcb),max(rmsrh)])
  minrms = min([min(rmsrcb),min(rmsrh)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(rmsrh))
  print "rmsrcb = "+str(len(rmsrcb))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+20",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcb],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  
  title = "Increaseby1"+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant Condition Number"
  maxcn = max(cn)
  mincn = min(cn)
  print "min = "+str(mincn)
  print "max = "+str(maxcn)
  print "n = "+str(n)
  si = Sampling(n,1.0,nc)
  color=None
  vlabel,vminmax,vint = "Condition Number",[mincn,maxcn],None
  hlabel,hminmax,hint = "nh=nb+20",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [cn],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  """
  
def showScaleFactorAmbiguity():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
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

  #Get know wavelets
  ww = WaveletWarpingCBGN()
  nb,kb = 5,-2
  nc,kc = 81,-40
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #No Scale Processing
  warp = Warper()
  bg = ww.applyH(nb,kb,bk,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,ck,sbg)

  #Scale Processing
  bkscale10 = mul(bk,10.0)
  ckscale10 = mul(ck,0.1)
  warp = Warper()
  bgscale10 = ww.applyH(nb,kb,bk,g)
  sbgscale10 = warp.applyS(u,bgscale10)
  csbgscale10 = ww.applyH(nc,kc,ck,sbgscale10)

  #Scale Processing
  bkscale100 = mul(bk,100.0)
  ckscale100 = mul(ck,0.01)
  warp = Warper()
  bgscale100 = ww.applyH(nb,kb,bk,g)
  sbgscale100 = warp.applyS(u,bgscale100)
  csbgscale100 = ww.applyH(nc,kc,ck,sbgscale100)

  #Plotting
    #Data
  pngDir = None
  title= "csbg csbgscale10 csbgscale100"
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-8.0
  hmax1,hmax2 = 1.0,8.0
  hlabel = ["csbg","csbgscale10","csbgscale100"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[csbg,csbgscale10,csbgscale100],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

def showShiftAmibiguity():
  #Synthetic parameters
  nt,ni,randomi = 1081,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
  freqc,decayc = 0.08,0.05
  freqd,decayd = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf

  #Get know wavelets
  ww = WaveletWarpingCBGN()
  nb,kb = 5,-2
  nc,kc = 81,-40
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  r0,r1 = 1.6,1.4#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #No shift Processing
  warp = Warper()
  bg = ww.applyH(nb,kb,bk,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,ck,sbg)

  #Shift Processing
  bkshift = ww.delay(1,bk)
  ckshift = ww.delay(1,ck)

  #shift Processing
  warp = Warper()
  bgshift = ww.applyH(nb,kb,bkshift,g)
  sbgshift = warp.applyS(u,bgshift)
  csbgshift = ww.applyH(nc,kc,ckshift,sbgshift)

  #Residual computation
  csbgmf = sub(csbg,f)
  csbgshiftmf = sub(csbgshift,f)
  print "rms csbgmf = "+str(rms(csbgmf))
  print "rms csbgmfshift = "+str(rms(csbgshiftmf))

  #Plotting
    #Data
  pngDir = None
  title= "csbg csbgshift"
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -200.0,-8.0
  hmax1,hmax2 = 200.0,8.0
  hlabel = ["f","csbg","csbgshift","csbg-f","csbgshift-f"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,csbg,csbgshift,csbgmf,csbgshiftmf],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

def showUniqueness():
  #Synthetic parameters
  nt,ni,randomi = 1081,10,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  v = 0.0#The amount of shift between p and q.
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  freqc,decayc = [0.08,0.16],[0.07,0.07]
  freqd,decayd = [0.08,0.08],[0.07,0.07]
  r0,r1 = [3.15,1.6,2.0],[1.55,1.4,2.0]#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  nrmsf = [0.0,0.2,0.5]
  nrmsg = nrmsf

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 21,-10# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb
  alphab = 0.00
  alphac = 0.00
  maxpc = 0.0001
  niter = 500

  nw = len(freqc)
  nnoise = len(nrmsf)
  nr = len(r0)
  ntot = nw*nnoise*nr
  count = 0

  freqcs = zerofloat(ntot)
  freqds = zerofloat(ntot)
  r0s = zerofloat(ntot)
  r1s = zerofloat(ntot)
  nrmsfs = zerofloat(ntot)
  firstcondnums = zerofloat(ntot)
  lastcondnums = zerofloat(ntot)
  maxcondnums = zerofloat(ntot)
  diffcondnums = zerofloat(ntot)
  percchangecondnums = zerofloat(ntot)
  totalchangecondnums = zerofloat(ntot)
  twonormbs = zerofloat(ntot)
  twonormcs = zerofloat(ntot)

  for iw in range(0,nw):
    for inoise in range(0,nnoise):
      for ir in range(0,nr):
        #Create synthetic f and g.
        p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc[iw],\
        decayc[iw],mpc,freqd[iw],decayd[iw],mpd,r0[ir],r1[ir],v,nrmsf[inoise],nrmsg[inoise],\
        nt,ni,randomi,moreps)

        #set tmin and tmax 
        tmin = tmin
        tmax = tmax

        ww = WaveletWarpingCBGN()
        ww.setMaxPercentChange(maxpc)#units are percentage.
        ww.setTimeRange(tmin,tmax)
        ww.setLineSearchMinScale(0.0000)
        ww.setPenalize(alphab,alphac)

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
        lastiter = ww.getLastIter()
        freqcs[count] = freqc[iw]
        freqds[count] = freqd[iw]
        r0s[count] = r0[ir]
        r1s[count] = r1[ir]
        nrmsfs[count] = nrmsf[inoise]
        firstcondnums[count] = condnum[0]
        lastcondnums[count] = condnum[lastiter]
        maxcondnums[count] = max(condnum)
        diffcondnums[count] = condnum[lastiter]-condnum[lastiter-1]
        percchangecondnums[count] = (diffcondnums[count]/condnum[lastiter-1])*100
        totalchangecondnums[count] = ((condnum[lastiter] - condnum[0])/condnum[0])*100
        twonormbs[count] = twonormb[lastiter]
        twonormcs[count] = twonormc[lastiter]
        count = count+1
  print "freqc"
  dump(freqcs)
  print "freqd"
  dump(freqds)
  print "r0s"
  dump(r0s)
  print "r1s"
  dump(r1s)
  print "nrmsfs"
  dump(nrmsfs)
  print"last change condnums"
  dump(diffcondnums)
  print "1st condnums"
  dump(firstcondnums)
  print "last condnums"
  dump(lastcondnums)
  print "max condnums"
  dump(maxcondnums)
  print "percent change"
  dump(percchangecondnums)
  print "total change"
  dump(totalchangecondnums)
  print "twonormbs"
  dump(twonormbs)
  print "twonormcs"
  dump(twonormcs)

def showCaseOfNotFindingWaveletsWithConstantSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 2.0,2.0#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0
  nrmsg = 0.0

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Wavelet estimation parameters
  nb,kb = 4,0#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  alphab = 0.0
  alphac = 0.0
  maxpc = 0.00
  niter = 1000
  ww = WaveletWarpingCBGN()
  ww.setMaxPercentChange(maxpc)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setLineSearchMinScale(0.0000)
  ww.setPenalize(alphab,alphac)
  
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  #bguess = [1.000000,  0.233888]
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
  lastiter = ww.getLastIter()

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
  title = "RMS ri (Black) RMS rf (Blue)"
  maxrmsri = max(rmsri)
  minrmsrf = min(rmsrf)
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastiter],None
  plotting.plotMeasInSamePlot(siter, [rmsri,rmsrf],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)

  title = "Step Length and Two Norm of Gradient"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,[steplength,twonormgrad],\
  color=None,\
  vlabel=["StepLength","2NormGrad"],\
  vminmax=[None,None],\
  vint=[None,None],\
  hlabel="Iterations",hminmax=[0,lastiter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  title = "2NormDeltaB,2NormDeltaC,2NormB,2NormC,ConditionNumber,RMS Percent Change"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,[twonormdeltab,twonormdeltac,\
  twonormb,twonormc,condnum,rmspercentchange],\
  color=None,\
  vlabel=["DeltMagb","DeltaMagc","2Normb","2Normc","CondNum","RMSPercChange"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0,lastiter],hint=None,\
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




def whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedB1pt6To1pt4():
#Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  maxpc = 0.0
  niter = 500
  nhfinal = 88
  nhinitial = 81
  khinitial = -40
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcbcon = zerofloat(n)#Starts with the previous b and c solved for.
  rmsrcbso = zerofloat(n)#Starts with an impulse and a shaping filter.
  cncon= zerofloat(n)
  cnso= zerofloat(n)
  wwh = WaveletWarpingCBGN()
  wwh.setTimeRange(tmin,tmax)
  wwh.setMaxPercentChange(maxpc)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
  h1 = wwh.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cwcon = h1
  bwcon = bone
  for i in range(0,n):
    print "################ "+str(i)+" out of "+str(n)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    alpha = 0.00
    sfac = 0.00001
    penb,penc,simpleStab = False,False,True
    warp = Warper()
      #Continued and Start Over
    wwcon = WaveletWarpingCBGN()
    wwso = WaveletWarpingCBGN()
    wwcon.setTimeRange(tmin,tmax)
    wwso.setTimeRange(tmin,tmax)
    wwcon.setMaxPercentChange(maxpc)
    wwso.setMaxPercentChange(maxpc)
    wwcon.setPenalize(alpha,sfac,penb,penc,simpleStab)
    wwso.setPenalize(alpha,sfac,penb,penc,simpleStab)
    #Shaping filter
    bone = zerofloat(nb)
    bone[-kb] = 1.0
    print "nc = "+str(nc)
    print "nb = "+str(nb)
    h = wwcon.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)

    #First guesses of c and b. 
      #Continued and Start Over
    cwcon = addZeros(cwcon,nc)
    cwso = wwso.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    bwcon = addZeros(bwcon,nb)
    bwso = bone
    print "start c guess Continued"
    dump(cwcon)
    print "start c guess Start Over"
    dump(cwso)
    print "start b guess Continued"
    dump(bwcon)
    print "start b guess Start over"
    dump(bwso)
    print "Continued:"
    cbwcon = wwcon.getWaveletCInverseB(nb,kb,bwcon,nc,kc,cwcon,u,f,g,niter)
    print "Start Over:"
    cbwso = wwso.getWaveletCInverseB(nb,kb,bwso,nc,kc,cwso,u,f,g,niter)
    cwcon = cbwcon[0]
    cwso = cbwso[0]
    bwcon = cbwcon[1]
    bwso = cbwso[1]
    print "final c Continued"
    dump(cwcon)
    print "final c Start Over"
    dump(cwso)
    print "final b Continued"
    dump(bwcon)
    print "final b Start over"
    dump(bwso)

    #Get iteration information
    rmsrh[i] = wwh.rms(wwh.computeDataResidual(nh,kh,h,nb,kb,bone,u,f,g))
    print "rmsrhAAA = "+str(rmsrh[i])
      #Continued and Start Over
    lastitercon = wwcon.getLastIter()
    lastiterso = wwso.getLastIter()
    rmsrfcon = wwcon.getRMSRf()
    rmsrfso = wwso.getRMSRf()
    condnumcon = wwcon.getCondNum()
    condnumso = wwso.getCondNum()
    cncon[i] = condnumcon[lastitercon]
    cnso[i] = condnumso[lastiterso]
    rmsrcbcon[i] = rmsrfcon[lastitercon]
    rmsrcbso[i] = rmsrfso[lastiterso]
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb Continued = "+str(rmsrcbcon)
    print "rmsrcb Start Over = "+str(rmsrcbso)
    print "Final Condition Number Continued = "+str(condnumcon[lastitercon])
    print "Final Condition Number Continued = "+str(condnumso[lastiterso])
    nb = nb+1
    nh = nh+1
    kh = -nh/2+1

  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  #pngDir = "./increaseby1_500iter/"
  pngDir = None
  title = "Increaseby1_1pt6To1pt4 "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max([rmsrcbcon,rmsrcbso]),max(rmsrh)])
  minrms = min([min([rmsrcbcon,rmsrcbso]),min(rmsrh)])
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcbcon,rmsrcbso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  
  title = "Increaseby1_1pt6To1pt4CN "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant Condition Number"
  maxcn = max([cncon,cnso])
  mincn = min([cncon,cnso])
  si = Sampling(n,1.0,nc)
  color=[Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "Condition Number",[mincn,maxcn],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [cncon,cnso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)


def whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedB3pt15To1pt55():
#Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  maxpc = 0.0
  niter = 500
  sfac = 0.00001
  nhfinal = 88
  nhinitial = 81
  khinitial = -40
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcbcon = zerofloat(n)#Starts with the previous b and c solved for.
  rmsrcbso = zerofloat(n)#Starts with an impulse and a shaping filter.
  cncon= zerofloat(n)
  cnso= zerofloat(n)
  wwh = WaveletWarpingCBGN()
  wwh.setTimeRange(tmin,tmax)
  wwh.setMaxPercentChange(maxpc)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
  h1 = wwh.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cwcon = h1
  bwcon = bone
  for i in range(0,n):
    print "################ "+str(i)+" out of "+str(n)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    alpha = 0.00
    penb,penc,pendeltab,pendeltac = False,False,False,False
    warp = Warper()
      #Continued and Start Over
    wwcon = WaveletWarpingCBGN()
    wwso = WaveletWarpingCBGN()
    wwcon.setTimeRange(tmin,tmax)
    wwso.setTimeRange(tmin,tmax)
    wwcon.setMaxPercentChange(maxpc)
    wwso.setMaxPercentChange(maxpc)
    wwcon.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    wwso.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    #Shaping filter
    bone = zerofloat(nb)
    bone[-kb] = 1.0
    print "nc = "+str(nc)
    print "nb = "+str(nb)
    h = wwcon.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)

    #First guesses of c and b. 
      #Continued and Start Over
    cwcon = addZeros(cwcon,nc)
    cwso = wwso.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    bwcon = addZeros(bwcon,nb)
    bwso = bone
    print "start c guess Continued"
    dump(cwcon)
    print "start c guess Start Over"
    dump(cwso)
    print "start b guess Continued"
    dump(bwcon)
    print "start b guess Start over"
    dump(bwso)
    cbwcon = wwcon.getWaveletCInverseB(nb,kb,bwcon,nc,kc,cwcon,u,f,g,niter)
    cbwso = wwso.getWaveletCInverseB(nb,kb,bwso,nc,kc,cwso,u,f,g,niter)
    cwcon = cbwcon[0]
    cwso = cbwso[0]
    bwcon = cbwcon[1]
    bwso = cbwso[1]
    print "final c Continued"
    dump(cwcon)
    print "final c Start Over"
    dump(cwso)
    print "final b Continued"
    dump(bwcon)
    print "final b Start over"
    dump(bwso)

    #Get iteration information
    rmsrh[i] = wwh.rms(wwh.computeDataResidual(nh,kh,h,nb,kb,bone,u,f,g))
    print "rmsrhAAA = "+str(rmsrh[i])
      #Continued and Start Over
    lastitercon = wwcon.getLastIter()
    lastiterso = wwso.getLastIter()
    rmsrfcon = wwcon.getRMSRf()
    rmsrfso = wwso.getRMSRf()
    condnumcon = wwcon.getCondNum()
    condnumso = wwso.getCondNum()
    cncon[i] = condnumcon[lastitercon]
    cnso[i] = condnumso[lastiterso]
    rmsrcbcon[i] = rmsrfcon[lastitercon]
    rmsrcbso[i] = rmsrfso[lastiterso]
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb Continued = "+str(rmsrcbcon)
    print "rmsrcb Start Over = "+str(rmsrcbso)
    print "Final Condition Number Continued = "+str(condnumcon[lastitercon])
    print "Final Condition Number Continued = "+str(condnumso[lastiterso])
    nb = nb+1
    nh = nh+1
    kh = -nh/2+1

  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  pngDir = "./increaseby1_500iter/"
  #pngDir = None
  title = "Increaseby1_3pt15To1pt55 "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max([rmsrcbcon,rmsrcbso]),max(rmsrh)])
  minrms = min([min([rmsrcbcon,rmsrcbso]),min(rmsrh)])
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcbcon,rmsrcbso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  
  title = "Increaseby1_3pt15To1pt55CN "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant Condition Number"
  maxcn = max([cncon,cnso])
  mincn = min([cncon,cnso])
  si = Sampling(n,1.0,nc)
  color=[Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "Condition Number",[mincn,maxcn],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [cncon,cnso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)

def whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedB1pt6To1pt4Noise():
#Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  maxpc = 0.0
  niter = 500
  sfac = 0.00001
  nhfinal = 88
  nhinitial = 81
  khinitial = -40
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcbcon = zerofloat(n)#Starts with the previous b and c solved for.
  rmsrcbso = zerofloat(n)#Starts with an impulse and a shaping filter.
  cncon= zerofloat(n)
  cnso= zerofloat(n)
  wwh = WaveletWarpingCBGN()
  wwh.setTimeRange(tmin,tmax)
  wwh.setMaxPercentChange(maxpc)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
  h1 = wwh.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cwcon = h1
  bwcon = bone
  for i in range(0,n):
    print "################ "+str(i)+" out of "+str(n)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    alpha = 0.00
    penb,penc,pendeltab,pendeltac = False,False,False,False
    warp = Warper()
      #Continued and Start Over
    wwcon = WaveletWarpingCBGN()
    wwso = WaveletWarpingCBGN()
    wwcon.setTimeRange(tmin,tmax)
    wwso.setTimeRange(tmin,tmax)
    wwcon.setMaxPercentChange(maxpc)
    wwso.setMaxPercentChange(maxpc)
    wwcon.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    wwso.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    #Shaping filter
    bone = zerofloat(nb)
    bone[-kb] = 1.0
    print "nc = "+str(nc)
    print "nb = "+str(nb)
    h = wwcon.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)

    #First guesses of c and b. 
      #Continued and Start Over
    cwcon = addZeros(cwcon,nc)
    cwso = wwso.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    bwcon = addZeros(bwcon,nb)
    bwso = bone
    print "start c guess Continued"
    dump(cwcon)
    print "start c guess Start Over"
    dump(cwso)
    print "start b guess Continued"
    dump(bwcon)
    print "start b guess Start over"
    dump(bwso)
    print "Continued:"
    cbwcon = wwcon.getWaveletCInverseB(nb,kb,bwcon,nc,kc,cwcon,u,f,g,niter)
    print "Start Over:"
    cbwso = wwso.getWaveletCInverseB(nb,kb,bwso,nc,kc,cwso,u,f,g,niter)
    cwcon = cbwcon[0]
    cwso = cbwso[0]
    bwcon = cbwcon[1]
    bwso = cbwso[1]
    print "final c Continued"
    dump(cwcon)
    print "final c Start Over"
    dump(cwso)
    print "final b Continued"
    dump(bwcon)
    print "final b Start over"
    dump(bwso)

    #Get iteration information
    rmsrh[i] = wwh.rms(wwh.computeDataResidual(nh,kh,h,nb,kb,bone,u,f,g))
    print "rmsrhAAA = "+str(rmsrh[i])
      #Continued and Start Over
    lastitercon = wwcon.getLastIter()
    lastiterso = wwso.getLastIter()
    rmsrfcon = wwcon.getRMSRf()
    rmsrfso = wwso.getRMSRf()
    condnumcon = wwcon.getCondNum()
    condnumso = wwso.getCondNum()
    cncon[i] = condnumcon[lastitercon]
    cnso[i] = condnumso[lastiterso]
    rmsrcbcon[i] = rmsrfcon[lastitercon]
    rmsrcbso[i] = rmsrfso[lastiterso]
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb Continued = "+str(rmsrcbcon)
    print "rmsrcb Start Over = "+str(rmsrcbso)
    print "Final Condition Number Continued = "+str(condnumcon[lastitercon])
    print "Final Condition Number Continued = "+str(condnumso[lastiterso])
    nb = nb+1
    nh = nh+1
    kh = -nh/2+1

  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  pngDir = "./increaseby1_500iter/"
  #pngDir = None
  title = "Increaseby1_1pt6To1pt4Noise "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max([rmsrcbcon,rmsrcbso]),max(rmsrh)])
  minrms = min([min([rmsrcbcon,rmsrcbso]),min(rmsrh)])
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcbcon,rmsrcbso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  
  title = "Increaseby1_1pt6To1pt4NoiseCN "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant Condition Number"
  maxcn = max([cncon,cnso])
  mincn = min([cncon,cnso])
  si = Sampling(n,1.0,nc)
  color=[Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "Condition Number",[mincn,maxcn],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [cncon,cnso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)


def whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedB3pt15To1pt55Noise():
#Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  maxpc = 0.0
  niter = 500
  sfac = 0.00001
  nhfinal = 88
  nhinitial = 81
  khinitial = -40
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcbcon = zerofloat(n)#Starts with the previous b and c solved for.
  rmsrcbso = zerofloat(n)#Starts with an impulse and a shaping filter.
  cncon= zerofloat(n)
  cnso= zerofloat(n)
  wwh = WaveletWarpingCBGN()
  wwh.setTimeRange(tmin,tmax)
  wwh.setMaxPercentChange(maxpc)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
  h1 = wwh.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cwcon = h1
  bwcon = bone
  for i in range(0,n):
    print "################ "+str(i)+" out of "+str(n)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    alpha = 0.00
    penb,penc,pendeltab,pendeltac = False,False,False,False
    warp = Warper()
      #Continued and Start Over
    wwcon = WaveletWarpingCBGN()
    wwso = WaveletWarpingCBGN()
    wwcon.setTimeRange(tmin,tmax)
    wwso.setTimeRange(tmin,tmax)
    wwcon.setMaxPercentChange(maxpc)
    wwso.setMaxPercentChange(maxpc)
    wwcon.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    wwso.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    #Shaping filter
    bone = zerofloat(nb)
    bone[-kb] = 1.0
    print "nc = "+str(nc)
    print "nb = "+str(nb)
    h = wwcon.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)

    #First guesses of c and b. 
      #Continued and Start Over
    cwcon = addZeros(cwcon,nc)
    cwso = wwso.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    bwcon = addZeros(bwcon,nb)
    bwso = bone
    print "start c guess Continued"
    dump(cwcon)
    print "start c guess Start Over"
    dump(cwso)
    print "start b guess Continued"
    dump(bwcon)
    print "start b guess Start over"
    dump(bwso)
    cbwcon = wwcon.getWaveletCInverseB(nb,kb,bwcon,nc,kc,cwcon,u,f,g,niter)
    cbwso = wwso.getWaveletCInverseB(nb,kb,bwso,nc,kc,cwso,u,f,g,niter)
    cwcon = cbwcon[0]
    cwso = cbwso[0]
    bwcon = cbwcon[1]
    bwso = cbwso[1]
    print "final c Continued"
    dump(cwcon)
    print "final c Start Over"
    dump(cwso)
    print "final b Continued"
    dump(bwcon)
    print "final b Start over"
    dump(bwso)

    #Get iteration information
    rmsrh[i] = wwh.rms(wwh.computeDataResidual(nh,kh,h,nb,kb,bone,u,f,g))
    print "rmsrhAAA = "+str(rmsrh[i])
      #Continued and Start Over
    lastitercon = wwcon.getLastIter()
    lastiterso = wwso.getLastIter()
    rmsrfcon = wwcon.getRMSRf()
    rmsrfso = wwso.getRMSRf()
    condnumcon = wwcon.getCondNum()
    condnumso = wwso.getCondNum()
    cncon[i] = condnumcon[lastitercon]
    cnso[i] = condnumso[lastiterso]
    rmsrcbcon[i] = rmsrfcon[lastitercon]
    rmsrcbso[i] = rmsrfso[lastiterso]
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb Continued = "+str(rmsrcbcon)
    print "rmsrcb Start Over = "+str(rmsrcbso)
    print "Final Condition Number Continued = "+str(condnumcon[lastitercon])
    print "Final Condition Number Continued = "+str(condnumso[lastiterso])
    nb = nb+1
    nh = nh+1
    kh = -nh/2+1

  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  pngDir = "./increaseby1_500iter/"
  #pngDir = None
  title = "Increaseby1_3pt15To1pt55Noise "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max([rmsrcbcon,rmsrcbso]),max(rmsrh)])
  minrms = min([min([rmsrcbcon,rmsrcbso]),min(rmsrh)])
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcbcon,rmsrcbso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  
  title = "Increaseby1_3pt15To1pt55NoiseCN "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant Condition Number"
  maxcn = max([cncon,cnso])
  mincn = min([cncon,cnso])
  si = Sampling(n,1.0,nc)
  color=[Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "Condition Number",[mincn,maxcn],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [cncon,cnso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)

def whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedBSino():
#Synthetic parameters
  #get sino trace
  x0 = 300
  f,g,u = getSinoTrace(x0)

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  maxpc = 0.0
  niter = 500
  sfac = 0.00001
  nhfinal = 30
  nhinitial = 21
  khinitial = -10
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcbcon = zerofloat(n)#Starts with the previous b and c solved for.
  rmsrcbso = zerofloat(n)#Starts with an impulse and a shaping filter.
  cncon= zerofloat(n)
  cnso= zerofloat(n)
  wwh = WaveletWarpingCBGN()
  wwh.setTimeRange(tmin,tmax)
  wwh.setMaxPercentChange(maxpc)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
  h1 = wwh.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cwcon = h1
  bwcon = bone
  for i in range(0,n):
    print "################ "+str(i)+" out of "+str(n)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    alpha = 0.00
    penb,penc,pendeltab,pendeltac = False,False,False,False
    warp = Warper()
      #Continued and Start Over
    wwcon = WaveletWarpingCBGN()
    wwso = WaveletWarpingCBGN()
    wwcon.setTimeRange(tmin,tmax)
    wwso.setTimeRange(tmin,tmax)
    wwcon.setMaxPercentChange(maxpc)
    wwso.setMaxPercentChange(maxpc)
    wwcon.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    wwso.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    #Shaping filter
    bone = zerofloat(nb)
    bone[-kb] = 1.0
    print "nc = "+str(nc)
    print "nb = "+str(nb)
    h = wwcon.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)

    #First guesses of c and b. 
      #Continued and Start Over
    cwcon = addZeros(cwcon,nc)
    cwso = wwso.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    bwcon = addZeros(bwcon,nb)
    bwso = bone
    print "start c guess Continued"
    dump(cwcon)
    print "start c guess Start Over"
    dump(cwso)
    print "start b guess Continued"
    dump(bwcon)
    print "start b guess Start over"
    dump(bwso)
    cbwcon = wwcon.getWaveletCInverseB(nb,kb,bwcon,nc,kc,cwcon,u,f,g,niter)
    cbwso = wwso.getWaveletCInverseB(nb,kb,bwso,nc,kc,cwso,u,f,g,niter)
    cwcon = cbwcon[0]
    cwso = cbwso[0]
    bwcon = cbwcon[1]
    bwso = cbwso[1]
    print "final c Continued"
    dump(cwcon)
    print "final c Start Over"
    dump(cwso)
    print "final b Continued"
    dump(bwcon)
    print "final b Start over"
    dump(bwso)

    #Get iteration information
    rmsrh[i] = wwh.rms(wwh.computeDataResidual(nh,kh,h,nb,kb,bone,u,f,g))
    print "rmsrhAAA = "+str(rmsrh[i])
      #Continued and Start Over
    lastitercon = wwcon.getLastIter()
    lastiterso = wwso.getLastIter()
    rmsrfcon = wwcon.getRMSRf()
    rmsrfso = wwso.getRMSRf()
    condnumcon = wwcon.getCondNum()
    condnumso = wwso.getCondNum()
    cncon[i] = condnumcon[lastitercon]
    cnso[i] = condnumso[lastiterso]
    rmsrcbcon[i] = rmsrfcon[lastitercon]
    rmsrcbso[i] = rmsrfso[lastiterso]
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb Continued = "+str(rmsrcbcon)
    print "rmsrcb Start Over = "+str(rmsrcbso)
    print "Final Condition Number Continued = "+str(condnumcon[lastitercon])
    print "Final Condition Number Continued = "+str(condnumso[lastiterso])
    nb = nb+1
    nh = nh+1
    kh = -nh/2+1

  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  pngDir = "./increaseby1_500iter/"
  #pngDir = None
  title = "Increaseby1Sino "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max([rmsrcbcon,rmsrcbso]),max(rmsrh)])
  minrms = min([min([rmsrcbcon,rmsrcbso]),min(rmsrh)])
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcbcon,rmsrcbso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  
  title = "Increaseby1SinoCN "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant Condition Number"
  maxcn = max([cncon,cnso])
  mincn = min([cncon,cnso])
  si = Sampling(n,1.0,nc)
  color=[Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "Condition Number",[mincn,maxcn],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [cncon,cnso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)

def whatDoesAddingBDoToPreviousReasonableShapingFilterDoWithThePreviousEstimatedBSinoNoiseReduction():
#Synthetic parameters
  #get sino trace
  x0 = 300
  halfwidth = 3
  f,g,u = getSinoTrace(x0)
  ref = RecursiveExponentialFilter(halfwidth)
  ref.apply1(g,g)


  #set tmin and tmax 
  tmin = 100
  tmax = 500

  maxpc = 0.0
  niter = 500
  sfac = 0.00001
  nhfinal = 30
  nhinitial = 21
  khinitial = -10
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcbcon = zerofloat(n)#Starts with the previous b and c solved for.
  rmsrcbso = zerofloat(n)#Starts with an impulse and a shaping filter.
  cncon= zerofloat(n)
  cnso= zerofloat(n)
  wwh = WaveletWarpingCBGN()
  wwh.setTimeRange(tmin,tmax)
  wwh.setMaxPercentChange(maxpc)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
  h1 = wwh.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cwcon = h1
  bwcon = bone
  for i in range(0,n):
    print "################ "+str(i)+" out of "+str(n)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    alpha = 0.00
    penb,penc,pendeltab,pendeltac = False,False,False,False
    warp = Warper()
      #Continued and Start Over
    wwcon = WaveletWarpingCBGN()
    wwso = WaveletWarpingCBGN()
    wwcon.setTimeRange(tmin,tmax)
    wwso.setTimeRange(tmin,tmax)
    wwcon.setMaxPercentChange(maxpc)
    wwso.setMaxPercentChange(maxpc)
    wwcon.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    wwso.setPenalize(alpha,penb,penc,pendeltab,pendeltac)
    #Shaping filter
    bone = zerofloat(nb)
    bone[-kb] = 1.0
    print "nc = "+str(nc)
    print "nb = "+str(nb)
    h = wwcon.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)

    #First guesses of c and b. 
      #Continued and Start Over
    cwcon = addZeros(cwcon,nc)
    cwso = wwso.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    bwcon = addZeros(bwcon,nb)
    bwso = bone
    print "start c guess Continued"
    dump(cwcon)
    print "start c guess Start Over"
    dump(cwso)
    print "start b guess Continued"
    dump(bwcon)
    print "start b guess Start over"
    dump(bwso)
    cbwcon = wwcon.getWaveletCInverseB(nb,kb,bwcon,nc,kc,cwcon,u,f,g,niter)
    cbwso = wwso.getWaveletCInverseB(nb,kb,bwso,nc,kc,cwso,u,f,g,niter)
    cwcon = cbwcon[0]
    cwso = cbwso[0]
    bwcon = cbwcon[1]
    bwso = cbwso[1]
    print "final c Continued"
    dump(cwcon)
    print "final c Start Over"
    dump(cwso)
    print "final b Continued"
    dump(bwcon)
    print "final b Start over"
    dump(bwso)

    #Get iteration information
    rmsrh[i] = wwh.rms(wwh.computeDataResidual(nh,kh,h,nb,kb,bone,u,f,g))
    print "rmsrhAAA = "+str(rmsrh[i])
      #Continued and Start Over
    lastitercon = wwcon.getLastIter()
    lastiterso = wwso.getLastIter()
    rmsrfcon = wwcon.getRMSRf()
    rmsrfso = wwso.getRMSRf()
    condnumcon = wwcon.getCondNum()
    condnumso = wwso.getCondNum()
    cncon[i] = condnumcon[lastitercon]
    cnso[i] = condnumso[lastiterso]
    rmsrcbcon[i] = rmsrfcon[lastitercon]
    rmsrcbso[i] = rmsrfso[lastiterso]
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb Continued = "+str(rmsrcbcon)
    print "rmsrcb Start Over = "+str(rmsrcbso)
    print "Final Condition Number Continued = "+str(condnumcon[lastitercon])
    print "Final Condition Number Continued = "+str(condnumso[lastiterso])
    nb = nb+1
    nh = nh+1
    kh = -nh/2+1

  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  pngDir = "./increaseby1_500iter/"
  #pngDir = None
  title = "NRIncreaseby1SinoNR "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max([rmsrcbcon,rmsrcbso]),max(rmsrh)])
  minrms = min([min([rmsrcbcon,rmsrcbso]),min(rmsrh)])
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcbcon,rmsrcbso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  
  title = "NRIncreaseby1SinoNRCN "+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant Condition Number"
  maxcn = max([cncon,cnso])
  mincn = min([cncon,cnso])
  si = Sampling(n,1.0,nc)
  color=[Color.BLUE,Color.RED]
  vlabel,vminmax,vint = "Condition Number",[mincn,maxcn],None
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],None
  plotting.plotMeasInSamePlot(si, [cncon,cnso],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)





def choosePenalizationSynthetic3pt15To1pt55():
#Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  maxpc = 0.01
  niter = 1000
  nhfinal = 31
  nhinitial = 21
  khinitial = -10

  sfac = [0.0001,0.0,0.0,0.0,0.0,0.0,0.0]
  #sfac = [0.0,0.0,0.0,0.0,0.0,0.0]
  penb = [False,True,False,True,False,False,False]
  penc = [False,False,True,True,False,False,False]
  pendeltab = [False,False,False,False,True,False,True]
  pendeltac = [False,False,False,False,False,True,True]
  #penb = [False,False,False]
  #penc = [False,False,False]
  #pendeltab = [True,False,True]
  #pendeltac = [False,True,True]
  alpha = [0.000001,0.00001,0.0001,0.001,0.01]

  np = len(penb)
  na = len(alpha)
  for ip in range(0,np):
    for ia in range(0,na):
      print "ip = "+str(ip)
      print "ia = "+str(ia)
      nb = 1
      kb = 0
      nc = nhinitial
      kc = khinitial
      nh = nhinitial
      kh = khinitial
      n = (nhfinal-nh)+1
      rmsrh = zerofloat(n)
      rmsrcb = zerofloat(n)
      cn = zerofloat(n)
      for i in range(0,n):
        print "################ "+str(i)+" out of "+str(n)+" #####################"
        print "nb = "+str(nb)
        print "kb = "+str(kb)
        print "nc = "+str(nc)
        print "kc = "+str(kc)
        print "nh = "+str(nh)
        print "kh = "+str(kh)
        #Estimate wavelet
        warp = Warper()
        ww = WaveletWarpingCBGN()
        ww.setTimeRange(tmin,tmax)
        ww.setMaxPercentChange(maxpc)
        ww.setPenalize(alpha[ia],penb[ip],penc[ip],pendeltab[ip],pendeltac[ip])
        #Shaping filter
        hstabfact = 0.0
        bone = zerofloat(nb)
        bone[-kb] = 1.0
        h = ww.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)

        #First guesses of c and b. 
        print "nb = "+str(nb)
        print "kb = "+str(kb)
        print "nc = "+str(nc)
        print "kc = "+str(kc)
        c = ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
        cbw = ww.getWaveletCInverseB(nb,kb,bone,nc,kc,c,u,f,g,niter)
        cw = cbw[0]
        bw = cbw[1]

        #Get iteration information
        lastiter = ww.getLastIter()
        rmsrf = ww.getRMSRf()
        condnum = ww.getCondNum()
        cn[i] = condnum[lastiter]
        rmsrcb[i] = rmsrf[lastiter]
        rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,bone,u,f,g))
        print "rmsrh = "+str(rmsrh)
        print "rmsrcb = "+str(rmsrcb)
        print "Final Condition Number = "+str(condnum[lastiter])
        nb = nb+1
        kb = -nb/2+1
        nh = nh+1
        kh = -nh/2+1

      #Plotting
      #plot rms of residuals with shaping filter (rmsrh) 
      #and rms of residuals with c and b (rmsrcb)
      #############################################################
      pngDir = "./SyntheticPenalizeTests3pt15To1pt55/"
      #pngDir = None 
      title = "ip = "+str(ip)+" alpha "+str(alpha[ia])+" penb "+str(penb[ip])+" penc "+str(penc[ip])+" pendeltab "+str(pendeltab[ip])+" pendeltac "+str(pendeltac[ip])+" Increaseby1"+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
      maxrms = max([max(rmsrcb),max(rmsrh)])
      minrms = min([min(rmsrcb),min(rmsrh)])
      print "min = "+str(minrms)
      print "max = "+str(maxrms)
      print "n = "+str(n)
      print "rmsrh = "+str(len(rmsrh))
      print "rmsrcb = "+str(len(rmsrcb))
      si = Sampling(n,1.0,nc)
      color=[Color.BLACK,Color.BLUE]
      vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
      hlabel,hminmax,hint = "nh=nb+20",[nhinitial,nhfinal],None
      plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcb],\
      color=color,\
      vlabel=vlabel, vminmax=vminmax, vint=vint,\
      hlabel=hlabel, hminmax=hminmax, hint=hint,\
      title=title, pngDir=pngDir,\
      paper=True, onecol=True, twocol=None)
      
      #############################################################

def choosePenalizationSynthetic1pt6To1pt4():
#Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  maxpc = 0.01
  niter = 1000
  nhfinal = 31
  nhinitial = 21
  khinitial = -10

  sfac = [0.0001,0.0,0.0,0.0,0.0,0.0,0.0]
  #sfac = [0.0,0.0,0.0,0.0,0.0,0.0]
  penb = [False,True,False,True,False,False,False]
  penc = [False,False,True,True,False,False,False]
  pendeltab = [False,False,False,False,True,False,True]
  pendeltac = [False,False,False,False,False,True,True]
  #penb = [False,False,False]
  #penc = [False,False,False]
  #pendeltab = [True,False,True]
  #pendeltac = [False,True,True]
  alpha = [0.00001,0.0001,0.001,0.01]

  np = len(penb)
  na = len(alpha)
  for ip in range(0,np):
    for ia in range(0,na):
      print "ip = "+str(ip)
      print "ia = "+str(ia)
      nb = 1
      kb = 0
      nc = nhinitial
      kc = khinitial
      nh = nhinitial
      kh = khinitial
      n = (nhfinal-nh)+1
      rmsrh = zerofloat(n)
      rmsrcb = zerofloat(n)
      cn = zerofloat(n)
      for i in range(0,n):
        print "################ "+str(i)+" out of "+str(n)+" #####################"
        print "nb = "+str(nb)
        print "kb = "+str(kb)
        print "nc = "+str(nc)
        print "kc = "+str(kc)
        print "nh = "+str(nh)
        print "kh = "+str(kh)
        #Estimate wavelet
        warp = Warper()
        ww = WaveletWarpingCBGN()
        ww.setTimeRange(tmin,tmax)
        ww.setMaxPercentChange(maxpc)
        ww.setPenalize(alpha[ia],penb[ip],penc[ip],pendeltab[ip],pendeltac[ip])
        #Shaping filter
        hstabfact = 0.0
        bone = zerofloat(nb)
        bone[-kb] = 1.0
        h = ww.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)

        #First guesses of c and b. 
        print "nb = "+str(nb)
        print "kb = "+str(kb)
        print "nc = "+str(nc)
        print "kc = "+str(kc)
        c = ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
        cbw = ww.getWaveletCInverseB(nb,kb,bone,nc,kc,c,u,f,g,niter)
        cw = cbw[0]
        bw = cbw[1]

        #Get iteration information
        lastiter = ww.getLastIter()
        rmsrf = ww.getRMSRf()
        condnum = ww.getCondNum()
        cn[i] = condnum[lastiter]
        rmsrcb[i] = rmsrf[lastiter]
        rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,bone,u,f,g))
        print "rmsrh = "+str(rmsrh)
        print "rmsrcb = "+str(rmsrcb)
        print "Final Condition Number = "+str(condnum[lastiter])
        nb = nb+1
        kb = -nb/2+1
        nh = nh+1
        kh = -nh/2+1

      #Plotting
      #plot rms of residuals with shaping filter (rmsrh) 
      #and rms of residuals with c and b (rmsrcb)
      #############################################################
      pngDir = "./SyntheticPenalizeTests1pt6To1pt4/"
      #pngDir = None 
      title = "ip = "+str(ip)+" alpha "+str(alpha[ia])+" penb "+str(penb[ip])+" penc "+str(penc[ip])+" pendeltab "+str(pendeltab[ip])+" pendeltac "+str(pendeltac[ip])+" Increaseby1"+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
      maxrms = max([max(rmsrcb),max(rmsrh)])
      minrms = min([min(rmsrcb),min(rmsrh)])
      print "min = "+str(minrms)
      print "max = "+str(maxrms)
      print "n = "+str(n)
      print "rmsrh = "+str(len(rmsrh))
      print "rmsrcb = "+str(len(rmsrcb))
      si = Sampling(n,1.0,nc)
      color=[Color.BLACK,Color.BLUE]
      vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
      hlabel,hminmax,hint = "nh=nb+20",[nhinitial,nhfinal],None
      plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcb],\
      color=color,\
      vlabel=vlabel, vminmax=vminmax, vint=vint,\
      hlabel=hlabel, hminmax=hminmax, hint=hint,\
      title=title, pngDir=pngDir,\
      paper=True, onecol=True, twocol=None)
      
      #############################################################



def syntheticNoiseReduction():
  #Synthetic parameters
  nt,ni,randomi = 1081,30,True# number of time samples in p and q; number of random impulses in p and q.
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
  p,q,f,g,noiseinf,noiseing,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #set tmin and tmax 
  dt = 0.004

  halfwidth = [1,3,5,10]
  n = len(halfwidth)
  for i in range(0,n):

    ref = RecursiveExponentialFilter(halfwidth[i])
    nrg = zerofloat(nt)
    ref.apply1(g,nrg)

    #Make RMS of image equal to 1
    ########################################################
    ww = WaveletWarpingCBGN()
    f = ww.makerms1(tmin,tmax,f)
    nonoiseg = ww.makerms1(tmin,tmax,sub(g,noiseing))
    noiseg = ww.makerms1(tmin,tmax,g)
    nrg = ww.makerms1(tmin,tmax,nrg)
    
    #plotting
    pngDir = None
    title= "g halfwidth = "+str(halfwidth[i])
    st = Sampling(nt,dt,0.0)
    vmin,vmax = 0*dt,nt*dt
    vlabel,vminmax,vint = "Time (s)",None,None
    hint1 = 2.5
    hmin1 = -5.0
    hmax1 = 5.0
    hlabel = ["g","No noise g","Noise Reduced g"]
    hminmax = [[hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1]]
    hint = [None,None,None]
    hsize,vsize = 960,560
    tilespacing = 5
    plotting.plotTracesSideBySide(st,[noiseg,nonoiseg,nrg],\
    vlabel=vlabel,vminmax=vminmax,vint=vint,\
    hlabel=hlabel,hminmax=hminmax,hint=hint,\
    tilespacing = tilespacing,\
    title=title,pngDir=pngDir,\
    paper=True,onecol=True)
    

def sinoNoiseReduction():
  #get sino trace
  x0 = 300
  f,g,u = getSinoTrace(x0)
  ng = len(g)
  nu = len(u)

  #set tmin and tmax 
  dt = 0.004
  itminf = 100
  itmaxf = 500
  itming = 225
  itmaxg = 800
  tminf = itminf*dt
  tmaxf = itmaxf*dt
  tming = itming*dt
  tmaxg = itmaxg*dt
  ntg = len(g)
  ntf = len(f)

  halfwidth = [3,5,10]
  n = len(halfwidth)
  for i in range(0,n):
    ref = RecursiveExponentialFilter(halfwidth[i])
    nrg = zerofloat(ntg)
    ref.apply1(g,nrg)

    #Make RMS of image equal to 1
    ########################################################
    ww = WaveletWarpingCBGN()
    f = ww.makerms1(itminf,itmaxf,f)
    g = ww.makerms1(itming,itmaxg,g)
    nrg = ww.makerms1(itming,itmaxg,nrg)
    
    #plotting
    pngDir = None
    title= "g halfwidth = "+str(halfwidth[i])
    ntg = len(g)
    st = Sampling(ntg,dt,0.0)
    vmin,vmax = 0*dt,ntg*dt
    vlabel,vminmax,vint = "Time (s)",[tming,tmaxg],None
    hint1 = 2.5
    hmin1 = -5.0
    hmax1 = 5.0
    hlabel = ["g","Noise Reduced g"]
    hminmax = [[hmin1,hmax1],[hmin1,hmax1]]
    hint = [None,None]
    hsize,vsize = 960,560
    tilespacing = 5
    plotting.plotTracesSideBySide(st,[g,nrg],\
    vlabel=vlabel,vminmax=vminmax,vint=vint,\
    hlabel=hlabel,hminmax=hminmax,hint=hint,\
    tilespacing = tilespacing,\
    title=title,pngDir=pngDir,\
    paper=True,onecol=True)
    
    """
    title= "f halfwidth = "+str(halfwidth[i])
    ntf = len(f)
    st = Sampling(ntf,dt,0.0)
    vmin,vmax = 0*dt,ntf*dt
    vlabel,vminmax,vint = "Time (s)",[tminf,tmaxf],None
    hint1 = 2.5
    hmin1 = -5.0
    hmax1 = 5.0
    hlabel = ["f"]
    hminmax = [[hmin1,hmax1]]
    hint = [None]
    hsize,vsize = 960,560
    tilespacing = 5
    plotting.plotTracesSideBySide(st,[f],\
    vlabel=vlabel,vminmax=vminmax,vint=vint,\
    hlabel=hlabel,hminmax=hminmax,hint=hint,\
    tilespacing = tilespacing,\
    title=title,pngDir=pngDir,\
    paper=True,onecol=True)
    """







 
#81 coefficients are supplied to the b and c. In this test, we will shift coefficients 
#from c to b. The initial starting point will be b having 1 coefficient and c having 81 
#coefficients. We will only look at filters that have an odd number of wavelet coefficients.
#The reasoning for this is that we will not be forcing a zero at zero frequency. 
#We will plot the final rms for each pair of coefficients. The only parameter that will
#change is the number of coefficients in c and b. This test will only be for 1D data.
def goTestNcNbTradeOff():
  #Synthetic parameters
  nt,ni,randomi = 1081,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf
  sfac = 0.0001

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  du = computeBackDiff(u)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)
  tmin = tmin
  tmax = tmax
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)

  nbfinal = 81
  nb = 1
  kb = 0
  nc = 81
  kc = -40
  n = nbfinal/2+1
  firstriri = zerofloat(n)
  lastrfrf = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    #Estimate wavelet
    niter = 500
    warp = Warper()
    ww = WaveletWarpingCBGN()
    ww.setTimeRange(tmin,tmax)
    ww.setMaxPercentChange(0.001)
    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    dump(bw)

    #Get iteration information
    riri = ww.getRiRi()
    rfrf = ww.getRfRf()
    vv = ww.getVV()
    stepl = ww.getStepLength()
    deltamag = ww.getDeltaMag()
    condnum = ww.getCondNum()
    lastiter = ww.getLastIter()
    firstriri[i] = riri[0]
    lastrfrf[i] = rfrf[lastiter]
    nb = nb+2
    kb = -nb/2+1
    nc = nc - 2
    kc = -nc/2+1

  pngDir = "./tests/"
  maxriri = max(firstriri)
  title = "first rms ri (black) last rms rf (blue)"+str(nrmsf)
  si = Sampling(n,2.0,1.0)
  plotSequenceMeas(si,[firstriri,lastrfrf],amax=maxriri,hlabel="nb=81-nc",labels="Sum of square diff",pngDir=pngDir,title=title)

  maxf = max(f)
  maxfd2 = maxf/2.0
  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,g],amax=amax,tmark=tmark,\
  labels=["f","g"],pngDir=pngDir,\
  title=title)

#81 coefficients are supplied to the b and c. In this test, we will shift coefficients 
#from c to b. The initial starting point will be b having 1 coefficient and c having 81 
#coefficients. We will only look at filters that have an odd number of wavelet coefficients.
#The reasoning for this is that we will not be forcing a zero at zero frequency. 
#We will plot the final rms for each pair of coefficients. The only parameter that will
#change is the number of coefficients in c and b. This test will only be for 1D data.
def goTestNcNbLessSqueezeNoiseTradeOff():
  #Synthetic parameters
  nt,ni,randomi = 1081,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  r0,r1 = 1.6,1.4#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf
  sfac = 0.0001

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  du = computeBackDiff(u)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)
  tmin = tmin
  tmax = tmax
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)

  nbfinal = 81
  nb = 1
  kb = 0
  nc = 81
  kc = -40
  n = nbfinal/2+1
  firstriri = zerofloat(n)
  lastrfrf = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    #Estimate wavelet
    niter = 500
    warp = Warper()
    ww = WaveletWarpingCBGN()
    ww.setTimeRange(tmin,tmax)
    ww.setMaxPercentChange(0.001)
    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    dump(bw)

    #Get iteration information
    riri = ww.getRiRi()
    rfrf = ww.getRfRf()
    vv = ww.getVV()
    stepl = ww.getStepLength()
    deltamag = ww.getDeltaMag()
    condnum = ww.getCondNum()
    lastiter = ww.getLastIter()
    firstriri[i] = riri[0]
    lastrfrf[i] = rfrf[lastiter]
    nb = nb+2
    kb = -nb/2+1
    nc = nc - 2
    kc = -nc/2+1

  pngDir = "./tests/"
  maxriri = max(firstriri)
  title = "LessSqueeze first rms ri (black) last rms rf (blue)"+str(nrmsf)
  si = Sampling(n,2.0,1.0)
  plotSequenceMeas(si,[firstriri,lastrfrf],amax=maxriri,hlabel="nb=81-nc",labels="Sum of square diff",pngDir=pngDir,title=title)

  maxf = max(f)
  maxfd2 = maxf/2.0
  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,g],amax=amax,tmark=tmark,\
  labels=["f","g"],pngDir=pngDir,\
  title=title)

def goTestTradeOffMinLessSqueezeNoise():
  #####1D##########
  #Synthetic parameters
  nt,ni,randomi = 1081,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  nb,kb = 21,-10#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 61,-30# sampling for wavelet H 
  niter = 500
  r0,r1 = 1.6,1.4#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf
  sfac = 0.0001

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  du = computeBackDiff(u)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)
  tmin = tmin
  tmax = tmax
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)

  #Estimate wavelet
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
  ww.setMaxPercentChange(0.001)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  dump(bw)

  #Get iteration information
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  vv = ww.getVV()
  stepl = ww.getStepLength()
  deltamag = ww.getDeltaMag()
  condnum = ww.getCondNum()

  #Get known wavelet
  dk = getWavelet(freqd,decayd,nc,kc,mpd)
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Normalize wavelets
  ncw = normalizeMAAWOS(cw)
  nck = normalizeMAAWOS(ck)
  nbw = normalizeMAAWOS(bw)
  nbk = normalizeMAAWOS(bk)
  ndw = normalizeMAAWOS(dw)
  ndk = normalizeMAAWOS(dk)
  #ncw = cw
  #nck = ck
  #nbw = bw
  #nbk = bk
  #ndk = dk

  #Create simply warped g
  sg = warp.applyS(u,g)
  sq = warp.applyS(u,q)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletC(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyH(nc,kc,h,sg)
  nhw = normalizeMAAWOS(h)

  #deconvolve f
  bf = ww.applyH(nb,kb,bw,f)
  
  #Create warped g
  bg = ww.applyH(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,cw,sbg)

  #plotting
  #pngDir = None
  pngDir = "./tests/"

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title1 = "LessSqueezeNoiseWaveletc"
  title2 = "LessSqueezeNoiseWaveletd"
  title3 = "LessSqueezeNoiseWaveletb"
  title4 = "LessSqueezeNoiseWaveleth"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[nhw],title=title4,pngDir=pngDir,onecol=True)

  if len(rfrf)>0:
    maxriri = max(riri)
    lastiter = ww.getLastIter()
    minrfrf = rfrf[lastiter]
    title = "LessSqueezeNoiseRMS ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="RMS of residuals",pngDir=pngDir,title=title)

    maxstepl = max(stepl)
    title = "LessSqueezeNoiseStepLength "
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[stepl],amax=maxstepl,hlabel="Iterations",labels="Step Length",pngDir=pngDir,title=title)

    maxvv= max(vv)
    title = "Noise 2 norm of neg. gradient square"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[vv],amax=maxvv,tmark=maxvv/100.0,hlabel="Iterations",labels="V'V",pngDir=None,title=title)


    maxdm = max(deltamag)
    title = "Noise Delta Mag"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[deltamag],amax=maxdm,tmark=maxdm/100.0,hlabel="Iterations",labels="Delta Mag",pngDir=None,title=title)

    maxcn = max(condnum)
    title = "Noise X Condition Number"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[condnum],amax=maxcn,tmark=maxcn/100.0,hlabel="Iterations",labels="X Condition Number",pngDir=None,title=title)
  #ssdzoom,fiter = zoomConverge(ssd,1.0)
  #nsubiter = len(ssdzoom)
  #ssub = Sampling(nsubiter,1.0,fiter)
  #med = getMedian(ssdzoom)
  #ssdzoom = sub(ssdzoom,med)
  #sp = SimplePlot()
  #sp.addPoints(ssub,ssdzoom)
  #sp.addTitle("Sum of square differences between f and CSBg")

  maxf = max(f)
  maxfd2 = maxf/2.0
  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title = "LessSqueezeNoisefg nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,csbg,hsg,g,sg],amax=amax,tmark=tmark,\
  labels=["f","CSBg","HSg","g","Sg"],pngDir=pngDir,\
  title=title)
  #plotSequences(st,[bf,sbg],amax=[1,1],tmark=[0.1,0.1],\
  #labels=["Bf","SBg"],pngDir=pngDir,\
  #title=title)



def goTestTradeOffMinNoNoise():
  #####1D##########
  #Synthetic parameters
  nt,ni,randomi = 1081,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 77,-38# sampling for wavelet H 
  niter = 500
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 0.0001

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  du = computeBackDiff(u)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)
  tmin = tmin
  tmax = tmax
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)

  #Estimate wavelet
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
  ww.setMaxPercentChange(0.001)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  dump(bw)

  #Get iteration information
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  vv = ww.getVV()
  stepl = ww.getStepLength()
  deltamag = ww.getDeltaMag()
  condnum = ww.getCondNum()

  #Get known wavelet
  dk = getWavelet(freqd,decayd,nc,kc,mpd)
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Normalize wavelets
  ncw = normalizeMAAWOS(cw)
  nck = normalizeMAAWOS(ck)
  nbw = normalizeMAAWOS(bw)
  nbk = normalizeMAAWOS(bk)
  ndw = normalizeMAAWOS(dw)
  ndk = normalizeMAAWOS(dk)
  #ncw = cw
  #nck = ck
  #nbw = bw
  #nbk = bk
  #ndk = dk

  #Create simply warped g
  sg = warp.applyS(u,g)
  sq = warp.applyS(u,q)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletC(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)
  nhw = normalizeMAAWOS(h)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyH(nc,kc,h,sg)

  #deconvolve f
  bf = ww.applyH(nb,kb,bw,f)
  
  #Create warped g
  bg = ww.applyH(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,cw,sbg)

  #plotting
  #pngDir = None
  pngDir = "./tests/"

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title1 = "NoNoiseWaveletc"
  title2 = "NoNoiseWaveletd"
  title3 = "NoNoiseWaveletb"
  title4 = "NoNoiseWaveleth"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[nhw],title=title4,pngDir=pngDir,onecol=True)

  if len(rfrf)>0:
    maxriri = max(riri)
    lastiter = ww.getLastIter()
    minrfrf = rfrf[lastiter]
    title = "NoNoiseRMS rms ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,amin=minrfrf,hlabel="Iterations",labels="RMS of residuals",pngDir=pngDir,title=title)

    maxstepl = max(stepl)
    title = "NoNoiseStepLength"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[stepl],amax=maxstepl,hlabel="Iterations",labels="Step Length",pngDir=pngDir,title=title)

    maxvv= max(vv)
    title = "2 norm of neg. gradient square"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[vv],amax=maxvv,tmark=maxvv/100.0,hlabel="Iterations",labels="V'V",pngDir=None,title=title)


    maxdm = max(deltamag)
    title = "Delta Mag"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[deltamag],amax=maxdm,tmark=maxdm/100.0,hlabel="Iterations",labels="Delta Mag",pngDir=None,title=title)

    maxcn = max(condnum)
    title = "X Condition Number"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[condnum],amax=maxcn,tmark=maxcn/100.0,hlabel="Iterations",labels="X Condition Number",pngDir=None,title=title)
  #ssdzoom,fiter = zoomConverge(ssd,1.0)
  #nsubiter = len(ssdzoom)
  #ssub = Sampling(nsubiter,1.0,fiter)
  #med = getMedian(ssdzoom)
  #ssdzoom = sub(ssdzoom,med)
  #sp = SimplePlot()
  #sp.addPoints(ssub,ssdzoom)
  #sp.addTitle("Sum of square differences between f and CSBg")

  maxf = max(f)
  maxfd2 = maxf/2.0
  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title = "NoNoisefg nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,csbg,hsg,g,sg],amax=amax,tmark=tmark,\
  labels=["f","CSBg","HSg","g","Sg"],pngDir=pngDir,\
  title=title)
  #plotSequences(st,[bf,sbg],amax=[1,1],tmark=[0.1,0.1],\
  #labels=["Bf","SBg"],pngDir=pngDir,\
  #title=title)

def goTestTradeOffMinNoise():
  #####1D##########
  #Synthetic parameters
  nt,ni,randomi = 1081,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  nb,kb = 45,-22#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 37,-18# sampling for wavelet H 
  niter = 500
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf
  sfac = 0.0001

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  du = computeBackDiff(u)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)
  tmin = tmin
  tmax = tmax
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)

  #Estimate wavelet
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
  ww.setMaxPercentChange(0.001)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  dump(bw)

  #Get iteration information
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  vv = ww.getVV()
  stepl = ww.getStepLength()
  deltamag = ww.getDeltaMag()
  condnum = ww.getCondNum()

  #Get known wavelet
  dk = getWavelet(freqd,decayd,nc,kc,mpd)
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Normalize wavelets
  ncw = normalizeMAAWOS(cw)
  nck = normalizeMAAWOS(ck)
  nbw = normalizeMAAWOS(bw)
  nbk = normalizeMAAWOS(bk)
  ndw = normalizeMAAWOS(dw)
  ndk = normalizeMAAWOS(dk)
  #ncw = cw
  #nck = ck
  #nbw = bw
  #nbk = bk
  #ndk = dk

  #Create simply warped g
  sg = warp.applyS(u,g)
  sq = warp.applyS(u,q)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletC(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyH(nc,kc,h,sg)
  nhw = normalizeMAAWOS(h)

  #deconvolve f
  bf = ww.applyH(nb,kb,bw,f)
  
  #Create warped g
  bg = ww.applyH(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,cw,sbg)

  #plotting
  #pngDir = None
  pngDir = "./tests/"

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title1 = "NoiseWaveletc"
  title2 = "NoiseWaveletd"
  title3 = "NoiseWaveletb"
  title4 = "NoiseWaveleth"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[nhw],title=title4,pngDir=pngDir,onecol=True)

  if len(rfrf)>0:
    maxriri = max(riri)
    lastiter = ww.getLastIter()
    minrfrf = rfrf[lastiter]
    title = "NoiseRMS ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="RMS of residuals",pngDir=pngDir,title=title)

    maxstepl = max(stepl)
    title = "NoiseStepLength "
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[stepl],amax=maxstepl,hlabel="Iterations",labels="Step Length",pngDir=pngDir,title=title)

    maxvv= max(vv)
    title = "Noise 2 norm of neg. gradient square"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[vv],amax=maxvv,tmark=maxvv/100.0,hlabel="Iterations",labels="V'V",pngDir=None,title=title)


    maxdm = max(deltamag)
    title = "Noise Delta Mag"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[deltamag],amax=maxdm,tmark=maxdm/100.0,hlabel="Iterations",labels="Delta Mag",pngDir=None,title=title)

    maxcn = max(condnum)
    title = "Noise X Condition Number"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[condnum],amax=maxcn,tmark=maxcn/100.0,hlabel="Iterations",labels="X Condition Number",pngDir=None,title=title)
  #ssdzoom,fiter = zoomConverge(ssd,1.0)
  #nsubiter = len(ssdzoom)
  #ssub = Sampling(nsubiter,1.0,fiter)
  #med = getMedian(ssdzoom)
  #ssdzoom = sub(ssdzoom,med)
  #sp = SimplePlot()
  #sp.addPoints(ssub,ssdzoom)
  #sp.addTitle("Sum of square differences between f and CSBg")

  maxf = max(f)
  maxfd2 = maxf/2.0
  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title = "Noisefg nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,csbg,hsg,g,sg],amax=amax,tmark=tmark,\
  labels=["f","CSBg","HSg","g","Sg"],pngDir=pngDir,\
  title=title)
  #plotSequences(st,[bf,sbg],amax=[1,1],tmark=[0.1,0.1],\
  #labels=["Bf","SBg"],pngDir=pngDir,\
  #title=title)


 
#81 coefficients are supplied to the b and c. In this test, we will shift coefficients 
#from c to b. The initial starting point will be b having 1 coefficient and c having 81 
#coefficients. We will only look at filters that have an odd number of wavelet coefficients.
#The reasoning for this is that we will not be forcing a zero at zero frequency. 
#We will plot the final rms for each pair of coefficients. The only parameter that will
#change is the number of coefficients in c and b. This test will only be for 1D data.
def goTestNcNb2DTradeOff():
  #get sino images
  x0,nx = 250,100
  f,g,u = getSinoImage(x0,nx)
  nx = len(f)
  ng = len(g[0])
  nu = len(u[0])

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  nbfinal = 41
  nb = 1
  kb = 0
  nc = 41
  kc = -20
  n = nbfinal/2+1
  firstriri = zerofloat(n)
  lastrfrf = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    #Estimate wavelet
    sfac = 0.0001
    niter = 500
    warp = Warper()
    ww = WaveletWarpingCBGN()
    ww.setTimeRange(tmin,tmax)
    ww.setMaxPercentChange(0.001)
    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    dump(bw)

    #Get iteration information
    riri = ww.getRiRi()
    rfrf = ww.getRfRf()
    vv = ww.getVV()
    stepl = ww.getStepLength()
    deltamag = ww.getDeltaMag()
    condnum = ww.getCondNum()
    lastiter = ww.getLastIter()
    firstriri[i] = riri[0]
    lastrfrf[i] = rfrf[lastiter]
    nb = nb+2
    kb = -nb/2+1
    nc = nc - 2
    kc = -nc/2+1

  pngDir = "./tests/"
  maxriri = max(firstriri)
  minrfrf = min(lastrfrf)
  title = "Sinopec first rms ri (black) last rms rf (blue)"
  si = Sampling(n,2.0,1.0)
  plotSequenceMeas(si,[firstriri,lastrfrf],amax=maxriri,amin=minrfrf,hlabel="nb=41-nc",labels="Sum of square diff",pngDir=pngDir,title=title)

#81 coefficients are supplied to the b and c. In this test, we will shift coefficients 
#from c to b. The initial starting point will be b having 1 coefficient and c having 81 
#coefficients. We will only look at filters that have an odd number of wavelet coefficients.
#The reasoning for this is that we will not be forcing a zero at zero frequency. 
#We will plot the final rms for each pair of coefficients. The only parameter that will
#change is the number of coefficients in c and b. This test will only be for 1D data.
def goTestNcNb2DTradeOffMin():
  #get sino images
  x0,nx = 250,100
  f,g,u = getSinoImage(x0,nx)
  nx = len(f)
  ng = len(g[0])
  nu = len(u[0])

  #Wavelet estimation parameters
  nb,kb = 7,-3#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 35,-17#sampling for wavelet H 
  nd,kd = 35,-17# sampling for wavelet H 

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  niter = 100
  sfac = 0.0001
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setParallel(False)
  ww.setMaxPercentChange(0.001)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setLineSearchMinScale(0.0)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  one = zerofloat(nb)
  one[-kb] = 1
  hw = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  vv = ww.getVV()
  stepl = ww.getStepLength()
  deltamag = ww.getDeltaMag()
  condnum = ww.getCondNum()

  #Normalize wavelets
  nhw = normalizeMAAWOS(hw)
  ncw = normalizeMAAWOS(cw)
  nbw = normalizeMAAWOS(bw)
  ndw = normalizeMAAWOS(dw)

  #Create warped g
  bg = ww.applyH(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,cw,sbg)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)
  cb = ww.applyH(nc,kc,cw,bw)

  #Create shaped Sg
  sg = warp.applyS(u,g)
  hsg = ww.applyH(nc,kc,hw,sg)

  #Simply warp g
  sg = warp.applyS(u,g)

  #rms scale
  rmsf = ww.rms(f)
  rmssg = ww.rms(sg)
  rmshsg = ww.rms(hsg)
  rmscsbg = ww.rms(csbg)
  f = div(f,rmsf)
  sg = div(sg,rmssg)
  hsg = div(hsg,rmshsg)
  csbg = div(csbg,rmscsbg)

  #plotting
  pngDir = None
  pngDir = "./tests/"

  
  dt = 0.004
  ft = 0.000
  nt = len(f[0])
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  plotUnknownWavelets(Sampling(nc,dt,kc*dt),[ncw],title="Sinoc",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(nb,dt,kb*dt),[nbw],title="Sinob",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(nd,dt,kd*dt),[ndw],title="Sinod",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(nc,dt,kc*dt),[nhw],title="Sinoh",pngDir=pngDir,onecol=True)
  du = computeBackDiff2D(u)

  maxf = max(f)
  maxfd2 = maxf/2.0
  ntg = len(g[0])

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  #suu = Sampling(nuu,dtgu*dt,0.0)
  #sgu = Sampling(ngu,dtgu*dt,0.0)
  #ssug = Sampling(nuu,dtgu*dt,0.0)
  su = Sampling(ntg,dt,ft)
  sx = Sampling(nx,0.015,0.0)
  ut = mul(u,dt)
  print "nf2 = "+str(len(f))
  print "nf1 = "+str(len(f[0]))
  print "ng2 = "+str(len(g))
  print "ng1 = "+str(len(g[0]))
  #plotImage(suu,sx,uu,title="uu")
  #plotImage(sgu,sx,ug,title="ug")
  #plotImage(ssug,sx,sug,title="sug")
  frms = ww.rms(f)
  frmsI = div(f,frms)
  sgrms = ww.rms(sg)
  sgrmsI = div(sg,sgrms)
  hsgrms = ww.rms(hsg)
  hsgrmsI = div(hsg,hsgrms)
  csbgrms = ww.rms(csbg)
  csbgrmsI = div(csbg,csbgrms)
  clip = 2
  title = "Sinofg"
  plot4ImagesSideBySide(st,sx, frmsI, csbgrmsI, hsgrmsI, sgrmsI, tmin*dt, tmax*dt, clip, 
  "PP time", title=title, pngDir=pngDir,paper=True, twocol=True)

  #plotImage(st,sx,sg,title="Sg")
  #plotImage(st,sx,du,title="du")
  #plotImage(st,sx,dlsug,title="DLSUg")
  #plotImage(st,sx,dlsug1,title="DLSUg1")

  #amax = [maxf,maxf,maxf,maxf]
  #tmark = [maxfd2,maxfd2,maxfd2,maxfd2]
  #plotSequences(st,[f[0],dlsug[0],sg[0]],amax=amax,tmark=tmark,\
  #labels=["f","HWg","Sg"])

  if len(rfrf)>0:
    maxriri = max(riri)
    lastiter = ww.getLastIter()
    minrfrf = rfrf[lastiter]
    title = "SinoRMS rms ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,amin=minrfrf,hlabel="Iterations",labels="RMS of residuals",pngDir=pngDir,title=title)

    maxstepl = max(stepl)
    minstepl = min(stepl)
    title = "SinoStepLength"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[stepl],amax=maxstepl,amin=minstepl,hlabel="Iterations",labels="StepLength",pngDir=pngDir,title=title)

    maxvv= max(vv)
    title = "2 norm of neg. gradient square"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[vv],amax=maxvv,tmark=maxvv/100.0,hlabel="Iterations",labels="V'V",pngDir=None,title=title)


    maxdm = max(deltamag)
    title = "Delta Mag"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[deltamag],amax=maxdm,tmark=maxdm/100.0,hlabel="Iterations",labels="Delta Mag",pngDir=None,title=title)

    maxcn = max(condnum)
    title = "X Condition Number"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[condnum],amax=maxcn,tmark=maxcn/100.0,hlabel="Iterations",labels="X Condition Number",pngDir=None,title=title)

def goTestSinopecNoiseReduction():
  #get sino images
  x0,nx = 250,100
  f,g,u = getSinoImage(x0,nx)
  SimplePlot.asPixels(f)
  SimplePlot.asPixels(g)
  ntg = len(g[0])
  ntf = len(f[0])
  tmin,tmax=100,500
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  halfwidth = 0
  ref = RecursiveExponentialFilter(halfwidth)
  nrg = zerofloat(ntg,nx)
  ref.apply1(g,nrg)
  #Plotting#
  dt = 0.004
  stg = Sampling(ntg,dt,0.0)
  stf = Sampling(ntf,dt,0.0)
  sx = Sampling(nx,0.0150,x0*0.0150)


  #Make RMS of image equal to 1
  ########################################################
  f = ww.makerms1(f)
  g = ww.makerms1(225,800,g)
  nrg = ww.makerms1(225,800,nrg)

  #Plot g with nrg (noise reduced g)
  #########################################################
  pngDir = "./tests/"
  #pngDir = None
  title= "Noise Reduction for g"
  vmin,vmax = 225*dt,800*dt
  hmin,hmax = None,None
  vlabel,vminmax,vint = "PS time (s)",None,1.0
  hlabel,hminmax,hint = ["g","Noise reduced g"],None,0.5
  clip = 2
  plotting.plotImagesSideBySide(stg,sx,[g,nrg],\
  vlabel=vlabel,vminmax=[vmin,vmax],vint=vint,\
  hlabel=hlabel,hminmax=None,hint=hint,\
  clip=clip,
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)
  #########################################################
  
  #Plot g with nrg (noise reduced g)
  #########################################################
  pngDir = "./tests/"
  #pngDir = None
  title= "Justg"
  vmin,vmax = 225*dt,800*dt
  hmin,hmax = None,None
  vlabel,vminmax,vint = "PS time (s)",None,1.0
  hlabel,hminmax,hint = ["g"],None,0.5
  clip = 2
  plotting.plotImagesSideBySide(stg,sx,[g],\
  vlabel=vlabel,vminmax=[vmin,vmax],vint=vint,\
  hlabel=hlabel,hminmax=None,hint=hint,\
  clip=clip,
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)
  #########################################################

  #Plot low noise f for comparison with g and nrg
  #########################################################
  pngDir = "./tests/"
  title= "f"
  vmin,vmax = tmin*dt,tmax*dt
  hmin,hmax = None,None
  vlabel,vminmax,vint = "PP time (s)",None,1.0
  hlabel,hminmax,hint = ["f"],None,0.5
  clip = 2
  plotting.plotImagesSideBySide(stf,sx,[f],\
  vlabel=vlabel,vminmax=[vmin,vmax],vint=vint,\
  hlabel=hlabel,hminmax=None,hint=hint,\
  clip=clip,
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)
  #########################################################

  itrace = 50
  #Plot trace 300 in g and nrg
  pngDir = None
  title= "g and ngr traces"
  vmin,vmax = 0*dt,1081*dt
  hmin,hmax = -6.0,6.0
  vlabel,vminmax,vint = "PS time (s)",None,1.0
  hlabel,hminmax,hint = ["g","nrg"],[[hmin,hmax],[hmin,hmax]],2.0
  hsize,vsize = 960,560
  plotting.plotTracesSideBySide(stg,[g[itrace],nrg[itrace]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

def goTestNcNb2DNoiseReductionTradeOff():
  #get sino images
  x0,nx = 250,100
  f,g,u = getSinoImage(x0,nx)
  halfwidth=3
  ref = RecursiveExponentialFilter(halfwidth)
  ref.apply1(g,g)
  nx = len(f)
  ng = len(g[0])
  nu = len(u[0])

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  nbfinal = 41
  nb = 1
  kb = 0
  nc = 41
  kc = -21
  n = nbfinal/2+1
  firstriri = zerofloat(n)
  lastrfrf = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    #Estimate wavelet
    sfac = 0.0001
    niter = 500
    warp = Warper()
    ww = WaveletWarpingCBGN()
    ww.setTimeRange(tmin,tmax)
    ww.setMaxPercentChange(0.001)
    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    dump(bw)

    #Get iteration information
    riri = ww.getRiRi()
    rfrf = ww.getRfRf()
    vv = ww.getVV()
    stepl = ww.getStepLength()
    deltamag = ww.getDeltaMag()
    condnum = ww.getCondNum()
    lastiter = ww.getLastIter()
    firstriri[i] = riri[0]
    lastrfrf[i] = rfrf[lastiter]
    nb = nb+2
    kb = -nb/2+1
    nc = nc - 2
    kc = -nc/2+1

  #Plotting
  #plot rms of residuals with shaping filter (firstriri) 
  #and rms of residuals with c and b (lastrfrf)
  #############################################################
  pngDir = "./tests/"
  title = "NR"+str(halfwidth)+str(nbfinal)+"Sinopec Noise Reduction first rms ri (black) last rms rf (blue)"
  maxriri = max(firstriri)
  minrfrf = min(lastrfrf)
  print "minrfrf = "+str(minrfrf)
  print "maxriri = "+str(maxriri)
  si = Sampling(n,2.0,1.0)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrfrf,maxriri],None
  hlabel,hminmax,hint = "nb=42-nc",[1.0,nbfinal],5.0
  plotting.plotMeasInSamePlot(si, [firstriri,lastrfrf],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  #############################################################

def goTestNbNhIncrease():
  #####1D##########
  #Synthetic parameters
  nt,ni,randomi = 1081,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  niter = 500
  sfac = 0.0001
  maxpc = 0.001

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  du = computeBackDiff(u)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)
  tmin = tmin
  tmax = tmax
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)

  nhfinal = 41
  nhinitial = 21
  khinitial = -10
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcb = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    warp = Warper()
    ww = WaveletWarpingCBGN()
    ww.setTimeRange(tmin,tmax)
    ww.setMaxPercentChange(maxpc)
    #Shaping filter
    b = zerofloat(nb)
    b[-kb] = 1.0
    h = ww.getWaveletC(nh,kh,nb,kb,b,u,f,g)
    rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,b,u,f,g))


    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    dump(bw)

    #Get iteration information
    rfrf = ww.getRfRf()
    vv = ww.getVV()
    stepl = ww.getStepLength()
    condnum = ww.getCondNum()
    lastiter = ww.getLastIter()
    rmsrcb[i] = rfrf[lastiter]
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb = "+str(rmsrcb)
    nb = nb+1
    kb = -nb/2+1
    nh = nh+1
    kh = -nh/2+1

  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  #pngDir = "./tests/"
  pngDir = None
  title = "Increaseby1"+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max(rmsrcb),max(rmsrh)])
  minrms = min([min(rmsrcb),min(rmsrh)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(rmsrh))
  print "rmsrcb = "+str(len(rmsrcb))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+20",[nhinitial,nhfinal],10.0
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcb],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  #############################################################

  #############1 Plot######################################
  pngDir = "./tests/"
  #pngDir = None
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,1081*dt
  hmin,hmax = -8.0,8.0
  vlabel,vminmax,vint = "time (s)",[vmin,vmax],None
  hlabel,hminmax,hint = ["f","g"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax]],2.0
  hsize,vsize = 960,560
  title= "Increaseby1 noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase fg nb and nh nc constant first rms ri (black) last rms rf (blue)"  
  plotting.plotTracesSideBySide(st,[f,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)
  ##########################################################

def goTest1SinoNbNhIncrease():
  #get sino images
  x0 = 250
  f,g,u = getSinoTrace(x0)
  ng = len(g)
  nu = len(u)

  #Wavelet estimation parameters
  nb,kb = 2,0#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 21,-10# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 

  #set tmin and tmax 
  tmin = 100
  tmax = 500
  sfac = 0.0001
  maxpc = 0.000
  minsl = 0.000
  niter = 200

  nhfinal = 41
  nhinitial = 21
  khinitial = -10
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcb = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    warp = Warper()
    ww = WaveletWarpingCBGN()

    ww.setMaxPercentChange(maxpc)
    ww.setTimeRange(tmin,tmax)
    ww.setLineSearchMinScale(minsl)
    #Shaping filter
    b = zerofloat(nb)
    b[-kb] = 1.0
    h = ww.getWaveletC(nh,kh,nb,kb,b,u,f,g)
    rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,b,u,f,g))

    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)
    dump(bw)

    #Get iteration information
    rfrf = ww.getRfRf()
    vv = ww.getVV()
    stepl = ww.getStepLength()
    condnum = ww.getCondNum()
    lastiter = ww.getLastIter()
    rmsrcb[i] = rfrf[lastiter]
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb = "+str(rmsrcb)
    nb = nb+1
    kb = -nb/2+1
    nh = nh+1
    kh = -nh/2+1

  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  #pngDir = "./tests/"
  pngDir = None
  title = "Increaseby1"+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" nhfinal"+str(nhfinal)+"Sino Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max(rmsrcb),max(rmsrh)])
  minrms = min([min(rmsrcb),min(rmsrh)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(rmsrh))
  print "rmsrcb = "+str(len(rmsrcb))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+20",[nhinitial,nhfinal],10.0
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcb],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  #############################################################


def goTestNbNh2DIncrease():
  #get sino images
  x0,nx = 250,1
  f,g,u = getSinoImage(x0,nx)
  halfwidth=0
  if halfwidth>0:
    ref = RecursiveExponentialFilter(halfwidth)
    ref.apply1(g,g)
  nx = len(f)
  ng = len(g[0])
  nu = len(u[0])

  #set tmin and tmax 
  tmin = 100
  tmax = 500
  SimplePlot.asPixels(g)
  SimplePlot.asPixels(f)

  nhfinal = 41
  nb = 1
  kb = 0
  nc = 21
  kc = -10
  nh = 21
  kh = -10
  n = (nhfinal-nh)+1
  rmsrh = zerofloat(n)
  rmsrcb = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    sfac = 0.0001
    niter = 500
    pc = 0.000001
    warp = Warper()
    ww = WaveletWarpingCBGN()
    ww.setTimeRange(tmin,tmax)
    ww.setMaxPercentChange(pc)
    #Shaping filter
    b = zerofloat(nb)
    b[-kb] = 1.0
    h = ww.getWaveletC(nh,kh,nb,kb,b,u,f,g)
    rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,b,u,f,g))


    """
    #First guesses of c and b. 
    if (i==0):
      b = zerofloat(nb)
      b[-kb] = 1.0
      c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
    else:
      c = cw
      b = usePreviousbw(nb,bw)
      dump(bw)
      dump(b)
      dump(cw)
      dump(c)
      title = "Increase start new"
    """
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
    title = "Increase start regular"

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get iteration information
    rfrf = ww.getRfRf()
    vv = ww.getVV()
    stepl = ww.getStepLength()
    condnum = ww.getCondNum()
    lastiter = ww.getLastIter()
    rmsrcb[i] = rfrf[lastiter]
    nb = nb+1
    kb = -nb/2+1
    nh = nh+1
    kh = -nh/2+1
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb = "+str(rmsrcb)


  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  #pngDir = "./tests/"
  pngDir = None
  title = title+"NewStopBy1"+"halfwidth = "+str(halfwidth)+"niter "+str(niter)+" pc"+str(pc)+" sfac"+str(sfac)+" nhfinal"+str(nhfinal)+"Sinopec Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"+" sfac"+str(sfac)+" pc"+str(pc)
  maxrms = max([max(rmsrcb),max(rmsrh)])
  minrms = min([min(rmsrcb),min(rmsrh)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(rmsrh))
  print "rmsrcb = "+str(len(rmsrcb))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+20",[nc,nhfinal],5.0
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcb],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  #############################################################

def goTestNbNhAll2DIncrease():
  #get sino images
  x0,nx = 0,721
  f,g,u = getSinoImage(x0,nx)
  halfwidth=0
  if halfwidth>0:
    ref = RecursiveExponentialFilter(halfwidth)
    ref.apply1(g,g)
  nx = len(f)
  ng = len(g[0])
  nu = len(u[0])

  #set tmin and tmax 
  tmin = 100
  tmax = 500
  SimplePlot.asPixels(g)
  SimplePlot.asPixels(f)

  nhfinal = 41
  nb = 1
  kb = 0
  nc = 21
  kc = -10
  nh = 21
  kh = -10
  n = (nhfinal-nh)/2+1
  rmsrh = zerofloat(n)
  rmsrcb = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    #Estimate wavelet
    sfac = 0.0001
    niter = 500
    warp = Warper()
    ww = WaveletWarpingCBGN()
    ww.setTimeRange(tmin,tmax)
    ww.setMaxPercentChange(0.001)
    #Shaping filter
    b = zerofloat(nb)
    b[-kb] = 1.0
    h = ww.getWaveletC(nh,kh,nb,kb,b,u,f,g)
    rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,b,u,f,g))


    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletC(nc,kc,nb,kb,b,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get iteration information
    rfrf = ww.getRfRf()
    vv = ww.getVV()
    stepl = ww.getStepLength()
    deltamag = ww.getDeltaMag()
    condnum = ww.getCondNum()
    lastiter = ww.getLastIter()
    rmsrcb[i] = rfrf[lastiter]
    nb = nb+2
    kb = -nb/2+1
    nh = nh+2
    kh = -nh/2+1
    print "rmsrh = "+str(rmsrh)
    print "rmsrcb = "+str(rmsrcb)


  #Plotting
  #plot rms of residuals with shaping filter (rmsrh) 
  #and rms of residuals with c and b (rmsrcb)
  #############################################################
  pngDir = "./tests/"
  #pngDir = None
  title = "Increase"+"halfwidth = "+str(halfwidth)+" nhfinal"+str(nhfinal)+"AllSinopec Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = max([max(rmsrcb),max(rmsrh)])
  minrms = min([min(rmsrcb),min(rmsrh)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(rmsrh))
  print "rmsrcb = "+str(len(rmsrcb))
  si = Sampling(n,2.0,21)
  color=[Color.BLACK,Color.BLUE]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+20",[21.0,nhfinal],5.0
  plotting.plotMeasInSamePlot(si, [rmsrh,rmsrcb],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True, onecol=True, twocol=None)
  #############################################################

    
def goSino1000ShapingFilter():
  #get sino images
  x0,nx = 250,100
  #f,g,u = getSinoTrace(x0)#getSinoImage(x0,nx)
  #nt = len(f)
  f,g,u = getSinoImage(x0,nx)
  nt = len(f[0])

  #Warp g to f using new method
  warp = Warper()
  sg = warp.applyS(u,g)

  #Estimate Wavelet
  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 1001,-500#sampling for wavelet H 
  nd,kd = 81,-40# sampling for wavelet H 
  niter = 300
  sfac = 0.1
  tmin = 100
  tmax = 500
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setLineSearchMinScale(0.001);
  ww.setTimeRange(tmin,tmax)
  #Create shaped squeezed g (HSg)
  one = zerofloat(nb)
  one[-kb] = 1.0
  hw = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyH(nc,kc,hw,sg)

  #Make images within specified window have rms of 1.
  frms = ww.rms(f)
  frmsI = div(f,frms)
  print "frms = "+str(frms)
  frmsIrms = ww.rms(frmsI)
  print "frmsI = "+str(frmsIrms)
  sgrms = ww.rms(sg)
  print "sgrms = "+str(sgrms)
  sgrmsI = div(sg,sgrms)
  sgrmsIrms = ww.rms(sgrmsI)
  print "sgrmsI = "+str(sgrmsIrms)
  hsgrms = ww.rms(hsg)
  print "hsgrms = "+str(hsgrms)
  hsgrmsI = div(hsg,hsgrms)
  hsgrmsIrms = ww.rms(hsgrmsI)
  print "hsgrmsI = "+str(hsgrmsIrms)


  #Normalize wavelets
  nhw = normalizeMAAWOS(hw)

  #Plotting
  tmin = 100 
  tmax = 400
  dt = 0.004
  hint = 3.0
  amax = 6.0
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,0.0150,x0*0.0150)
  clip = 2
  #############1 Plot######################################
  pngDir = None#"./report15/synthetics/"
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = tmin*dt,tmax*dt
  hmin,hmax = None,None
  vlabel,vminmax,vint = "time (s)",[vmin,vmax],1.0
  hlabel,hminmax,hint = ["f","hsg","sg"],None,0.5
  clip = 2
  title= "testsino"
  plotting.plotImagesSideBySide(st,sx,[frmsI,hsgrmsI,sgrmsI],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=None,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  title = "sinoh"
  plotting.plotWavelets(Sampling(nc,dt,kc*dt),[nhw],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)


def goSinopec1():
  #get sino images
  x0 = 400
  f,g,u = getSinoTrace(x0)
  st = Sampling(len(f),0.004,0.0)
  plotAmplitudeSpectrumT(st,f,0,499,"amp f")
  plotAmplitudeSpectrumT(st,g,0,700,"amp g")
  #ref = RecursiveExponentialFilter(10)
  #ref.apply1(g,g)
  #SimplePlot.asPoints(g)
  #ref = RecursiveExponentialFilter(3)
  #ref.apply1(f,f)
  #SimplePlot.asPoints(f)
  """
  ng = len(g)
  nu = len(u)

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 21,-10# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  maxpc = 0.01
  niter = 50
  ww = WaveletWarpingCBGN()
  #ww = Start()
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
  plotAmplitudeSpectrumT(Sampling(nb,0.004,0.004*kb), bw, 0, nb, "b Spectrum")

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
  rmsPercentChange = ww.getRmsPercentChangeS()
  lastIter = ww.getLastIter()

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))
  SimplePlot.asPoints(warp.applyS(u,g))

  #Print
  nt = len(f)
  dt = 0.004
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)

  #Plotting
    #Data
  pngDir = None
  title= "f csbg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[.4,2],None
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-8.0
  hmax1,hmax2 = 1.0,8.0
  hlabel = ["f","HSg","CSBg"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin1,hmax1],\
  [hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin2,hmax2]]
  hint = [None,None,None,None,None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

    #GN Measurements
  """
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
  [allRes2NormSqInitS,dataRes2NormSqInitS,bPenaltyRes2NormSqInitS,cPenaltyRes2NormSqInitS,alphaBS],\
  color=None,\
  vlabel=["allRes","dataRes","bPenRes","cPenRes","alphaB"],\
  vminmax=[None,None,None,None,None],\
  vint=[None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  title = "Final 2NormSq of individual and all residuals"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,
  [allRes2NormSqFinaS,dataRes2NormSqFinaS,bPenaltyRes2NormSqFinaS,cPenaltyRes2NormSqFinaS,alphaBS],\
  color=None,\
  vlabel=["allRes","dataRes","bPenRes","cPenRes","alphaB"],\
  vminmax=[None,None,None,None,None],\
  vint=[None,None,None,None,None],\
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
  """


  """  #Normalize
  ncguess = normalizeM(cguess)
  nbguess = normalizeM(bguess)
  ncw = normalizeM(cw)
  ndw = normalizeM(dw)
  nbw = normalizeM(bw)

    #Wavelets
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated c"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated d"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ndw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)
  """
  """

def goSinopec():
  #get sino images
  x0,nx = 0,721
  f,g,u = getSinoImage(x0,nx)
  nx = len(f)
  nt = len(f[0])
  ng = len(g[0])
  nu = len(u[0])

  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  niter = 100
  ww = WaveletWarpingCBGN()
  ww.setMaxPercentChange(0.0001)#units are percentage.
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

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  f = ww.makeRms1(f)
  sbg = ww.makeRms1(sbg)
  csbg = ww.makeRms1(csbg)
  hsg = ww.makeRms1(hsg)
  csbgmhsg = ww.makeRms1(sub(csbg,hsg))

  #Print
  dt = 0.004
  dx = 0.015
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)
  print "Number of time samples: "+str(nt)

  #Plotting
    #Data
  pngDir = None
  title= "f hsg csbg sbg"
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hlabel = ["f(km)","HSg(km)","CSBg(km)","SBg(km)","CSBg-HSg"]
  hminmax = None
  hint = None
  tilespacing = 5
  clip = 2.0
  plotting.plotImagesSideBySide(st,sx,[f,hsg,csbg,sbg,csbgmhsg],\
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

    #Wavelets
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated c"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated d"
  hint = None
  st = Sampling(nc,dt,kc*dt)
  plotting.plotWavelets(st,[ncw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbguess],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "Estimated b"
  hint = None
  st = Sampling(nb,dt,kb*dt)
  plotting.plotWavelets(st,[nbw],hint=hint,title=title,pngDir=pngDir,paper=True,
  onecol=True)




def goVarySyntheticTest():
  #Synthetic parameters
  nt,ni,randomi = 1081,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  freq,decay,= 0.08,0.05
  na,ka = 200,0 #sampling for inverse wavelet A 
  nc,kc = 181,-90 # sampling for wavelet H 
  mp = True#is wavelet in f minimum phase?
  r00,r01,dr0 = 2.0,2.0,0.3
  r10,r11,dr1 = 1.99999,1.99999,0.3
  c0,c1,dc = 0.0,0.0,10.0
  nr1 = int((r11-r10)/dr1)+1
  nr0 = int((r01-r00)/dr0)+1
  nc = int((c1-c0)/dc)+1
  noise = False
  nrmsf = 0.01
  nrmsg = 0.01

  na0 = 3
  na1 = 480
  dna = 10
  nna = (na1-na0)/dna+1
  maxna = ((na1-na0)/dna)*dna+na0
  l0ns = zerofloat(nna)
  l1ns= zerofloat(nna)
  measures1 = zerofloat(nna)
  measures2 = zerofloat(nna)
  naks = zerofloat(maxna,nna)
  naws = zerofloat(maxna,nna)
  ncks = zerofloat(nc,nna)
  ncws = zerofloat(nc,nna)
  print str(maxna)
  print str(nna)

  write = False 
  readplot = True
  if write:
    #measures = zerofloat(nc,nr1,nr0)
    measures = zerofloat(nna)
    for na in range(na0,na1,dna):
      print "na = "+str(na)
      for ir0 in range(nr0):
        for ir1 in range(nr1):
          for ic in range(nc):
            r0 = r00+dr0*ir0  
            r1 = r10+dr1*ir1  
            c = c0+dc*ic
            #Estimation and sampling parameters.
            sfac = 1.000
            p,q,f,g,u,tmin,tmax = createSyntheticLn1D(freq,decay,mp,r0,r1,c,\
            noise,nrmsf,nrmsg,nt,ni,randomi,moreps)
            #p,q,f,g,u,tmin,tmax = createSyntheticCos(freq,decay,mp,r0,r1,c,\
            #noise,nrmsf,nrmsg,nt,ni)
            #p,q,f,g,u,tmin,tmax = createSyntheticShift(freq,decay,mp,c,\
            #noise,nrmsf,nrmsg,nt,ni)
            du = computeBackDiff(u)
            nt = len(f)

            #main piece of wavelet estimation code
            ww = WaveletWarpingAEig()
            ww.setTimeRange(tmin,tmax)
            #ww.setGaussWeights(1.0,0.0,50,nt-1,nt)
            aw = ww.getInverseA(na,ka,u,f,g)
            hw = ww.getWaveletC(na,ka,aw,nc,kc)
            hk = getWavelet(freq,decay,nc,kc,mp) # known wavelet
            #ak = ww.getWaveletC(nc,kc,hk,na,ka)
            ak = getMinPhaseInverseWavelet(na,freq,decay)
            #dump(ak)
            naw = normalizeMAAWOS(aw)
            nak = normalizeMAAWOS(ak)
            ncw = normalizeMAAWOS(hw)
            nck = normalizeMAAWOS(hk)
            measures1[(na-na0)/dna] = ww.getUniqMeasure1()
            measures2[(na-na0)/dna] = ww.getUniqMeasure2()
            lambda0 = ww.getEigVal(0)
            lambda1 = ww.getEigVal(1)
            lambdan = ww.getEigVal(na-1)
            l0ns[(na-na0)/dna] = lambda0/lambdan
            l1ns[(na-na0)/dna] = lambda1/lambdan
            for i in range(0,na):
              naks[(na-na0)/dna][i] = nak[i]
              naws[(na-na0)/dna][i] = naw[i]
            for i in range(0,nc):
              ncks[(na-na0)/dna][i] = nck[i]
              ncws[(na-na0)/dna][i] = ncw[i]


    #Write out information to file for later
    title = "r0"+str(r0)+"r1"+str(r1)+"na0"+str(na0)+"na1"+str(na1)+"ka"+str(ka)+"ni"+str(ni)
    aos1 = ArrayOutputStream("./tests/"+title+"measures1")
    aos1.writeFloats(measures1)
    aos1.close()
    aos2 = ArrayOutputStream("./tests/"+title+"measures2")
    aos2.writeFloats(measures2)
    aos2.close()
    aos3 = ArrayOutputStream("./tests/"+title+"naks")
    aos3.writeFloats(naks)
    aos3.close()
    aos4 = ArrayOutputStream("./tests/"+title+"naws")
    aos4.writeFloats(naws)
    aos4.close()
    aos5 = ArrayOutputStream("./tests/"+title+"l0ns")
    aos5.writeFloats(l0ns)
    aos5.close()
    aos6 = ArrayOutputStream("./tests/"+title+"l1ns")
    aos6.writeFloats(l1ns)
    aos6.close()
    aos7 = ArrayOutputStream("./tests/"+title+"ncks")
    aos7.writeFloats(ncks)
    aos7.close()
    aos8 = ArrayOutputStream("./tests/"+title+"ncws")
    aos8.writeFloats(ncws)
    aos8.close()
    
  
  if readplot:
    dt = 0.004
    #Read and plot
    title = "r0"+str(r00)+"r1"+str(r10)+"na0"+str(na0)+"na1"+str(na1)+"ka"+str(ka)+"ni"+str(ni)
    ais1 = ArrayInputStream("./tests/"+title+"measures1")
    ais1.readFloats(measures1)
    ais1.close()
    ais2 = ArrayInputStream("./tests/"+title+"measures2")
    ais2.readFloats(measures2)
    ais2.close()
    ais3 = ArrayInputStream("./tests/"+title+"naks")
    ais3.readFloats(naks)
    ais3.close()
    ais4 = ArrayInputStream("./tests/"+title+"naws")
    ais4.readFloats(naws)
    ais4.close()
    ais5 = ArrayInputStream("./tests/"+title+"l0ns")
    ais5.readFloats(l0ns)
    ais5.close()
    ais6 = ArrayInputStream("./tests/"+title+"l1ns")
    ais6.readFloats(l1ns)
    ais6.close()
    ais7 = ArrayInputStream("./tests/"+title+"ncks")
    ais7.readFloats(ncks)
    ais7.close()
    ais8 = ArrayInputStream("./tests/"+title+"ncws")
    ais8.readFloats(ncws)
    ais8.close()


    #title = "r0"+str(r0)+"r1"+str(r1)+"na0"+str(na0)+"na1"+str(na1)+"ka"+str(ka)+"ni"+str(ni)
    title = " "
    snaks = zerofloat(nna,maxna)
    snaws = zerofloat(nna,maxna)
    for i in range(0,nna):
      for j in range(0,maxna):
        snaks[j][i] = naks[i][j]
        snaws[j][i] = naws[i][j]
    #dump(snaks)
    for i in range(0,nna):
      plotWavelets(Sampling(nc,1.0,kc*dt),[ncws[i],ncks[i]],title="na = "+str(i*dna+na0)+title,pngDir=pngDir,onecol=True)
    for i in range(0,13):
      plotWavelets(Sampling(nna,dna,na0),[snaws[i],snaks[i]],title="naws"+str(i)+" and naks"+str(i)+title,pngDir=pngDir,onecol=True)
    plotSequenceMeas(Sampling((nna),dna,na0),[measures1],labels=["Measure1"],pngDir=pngDir,title=title)
    plotSequenceMeas(Sampling((nna),dna,na0),[measures2],labels=["Measure2"],pngDir=pngDir,title=title)
    plotSequenceMeas(Sampling((nna),dna,na0),[l0ns],labels=["l0ns"],pngDir=pngDir,title=title)
    plotSequenceMeas(Sampling((nna),dna,na0),[l1ns],labels=["l1ns"],pngDir=pngDir,title=title)
  #3D plotting
  #sr0 = Sampling(nr0,dr0,r00)
  #sr1 = Sampling(nr1,dr1,r10)
  #sc = Sampling(nc,dc,c0)
  #ipg = ImagePanelGroup(sc,sr1,sr0,measures)
  #sf = SimpleFrame(AxesOrientation.XOUT_YRIGHT_ZUP)
  #sf.addImagePanels(ipg)

  


#Compares two different methods to create impulses that are not aliased. 
#The two methods are using the accumulate method in the SincInterpolator class in the
#MinesJTK and applying an ideal anti-alias filter to Dirac delta functions.
#Sequences created by the accumulate method will end with ac and sequences created by 
#the anti-alias filter will end with aa.
def goAccumulateVsAntiAliasingDeltaFunctions():
  #Synthetic parameters
  nt,ni,randomi = 1081,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  freq,decay,= 0.08,0.05
  na,ka = 40,0 #sampling for inverse wavelet A 
  nc,kc = 181,-90 # sampling for wavelet H 
  mp = True#is wavelet in f minimum phase?
  c = 0.0#The amount of shift between p and q.
  r0,r1 = 4.0,1.1#The first and last warping amounts that will be applied to trace p. 
                 #A line is used to calculate the warping amount for a particular time.
  nrmsf = 0.01
  nrmsg = 0.01

  #Estimation and sampling parameters.
  sfac = 1.000
  pac,qac,f,g,u,tmin,tmax = createSyntheticLn1D(freq,decay,mp,r0,r1,c,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  paa,qaa,u,tmin,tmax = deltas(r0,r1,nt,ni)#anti-alias filter
  du = computeBackDiff(u)
  diffp = sub(pac,paa)
  diffq = sub(qac,qaa)
  mdiffp = max(diffp) 
  mdiffq = max(diffq) 
  diffpd4 = mdiffp/2.0
  diffqd4 = mdiffq/2.0
  
  #plotting
  nt,dt,ft = nt,0.004,0.000 # used for plotting only
  st = Sampling(nt,dt,ft)
  amax = [1.0,1.0,mdiffp,1.0,1.0,mdiffq,.023,r0]
  tmark = [0.5,0.5,diffpd4,0.5,0.5,diffpd4,1.0]
  title = "na="+str(na)+"ka="+str(ka)+"nc="+str(nc)+"ka="+str(ka)
  plotSequences(st,[pac,paa,diffp,qac,qaa,diffq,du],labels=["pac","paa","diffp","qac","qaa","diffq","du"],amax=amax,tmark=tmark,title="pq"+title)

#Create a single array with zeros inbetween the vectors 
def createSingleArray(nfg,nf,f,ng,g):
  fg = zerofloat(nfg)
  ng0 = nfg-nf
  for i in range(0,nf):
    fg[i] = f[i]
  for i in range(ng0,nfg):
    fg[i] = g[i-ng0]
  return fg

def normalize(h):
  return div(h,rms(h))

def createUnitVector(h):
  nc = len(h) 
  sum = 0.0
  for i in range(0,nc): 
    sum=sum+(h[i]*h[i])
  return div(h,sqrt(sum))

def rms(h):
  return sqrt(sum(mul(h,h))/len(h))

def norm2sq(h):
  nc = len(h)
  sum = 0.0
  for i in range(0,nc):
    sum = sum+h[i]*h[i]
  return sum


def deltas(r0,r1,nt,ni):
  umax = nt-1
  dt = 1
  a = umax/log(r0/r1)
  b = r0*log(r0/r1)/umax
  p = zerofloat(nt)
  q = zerofloat(nt)
  if r0>r1:
    rmax = r0
  else:
    rmax = r1
  t = rampfloat(0.0,1.0,nt)
  u = mul(a,log(add(1,mul(b,t))))
  dq = umax/float(ni+1)
  print "dq = "+str(dq)
  ts = rampfloat(dq,dq,ni)  
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.40)
  tp = zerofloat(ni)
  tq = zerofloat(ni)
  for ji in range(ni):
    tq[ji] = ts[ji]#time u
    tp[ji] = (exp(tq[ji]/a)-1.0)/b#time t(u)
    print "tp = "+str(tp[ji])
    print "tq = "+str(tq[ji])

  rj = 1.0
  xq = 0
  xp = 0
  sumq = 0
  sump = 0
  for n in range(0,umax):
    sump = 0
    sumq = 0
    for i in range(0,ni):
      xq = PI*(n*dt-tq[i])/dt
      sumq = sumq + sinc(xq)
      xp = PI*(n*dt-tp[i])/dt
      sump = sump + sinc(xp)
    p[n] = (rj/dt)*sump
    q[n] = (rj/dt)*sumq
  #for n in range(0,umax):
    #xq0 = PI*(n*dt-tq[0])/dt
    #xq1 = PI*(n*dt-tq[1])/dt
    #xq2 = PI*(n*dt-tq[2])/dt
    #xp0 = PI*(n*dt-tp[0])/dt
    #xp1 = PI*(n*dt-tp[1])/dt
    #xp2 = PI*(n*dt-tp[2])/dt
    #p[n] = (rj/dt)*(sinc(xp0)+sinc(xp1)+sinc(xp2))
    #q[n] = (rj/dt)*(sinc(xq0)+sinc(xq1)+sinc(xq2))

  itmin = 0
  itmax = nt-1
  return p,q,u,itmin,itmax

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
