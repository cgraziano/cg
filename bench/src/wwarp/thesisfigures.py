#############################################################################
# Demo of 2 wavelet estimations from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.awt.ColorMap import *
from edu.mines.jtk.lapack import *
from wwarp import WaveletWarpingCBGN, WaveletWarpingCBCyclic, WaveletWarpingCA, Warper, Start
from wwarp import AmpSpectrum 
import synthetic
import plotting
from java.util import Random

############################################################################

def main(args):
  #impulseConstantSqueezing()
  #impulseVaryingSqueezing()
  #oneWaveletNoDistortion()
  #twoWaveletsNoDistortion()
  #estimateOneWavelet2t()
  #estimateOneWaveletlnt()
  #estimateTwoWaveletsCycliclnt()
  #estimateTwoWaveletsGNlnt()
  #estimateTwoWaveletsNoiseLowSqueezingCycliclnt()
  #estimateTwoWaveletsNoiseLowSqueezingGNlnt()
  #increaseBy1GNNoise()
  #increaseBy1PreviousSolutionGNNoise()
  #increaseBy1CyclicNoise()
  #increaseBy1PreviousSolutionCyclicNoise()
  #noiseComparisonGN()
  #noiseComparisonCyclic()
  goSinopecGN()
  #goSinopecCyclic()
  goGBCGN()

def impulseConstantSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 181,-90 # sampling for wavelet H 
  niter = 1
  r0,r1 = 2.0,2.0#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.07
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 0.00

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1DSimple(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,randomi,moreps)
  SimplePlot.asPoints(p)

  #Warp q to p using new method
  warp = Warper()
  sq = warp.applyS(u,q)
  ugWugLwugDlwug = warp.getWarpingStages()
  uq = ugWugLwugDlwug[0]
  wuq = ugWugLwugDlwug[1]
  lwuq = ugWugLwugDlwug[2]
  dlwuq = ugWugLwugDlwug[3]
  scaleR = warp.getSamplingIntervalScaleFactor()
  
  #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/2pt2/"
  pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  hmin,hmax = -1.5,1.5
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax]],[1.0,1.0,1.0]
  hsize,vsize = 960,560
  title= "p q sq 2t"
  plotting.plotTracesSideBySide(st,[p,q,sq],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  #Plot amplitude spectrums
  ampspec = AmpSpectrum()
  nq = len(q)
  nuq = len(uq)
  nwuq = len(wuq)
  nlwuq = len(lwuq)
  ndlwuq = len(dlwuq)
  dtf = 0.004
  dtg = 0.004
  dtug = dtf/scaleR
  ft = 0.0
  stq = Sampling(nq,dtg,ft)
  stuq = Sampling(nuq,dtug,ft)
  stwuq = Sampling(nwuq,dtug,ft)
  stlwuq = Sampling(nlwuq,dtug,ft)
  stdlwuq = Sampling(ndlwuq,dtg,ft)

  simplePlot = False
  forSlide = False
  forPrint = True
  pngDir = "./thesisFigures/2pt2/"
  fracWidth,fracHeight,aspectRatio = 0.8,0.9,16.0/9.0
  ampspec.plotAmplitudeSpectrum(stq,q,0,nq/2,simplePlot,forSlide,forPrint,pngDir,"amp q",fracWidth,fracHeight,aspectRatio);
  ampspec.plotAmplitudeSpectrum(stuq,uq,0,nq,simplePlot,forSlide,forPrint,pngDir,"amp uq",fracWidth,fracHeight,aspectRatio);
  ampspec.plotAmplitudeSpectrum(stwuq,wuq,0,nq/2,simplePlot,forSlide,forPrint,pngDir,"amp wuq",fracWidth,fracHeight,aspectRatio);
  ampspec.plotAmplitudeSpectrum(stlwuq,lwuq,0,nq/2,simplePlot,forSlide,forPrint,pngDir,"amp lwuq",fracWidth,fracHeight,aspectRatio);
  ampspec.plotAmplitudeSpectrum(stdlwuq,dlwuq,0,nq/4,simplePlot,forSlide,forPrint,pngDir,"amp dlwuq",fracWidth,fracHeight,aspectRatio);



def impulseVaryingSqueezing():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 181,-90 # sampling for wavelet H 
  niter = 1
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.07
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 0.00

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1DSimple(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,randomi,moreps)
  SimplePlot.asPoints(p)

  #Warp q to p using new method
  warp = Warper()
  sq = warp.applyS(u,q)
  SimplePlot.asPoints(sq)
  
  #Warp q to p using old method
  warp = Warper()
  wlq = warp.applyOldS(u,q)
  SimplePlot.asPoints(wlq)

  #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/2pt2/"
  pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -1.5,1.5
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax]],[1.0,1.0,1.0]
  hsize,vsize = 960,560
  title= "p q sq lnt"
  plotting.plotTracesSideBySide(st,[p,q,sq],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)


def oneWaveletNoDistortion():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 81,-40 # sampling for wavelet H 
  niter = 1
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.07
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 0.00

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1DSimple(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,randomi,moreps)

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  maxpc = 0.01
  niter = 315
  ww = WaveletWarpingCBGN()
  ww.setMinPercentChange(maxpc)#units are percentage.
  ww.setTimeRange(tmin,tmax)

  #Get know wavelets
  nd,kd = nc,kc
  na,ka = nb,kb
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  ak = ww.getWaveletC(nc,kc,ck,na,ka)
  bk = ww.getWaveletC(nd,kd,dk,nb,kb)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bk,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,ck,sbg)
  sg = warp.applyS(u,g)
  csbg = ww.applyC(nc,kc,ck,warp.applyS(u,ww.applyC(nb,kb,bk,g)))

  #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/2pt3/"
  pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -12.0,12.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax]],[5.0,5.0,5.0]
  hsize,vsize = 960,560
  title= "f g Sg lnt one wavelet"
  plotting.plotTracesSideBySide(st,[f,g,sg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  directory = "./thesisFigures/2pt3/"
  pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -12.0,12.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax]],[5.0,5.0,5.0]
  hsize,vsize = 960,560
  title= "f csbg sg lnt one wavelet"
  plotting.plotTracesSideBySide(st,[f,csbg,sg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  directory = "./thesisFigures/2pt3/"
  pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -5.5,5.5
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax]],[2.0,2.0,2.0,2.0]
  hsize,vsize = 960,560
  title= "g bg sbg csbg lnt one wavelet"
  plotting.plotTracesSideBySide(st,[g,bg,sbg,csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

    #Normalize
  nck = normalizeM(ck)
  #ndk = normalizeM(dk)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  #ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
  nc = nc2
  #nb = nb2
  dt = dt2
  #"""

  pngDir = directory
  #pngDir = None     
  title = "Known c one wavelet"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  lineStyle = [PointsView.Line.SOLID,PointsView.Line.SOLID,PointsView.Line.SOLID]
  lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  hsize = 960
  vsize = 350
  plotting.plotWavelets(st,[nck],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,linecolor=lineColor,markstyle=markStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

def twoWaveletsNoDistortion():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 81,-40 # sampling for wavelet H 
  niter = 1
  r0,r1 = 3.15,1.55#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.07
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 0.00

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1DSimple(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,randomi,moreps)

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  maxpc = 0.01
  niter = 300
  ww = WaveletWarpingCBGN()
  ww.setMinPercentChange(maxpc)#units are percentage.
  ww.setTimeRange(tmin,tmax)

  #Get know wavelets
  nd,kd = nc,kc
  na,ka = nb,kb
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  ak = ww.getWaveletC(nc,kc,ck,na,ka)
  bk = ww.getWaveletC(nd,kd,dk,nb,kb)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bk,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,ck,sbg)
  sg = warp.applyS(u,g)
  csbg = ww.applyC(nc,kc,ck,warp.applyS(u,ww.applyC(nb,kb,bk,g)))

  #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/2pt3/"
  pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -12.0,12.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax]],[5.0,5.0,5.0]
  hsize,vsize = 960,560
  title= "f g Sg lnt two wavelets"
  plotting.plotTracesSideBySide(st,[f,g,sg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  directory = "./thesisFigures/2pt3/"
  pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -12.0,12.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax]],[5.0,5.0,5.0]
  hsize,vsize = 960,560
  title= "f csbg sg lnt two wavelets"
  plotting.plotTracesSideBySide(st,[f,csbg,sg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  directory = "./thesisFigures/2pt3/"
  pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -5.5,5.5
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax]],[2.0,2.0,2.0,2.0]
  hsize,vsize = 960,560
  title= "g bg sbg csbg lnt two wavelets"
  plotting.plotTracesSideBySide(st,[g,bg,sbg,csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

    #Normalize
  nck = normalizeM(ck)
  ndk = normalizeM(dk)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
  nc = nc2
  #nb = nb2
  dt = dt2
  #"""

  lineStyle = [PointsView.Line.SOLID,PointsView.Line.SOLID,PointsView.Line.SOLID]
  lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known c two wavelets"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[nck],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,linecolor=lineColor,markstyle=markStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Known d two wavelets"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndk],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,linecolor=lineColor,markstyle=markStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

def estimateOneWavelet2t():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 2.0,2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0#0.4
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  na,ka = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  sfac = 0.0

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  ww = WaveletWarpingCA()
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  aw = ww.getInverseA(na,ka,u,f,g)
  cw = ww.getWaveletC(na,ka,aw,nc,kc)

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  ak = ww.getWaveletC(nc,kc,ck,na,ka)

  #Processing
  warp = Warper()
  ag = ww.applyC(na,ka,aw,g)
  sg = warp.applyS(u,g)
  sag = warp.applyS(u,ag)
  csag = ww.applyC(nc,kc,cw,sag)
  aone = zerofloat(na)
  aone[-ka] = 1.0
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,na,ka,aone,hstabfact,u,f,g)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(na,ka,aone,g)))
  SimplePlot.asPoints(warp.applyS(u,g))

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

 #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/3pt2/"
  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -13.0,13.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amount of squeezing"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[0,3]],[5.0,5.0,5.0,1.0]
  hsize,vsize = 960,560
  title= "f g sg dudt 2t one wavelet"
  plotting.plotTracesSideBySide(st,[f,g,sg,du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "[f,csag] [f,hsg] 2t"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,1.088
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint1,hint2 = 0.5,5.0
  hmin1,hmin2 = -1.0,-12.0
  hmax1,hmax2 = 1.0,12.0
  hlabel = ["Amplitude","Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [hint2,hint2,hint2]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,csag],[f,hsg]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing=tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

    #Normalize
  nck = normalizeM(ck)
  ncw = normalizeM(cw)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  #nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  nc = nc2
  #nb = nb2
  dt = dt2
  #"""

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated c one wavelets 2t"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

def estimateOneWaveletlnt():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0#0.4
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  na,ka = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  sfac = 0.0

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  ww = WaveletWarpingCA()
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  aw = ww.getInverseA(na,ka,u,f,g)
  cw = ww.getWaveletC(na,ka,aw,nc,kc)

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  ak = ww.getWaveletC(nc,kc,ck,na,ka)

  #Processing
  warp = Warper()
  ag = ww.applyC(na,ka,aw,g)
  sg = warp.applyS(u,g)
  sag = warp.applyS(u,ag)
  csag = ww.applyC(nc,kc,cw,sag)
  aone = zerofloat(na)
  aone[-ka] = 1.0
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,na,ka,aone,hstabfact,u,f,g)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(na,ka,aone,g)))
  SimplePlot.asPoints(warp.applyS(u,g))

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

 #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/3pt2/"
  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -13.0,13.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amount of squeezing"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[0,3.5]],[5.0,5.0,5.0,1.0]
  hsize,vsize = 960,560
  title= "f g sg dudt lnt one wavelet"
  plotting.plotTracesSideBySide(st,[f,g,sg,du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "[f,csag] [f,hsg] lnt"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint1,hint2 = 0.5,5.0
  hmin1,hmin2 = -1.0,-12.0
  hmax1,hmax2 = 1.0,12.0
  hlabel = ["Amplitude","Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [hint2,hint2,hint2]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,csag],[f,hsg]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing=tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

    #Normalize
  nck = normalizeM(ck)
  ncw = normalizeM(cw)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  #nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  nc = nc2
  #nb = nb2
  dt = dt2
  #"""

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated c one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

def estimateTwoWaveletsCycliclnt():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0#0.4
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  sfac = 0.0

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  niter = 500
  ww = WaveletWarpingCBCyclic()
  ww.setMinPercentChange(0.01)
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(0.000)
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)

  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  lastIter = ww.getLastIter()
  print "lastIter = "+str(lastIter)
  allResRmsAllS = ww.getAllResRmsAllS()
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sg = warp.applyS(u,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
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

 #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/3pt3/"
  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -13.0,13.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amount of squeezing"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[0,3.5]],[3.0,3.0,3.0,1.0]
  hsize,vsize = 960,560
  title= "f g sg dudt lnt one wavelet"
  plotting.plotTracesSideBySide(st,[f,g,sg,du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "[f,csbg] [f,hsg] lnt"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [hint2,hint2,hint2]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,csbg],[f,hsg]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing=tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

    #GN Meas
  print "All RMS"
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals"
  maxrmsri = 0.14#max(allResRmsAllS)#0.15
  minrmsrf = 0.0#allResRmsAllS[lastIter]
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],20.0
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

    #Normalize
  nck = normalizeMAAWOS(ck)
  ndk = normalizeMAAWOS(dk)
  ncw = normalizeMAAWOS(cw)
  ndw = normalizeMAAWOS(dw)
  dump(ncw)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  #nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndw = interpolate(nd,kc,ndw,dt,nc2,dt2,error,freq)
  nc = nc2
  #nb = nb2
  dt = dt2
  #"""

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated c one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated d one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw,ndk],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

def estimateTwoWaveletsGNlnt():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0#0.4
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  sfac = 0.0

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  maxpc = 0.01
  niter = 500
  ww = WaveletWarpingCBGN()
  ww.setMinPercentChange(maxpc)#units are percentage.
  ww.setTimeRange(tmin,tmax)

  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  lastIter = ww.getLastIter()
  print "lastIter = "+str(lastIter)
  allResRmsAllS = ww.getAllResRmsAllS()
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sg = warp.applyS(u,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
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

 #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/3pt4/"
  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -13.0,13.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amount of squeezing"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[0,3.5]],[3.0,3.0,3.0,1.0]
  hsize,vsize = 960,560
  title= "f g sg dudt lnt one wavelet"
  plotting.plotTracesSideBySide(st,[f,g,sg,du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "[f,csbg] [f,hsg] lnt"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [hint2,hint2,hint2]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,csbg],[f,hsg]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing=tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

    #GN Meas
  print "All RMS"
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals"
  maxrmsri = 0.14#max(allResRmsAllS)#0.15
  minrmsrf = 0.0#allResRmsAllS[lastIter]
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],20.0
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

    #Normalize
  nck = normalizeMAAWOS(ck)
  ndk = normalizeMAAWOS(dk)
  ncw = normalizeMAAWOS(cw)
  ndw = normalizeMAAWOS(dw)
  dump(ncw)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  #nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndw = interpolate(nd,kc,ndw,dt,nc2,dt2,error,freq)
  nc = nc2
  #nb = nb2
  dt = dt2
  #"""

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated c one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated d one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw,ndk],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

def estimateTwoWaveletsNoiseLowSqueezingCycliclnt():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 1.6,1.4
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0#0.4
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  sfac = 0.0

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  niter = 500
  ww = WaveletWarpingCBCyclic()
  ww.setMinPercentChange(0.01)
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(0.000)
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)

  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  lastIter = ww.getLastIter()
  print "lastIter = "+str(lastIter)
  allResRmsAllS = ww.getAllResRmsAllS()
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sg = warp.applyS(u,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
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

 #Plotting
  #############1 Plot######################################
  #directory = "./thesisFigures/3pt3/"
  directory = None
  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -13.0,13.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amount of squeezing"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[0,3.5]],[3.0,3.0,3.0,1.0]
  hsize,vsize = 960,560
  title= "f g sg dudt lnt one wavelet"
  plotting.plotTracesSideBySide(st,[f,g,sg,du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "[f,csbg] [f,hsg] lnt"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [hint2,hint2,hint2]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,csbg],[f,hsg]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing=tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

    #GN Meas
  print "All RMS"
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals"
  maxrmsri = 0.14#max(allResRmsAllS)#0.15
  minrmsrf = 0.0#allResRmsAllS[lastIter]
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],20.0
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

    #Normalize
  nck = normalizeMAAWOS(ck)
  ndk = normalizeMAAWOS(dk)
  ncw = normalizeMAAWOS(cw)
  ndw = normalizeMAAWOS(dw)
  dump(ncw)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  #nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndw = interpolate(nd,kc,ndw,dt,nc2,dt2,error,freq)
  nc = nc2
  #nb = nb2
  dt = dt2
  #"""

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated c one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated d one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw,ndk],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

def estimateTwoWaveletsNoiseLowSqueezingGNlnt():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 1.6,1.4
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0#0.4
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  sfac = 0.0

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  #Estimate wavelet
  maxpc = 0.01
  niter = 500
  ww = WaveletWarpingCBGN()
  ww.setMinPercentChange(maxpc)#units are percentage.
  ww.setTimeRange(tmin,tmax)

  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  lastIter = ww.getLastIter()
  print "lastIter = "+str(lastIter)
  allResRmsAllS = ww.getAllResRmsAllS()
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)

  #Get know wavelets
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sg = warp.applyS(u,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  hstabfact = 0.0
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

 #Plotting
  #############1 Plot######################################
  #directory = "./thesisFigures/3pt4/"
  directory = None
  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -13.0,13.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amount of squeezing"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[0,3.5]],[3.0,3.0,3.0,1.0]
  hsize,vsize = 960,560
  title= "f g sg dudt lnt one wavelet"
  plotting.plotTracesSideBySide(st,[f,g,sg,du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "[f,csbg] [f,hsg] lnt"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2]]
  hint = [hint2,hint2,hint2]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,csbg],[f,hsg]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing=tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

    #GN Meas
  print "All RMS"
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals"
  maxrmsri = 0.14#max(allResRmsAllS)#0.15
  minrmsrf = 0.0#allResRmsAllS[lastIter]
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],20.0
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

    #Normalize
  nck = normalizeMAAWOS(ck)
  ndk = normalizeMAAWOS(dk)
  ncw = normalizeMAAWOS(cw)
  ndw = normalizeMAAWOS(dw)
  dump(ncw)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  #nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndw = interpolate(nd,kc,ndw,dt,nc2,dt2,error,freq)
  nc = nc2
  #nb = nb2
  dt = dt2
  #"""

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated c one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw,nck],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)

  lineStyle = [PointsView.Line.NONE,PointsView.Line.SOLID,PointsView.Line.NONE]
  #lineColor = [Color.BLACK,Color.BLACK,Color.BLACK]
  #markStyle = [PointsView.Mark.NONE,PointsView.Mark.NONE,PointsView.Mark.NONE]
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "Known and Estimated d one wavelets lnt"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw,ndk],hint=hint,hsize=hsize,vsize=vsize,linestyle=lineStyle,title=title,pngDir=pngDir,paper=True,twocol=True)


def increaseBy1GNNoise():
  directory = "./thesisFigures/3pt5/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.5#0.4
  nrmsg = nrmsf
  maxpc = 0.01
  niter = 300

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  nhfinal = 91
  nhinitial = 81
  khinitial = -40
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  allRmsResH = zerofloat(n)
  allRmsResCB = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    ww = WaveletWarpingCBGN()
    ww.setMinPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
  
    #First guesses of c and b. 
    bone = zerofloat(nb)
    bone [-kb] = 1.0
    bguess = copy(bone)
    hstabfact = 0.0
    hw =  ww.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)
    cguess = ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get iteration information
    lastIter = ww.getLastIter()
    allResRmsAllS = ww.getAllResRmsAllS()
    allRmsResCB[i] = allResRmsAllS[lastIter]
    warp = Warper()
    sg = warp.applyS(u,g)
    hsg = ww.applyC(nh,kh,hw,sg)
    bpen = zerofloat(nb)
    cpen = zerofloat(nc)
    allRmsResH[i] = ww.rmsOfObjectiveFunction(sub(hsg,f),bpen,cpen)
    print "rmsrh = "+str(allRmsResH)
    print "rmsrcb = "+str(allRmsResCB)
    kb = -int(nb/2)
    nb = nb+1
    nh = nh+1

  pngDir = None
  pngDir = directory
  title = "Increaseby1GN"
  maxrms = 0.27#max([max(allRmsResCB),max(allRmsResH)])
  minrms = 0.20#min([min(allRmsResCB),min(allRmsResH)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(allRmsResH))
  print "rmsrcb = "+str(len(allRmsResCB))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.RED]
  hsize,vsize = 960,350
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],0.01
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],1.0
  plotting.plotMeasInSamePlotLargeMarks(si, [allRmsResH,allRmsResCB],\
  color=color,hsize=hsize,vsize=vsize,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

def increaseBy1PreviousSolutionGNNoise():
  directory = "./thesisFigures/3pt5/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.5#0.4
  nrmsg = nrmsf
  maxpc = 0.01
  niter = 300

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  nhfinal = 91
  nhinitial = 81
  khinitial = -40
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  allRmsResH = zerofloat(n)
  allRmsResCB = zerofloat(n)
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bw = copy(bone)
  hstabfact = 0.0
  ww = WaveletWarpingCBGN()
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cw = copy(hw)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    ww = WaveletWarpingCBGN()
    ww.setMinPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
    bw = addZeros(bw,nb)
    cw = addZeros(cw,nc)
    dump(bw)
    bone = zerofloat(nb)
    bone [-kb] = 1.0
    hw =  ww.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,bw,nc,kc,cw,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get iteration information
    lastIter = ww.getLastIter()
    allResRmsAllS = ww.getAllResRmsAllS()
    allRmsResCB[i] = allResRmsAllS[lastIter]
    warp = Warper()
    sg = warp.applyS(u,g)
    hsg = ww.applyC(nh,kh,hw,sg)
    bpen = zerofloat(nb)
    cpen = zerofloat(nc)
    allRmsResH[i] = ww.rmsOfObjectiveFunction(sub(hsg,f),bpen,cpen)
    print "rmsrh = "+str(allRmsResH)
    print "rmsrcb = "+str(allRmsResCB)
    kb = -int(nb/2)
    nb = nb+1
    nh = nh+1

  pngDir = None
  pngDir = directory
  title = "Increaseby1GNPreviousSol"
  maxrms = 0.27#max([max(allRmsResCB),max(allRmsResH)])
  minrms = 0.20#min([min(allRmsResCB),min(allRmsResH)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(allRmsResH))
  print "rmsrcb = "+str(len(allRmsResCB))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.RED]
  hsize,vsize = 960,350
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],0.01
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],1.0
  plotting.plotMeasInSamePlotLargeMarks(si, [allRmsResH,allRmsResCB],\
  color=color,hsize=hsize,vsize=vsize,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

def increaseBy1CyclicNoise():
  directory = "./thesisFigures/3pt5/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.5#0.4
  nrmsg = nrmsf
  maxpc = 0.01
  niter = 300

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  nhfinal = 91
  nhinitial = 81
  khinitial = -40
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  allRmsResH = zerofloat(n)
  allRmsResCB = zerofloat(n)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    ww = WaveletWarpingCBCyclic()
    ww.setMinPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
  
    #First guesses of c and b. 
    bone = zerofloat(nb)
    bone [-kb] = 1.0
    bguess = copy(bone)
    hstabfact = 0.0
    hw =  ww.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)
    cguess = ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get iteration information
    lastIter = ww.getLastIter()
    allResRmsAllS = ww.getAllResRmsAllS()
    allRmsResCB[i] = allResRmsAllS[lastIter]
    warp = Warper()
    sg = warp.applyS(u,g)
    hsg = ww.applyC(nh,kh,hw,sg)
    bpen = zerofloat(nb)
    cpen = zerofloat(nc)
    allRmsResH[i] = ww.rmsOfObjectiveFunction(sub(hsg,f),bpen,cpen)
    print "rmsrh = "+str(allRmsResH)
    print "rmsrcb = "+str(allRmsResCB)
    kb = -int(nb/2)
    nb = nb+1
    nh = nh+1

  pngDir = None
  pngDir = directory
  title = "Increaseby1Cyclic"
  maxrms = 0.27#max([max(allRmsResCB),max(allRmsResH)])
  minrms = 0.20#min([min(allRmsResCB),min(allRmsResH)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(allRmsResH))
  print "rmsrcb = "+str(len(allRmsResCB))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.RED]
  hsize,vsize = 960,350
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],0.01
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],1.0
  plotting.plotMeasInSamePlotLargeMarks(si, [allRmsResH,allRmsResCB],\
  color=color,hsize=hsize,vsize=vsize,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

def increaseBy1PreviousSolutionCyclicNoise():
  directory = "./thesisFigures/3pt5/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.5#0.4
  nrmsg = nrmsf
  maxpc = 0.01
  niter = 300

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

  #set tmin and tmax 
  tmin = tmin
  tmax = tmax

  nhfinal = 91
  nhinitial = 81
  khinitial = -40
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  allRmsResH = zerofloat(n)
  allRmsResCB = zerofloat(n)
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bw = copy(bone)
  hstabfact = 0.0
  ww = WaveletWarpingCBGN()
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  cw = copy(hw)
  for i in range(0,n):
    print "################i = "+str(i)+" #####################"
    print "nb = "+str(nb)
    print "kb = "+str(kb)
    print "nc = "+str(nc)
    print "kc = "+str(kc)
    print "nh = "+str(nh)
    print "kh = "+str(kh)
    #Estimate wavelet
    ww = WaveletWarpingCBCyclic()
    ww.setMinPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)
    bw = addZeros(bw,nb)
    cw = addZeros(cw,nc)
    dump(bw)
    bone = zerofloat(nb)
    bone [-kb] = 1.0
    hw =  ww.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,bw,nc,kc,cw,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get iteration information
    lastIter = ww.getLastIter()
    allResRmsAllS = ww.getAllResRmsAllS()
    allRmsResCB[i] = allResRmsAllS[lastIter]
    warp = Warper()
    sg = warp.applyS(u,g)
    hsg = ww.applyC(nh,kh,hw,sg)
    bpen = zerofloat(nb)
    cpen = zerofloat(nc)
    allRmsResH[i] = ww.rmsOfObjectiveFunction(sub(hsg,f),bpen,cpen)
    print "rmsrh = "+str(allRmsResH)
    print "rmsrcb = "+str(allRmsResCB)
    kb = -int(nb/2)
    nb = nb+1
    nh = nh+1

  pngDir = None
  pngDir = directory
  title = "Increaseby1CyclicPreviousSol"
  maxrms = 0.27#max([max(allRmsResCB),max(allRmsResH)])
  minrms = 0.20#min([min(allRmsResCB),min(allRmsResH)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(allRmsResH))
  print "rmsrcb = "+str(len(allRmsResCB))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.RED]
  hsize,vsize = 960,350
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],0.01
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],1.0
  plotting.plotMeasInSamePlotLargeMarks(si, [allRmsResH,allRmsResCB],\
  color=color,hsize=hsize,vsize=vsize,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

def noiseComparisonGN():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsfs = [0.0,0.2,0.4,0.6,0.8,1.0]
  nnrmsf = len(nrmsfs)
  rmsNoise = zerofloat(nnrmsf)
  fNoise = zerofloat(nt,nnrmsf)
  gNoise = zerofloat(nt,nnrmsf)
  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  sfac = 0.0
  ncwNoise = zerofloat(nc,nnrmsf)
  ndwNoise = zerofloat(nc,nnrmsf)
  nckNoise = zerofloat(nc,nnrmsf)
  ndkNoise = zerofloat(nc,nnrmsf)
  scale = 4
  nc2 = scale*(nc-1)+1
  incwNoise = zerofloat(nc2,nnrmsf)
  indwNoise = zerofloat(nc2,nnrmsf)
  inckNoise = zerofloat(nc2,nnrmsf)
  indkNoise = zerofloat(nc2,nnrmsf)
  for i in range(nnrmsf):
    nrmsf = nrmsfs[i]
    nrmsg = nrmsf
    #Create synthetic f and g.
    p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
    freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
    du = computeBackDiff(u)
    fNoise[i] = f
    gNoise[i] = g

    #Wavelet estimation parameters
    nc,kc = 81,-40# sampling for wavelet H 
    nd,kd = nc,kc
    nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
    sfac = 0.0

    #set tmin and tmax 
    tmin = tmin
    tmax = tmax

    #Estimate wavelet
    maxpc = 0.01
    niter = 500
    ww = WaveletWarpingCBGN()
    ww.setMinPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)

    #First guesses of c and b. 
    bone = zerofloat(nb)
    bone[-kb] = 1.0
    bguess = copy(bone)
    hstabfact = 0.0
    hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    cguess = copy(hw)
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    lastIter = ww.getLastIter()
    print "lastIter = "+str(lastIter)
    allResRmsAllS = ww.getAllResRmsAllS()
    rmsNoise[i] = allResRmsAllS[lastIter]
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get know wavelets
    ck = getWavelet(freqc,decayc,nc,kc,mpc)
    dk = getWavelet(freqd,decayd,nd,kd,mpd)
    bk = ww.getWaveletC(nc,kc,ck,nb,kb)

        #Normalize
    nck = normalizeMAAWOS(ck)
    ndk = normalizeMAAWOS(dk)
    ncw = normalizeMAAWOS(cw)
    ndw = normalizeMAAWOS(dw)
    #"""
    #Wavelet interpolation
    error = 0.001
    freq = 0.49
    dt = 0.004
    scale = 4
    nc2 = scale*(nc-1)+1
    #nb2 = scale*(nb-1)+1
    dt2 = dt/scale
    nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
    ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
    ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
    ndw = interpolate(nd,kc,ndw,dt,nc2,dt2,error,freq)
    nc = nc2
    #nb = nb2
    dt = dt2
    #"""
    dti = dt

    incwNoise[i] = ncw
    inckNoise[i] = nck
    indwNoise[i] = ndw
    indkNoise[i] = ndk

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

 #Plotting
  #############1 Plot######################################
  #directory = "./thesisFigures/3pt5/"
  directory = None
  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -3.0,3.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax]],[2.0,2.0,2.0,2.0,2.0,2.0]
  hsize,vsize = 960,350
  title= "f Noises"
  plotting.plotTracesSideBySide(st,[fNoise[0],fNoise[1],fNoise[2],fNoise[3],fNoise[4],fNoise[5]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -9.0,9.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax]],[4.0,4.0,4.0,4.0,4.0,4.0]
  hsize,vsize = 960,350
  title= "g Noises"
  plotting.plotTracesSideBySide(st,[gNoise[0],gNoise[1],gNoise[2],gNoise[3],gNoise[4],gNoise[5]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)


  print "GN rmsNoise"
  dump(rmsNoise)
  #GN Meas
  pngDir = directory
  #pngDir = None
  title = "RMS NoisesGN"
  maxrmsri = 0.5#max(rmsNoise)#0.15
  minrmsrf = 0.0#min(rmsNoise)
  siter = Sampling(nnrmsf,0.2,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Noise-to-signal ratio",None,2.0
  plotting.plotMeasInSamePlotNoise(siter, [rmsNoise],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "WaveletsCNoiseGN"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWaveletsNoise(st,[inckNoise[0],incwNoise[0],incwNoise[1],incwNoise[2],incwNoise[3],incwNoise[4],incwNoise[5]],hint=hint,hsize=hsize,vsize=vsize,title=title,pngDir=pngDir,paper=True,twocol=True)
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "WaveletsDNoiseGN"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWaveletsNoise(st,[indkNoise[0],indwNoise[0],indwNoise[1],indwNoise[2],indwNoise[3],indwNoise[4],indwNoise[5]],hint=hint,hsize=hsize,vsize=vsize,title=title,pngDir=pngDir,paper=True,twocol=True)

def noiseComparisonCyclic():
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsfs = [0.0,0.2,0.4,0.6,0.8,1.0]
  nnrmsf = len(nrmsfs)
  rmsNoise = zerofloat(nnrmsf)
  fNoise = zerofloat(nt,nnrmsf)
  gNoise = zerofloat(nt,nnrmsf)
  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  sfac = 0.0
  ncwNoise = zerofloat(nc,nnrmsf)
  ndwNoise = zerofloat(nc,nnrmsf)
  nckNoise = zerofloat(nc,nnrmsf)
  ndkNoise = zerofloat(nc,nnrmsf)
  scale = 4
  nc2 = scale*(nc-1)+1
  incwNoise = zerofloat(nc2,nnrmsf)
  indwNoise = zerofloat(nc2,nnrmsf)
  inckNoise = zerofloat(nc2,nnrmsf)
  indkNoise = zerofloat(nc2,nnrmsf)
  for i in range(nnrmsf):
    nrmsf = nrmsfs[i]
    nrmsg = nrmsf
    #Create synthetic f and g.
    p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
    freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
    du = computeBackDiff(u)
    fNoise[i] = f
    gNoise[i] = g

    #Wavelet estimation parameters
    nc,kc = 81,-40# sampling for wavelet H 
    nd,kd = nc,kc
    nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
    sfac = 0.0

    #set tmin and tmax 
    tmin = tmin
    tmax = tmax

    #Estimate wavelet
    maxpc = 0.01
    niter = 500
    ww = WaveletWarpingCBCyclic()
    ww.setMinPercentChange(maxpc)#units are percentage.
    ww.setTimeRange(tmin,tmax)

    #First guesses of c and b. 
    bone = zerofloat(nb)
    bone[-kb] = 1.0
    bguess = copy(bone)
    hstabfact = 0.0
    hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
    cguess = copy(hw)
    cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
    lastIter = ww.getLastIter()
    print "lastIter = "+str(lastIter)
    allResRmsAllS = ww.getAllResRmsAllS()
    rmsNoise[i] = allResRmsAllS[lastIter]
    #Estimated Wavelets
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get know wavelets
    ck = getWavelet(freqc,decayc,nc,kc,mpc)
    dk = getWavelet(freqd,decayd,nd,kd,mpd)
    bk = ww.getWaveletC(nc,kc,ck,nb,kb)

        #Normalize
    nck = normalizeMAAWOS(ck)
    ndk = normalizeMAAWOS(dk)
    ncw = normalizeMAAWOS(cw)
    ndw = normalizeMAAWOS(dw)
    #"""
    #Wavelet interpolation
    error = 0.001
    freq = 0.49
    dt = 0.004
    scale = 4
    nc2 = scale*(nc-1)+1
    #nb2 = scale*(nb-1)+1
    dt2 = dt/scale
    nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
    ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
    ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
    ndw = interpolate(nd,kc,ndw,dt,nc2,dt2,error,freq)
    nc = nc2
    #nb = nb2
    dt = dt2
    #"""
    dti = dt

    incwNoise[i] = ncw
    inckNoise[i] = nck
    indwNoise[i] = ndw
    indkNoise[i] = ndk

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

 #Plotting
  #############1 Plot######################################
  directory = "./thesisFigures/3pt5/"
  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -3.0,3.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax]],[2.0,2.0,2.0,2.0,2.0,2.0]
  hsize,vsize = 960,350
  title= "f Noises"
  plotting.plotTracesSideBySide(st,[fNoise[0],fNoise[1],fNoise[2],fNoise[3],fNoise[4],fNoise[5]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  #pngDir = None
  pngDir = directory
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -9.0,9.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude","Amplitude","Amplitude","Amplitude","Amplitude","Amplitude"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax]],[4.0,4.0,4.0,4.0,4.0,4.0]
  hsize,vsize = 960,350
  title= "g Noises"
  plotting.plotTracesSideBySide(st,[gNoise[0],gNoise[1],gNoise[2],gNoise[3],gNoise[4],gNoise[5]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)


  #GN Meas
  rmsNoiseGN = [0.00215383,0.10801400,0.20168807,0.28547618,0.36703083,0.45062494]
  pngDir = directory
  #pngDir = None
  title = "RMS NoisesCyclic"
  maxrmsri = 0.5#max(rmsNoise)#0.15
  minrmsrf = 0.0#min(rmsNoise)
  siter = Sampling(nnrmsf,0.2,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Noise-to-signal ratio",None,0.2
  plotting.plotMeasInSamePlotNoise(siter, [rmsNoiseGN,rmsNoise],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "WaveletsCNoiseCyclic"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWaveletsNoise(st,[inckNoise[0],incwNoise[0],incwNoise[1],incwNoise[2],incwNoise[3],incwNoise[4],incwNoise[5]],hint=hint,hsize=hsize,vsize=vsize,title=title,pngDir=pngDir,paper=True,twocol=True)
  pngDir = directory
  hsize = 960
  vsize = 350
  #pngDir = None     
  title = "WaveletsDNoiseCyclic"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWaveletsNoise(st,[indkNoise[0],indwNoise[0],indwNoise[1],indwNoise[2],indwNoise[3],indwNoise[4],indwNoise[5]],hint=hint,hsize=hsize,vsize=vsize,title=title,pngDir=pngDir,paper=True,twocol=True)

def goSinopecGN():
  directory = None
  #directory = "./thesisFigures/4pt2/"
  nb,kb = 23,-11#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/0sinopec/"
  #directory = "./slides/2sinopec/"
  #directory = "./slides/3sinopec/"
  #directory = "./slides/4sinopec/"
  #directory = "./slides/5sinopec/"
  #directory = "./slides/10sinopec/"
  #get sino images
  x0,nx = 260,100
  f,gNoNR,u = getSinoImage(x0,nx)
  du = computeBackDiff2D(u)
  SimplePlot.asPixels(u)
  SimplePlot.asPixels(du)

  #halfwidth = 0
  halfwidth = 2
  #halfwidth = 3
  #halfwidth = 4
  #halfwidth = 5
  #halfwidth = 10
  if halfwidth==0:
    print "No filtering"
    g = copy(gNoNR)
  else:
    ref = RecursiveExponentialFilter(halfwidth)
    g = zerofloat(len(gNoNR[0]),len(gNoNR))
    ref.apply2(gNoNR,g)

  nx = len(f)
  nt = len(f[0])
  ng = len(g[0])
  nu = len(u[0])

  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  niter = 300
  ww = WaveletWarpingCBGN()
  ww.setMinPercentChange(0.01)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  print "beginning h"
  dump(hw)
  cguess = copy(hw)
  print "cguess";
  dump(cguess);
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
  allResRmsAllS = ww.getAllResRmsAllS()
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
  sg = warp.applyS(u,g)
  bg = ww.applyC(nb,kb,bw,g)


  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  g = ww.makeRms1(225,800,g)
  gNoNR = ww.makeRms1(225,800,gNoNR)
  f = ww.makeRms1(f)
  fZoom = ww.makeRms1(100,200,f)
  sbgZoom = ww.makeRms1(100,200,sbg)
  sgZoom = ww.makeRms1(100,200,sg)
  bg = ww.makeRms1(225,800,bg)
  csbgZoom = ww.makeRms1(100,200,csbg)
  hsgZoom = ww.makeRms1(100,200,hsg)
  csbgmhsg = ww.makeRms1(sub(csbg,hsg))

  #Print
  dt = 0.004
  dx = 0.015
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)
  print "Number of time samples: "+str(nt)

  #Plotting
    #Data
  pngDir = directory
  #pngDir = None
  title= "f small" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  hsize,vsize = 480,560
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "g No Noise reduction SMall" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  hsize,vsize = 480,560
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[gNoNR],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "gNoNR gNR large" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[gNoNR,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "gNR large" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "f csbg csbg hsg GN" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,200*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)","Distance (km)","Distance (km)","Distance (km)"]
  hminmax = None
  hint = 1.0
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[fZoom,csbgZoom,csbgZoom,hsgZoom],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "f csbg csbg sg GN" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,200*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)","Distance (km)","Distance (km)"]
  hminmax = None
  hint = 1.0
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[fZoom,csbgZoom,csbgZoom,sgZoom],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)


  pngDir = directory
  #pngDir = None
  title= "du large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  clipmin,clipmax = 1.3,1.7
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySideDU(st,sx,[du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clipmin=clipmin,clipmax=clipmax,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "duabc"
  dump(du[0]);

  dusmall = copy(399,len(du),101,0,1,1,du)
  SimplePlot.asPixels(dusmall)
  print "maxdusmall = "+str(max(dusmall))
  print "mindusmall = "+str(min(dusmall))

  vpvs = sub(mul(du,2.0),1.0)
  pngDir = directory
  #pngDir = None
  title= "vpvs large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  clipmin,clipmax = 1.6,2.4
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySideDU(st,sx,[vpvs],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clipmin=clipmin,clipmax=clipmax,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)  
  print "duabc"
  dump(du[0]);

    #GN Meas
  print "All RMS "
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals GN"
  maxrmsri = 0.736#max(allResRmsAllS)#0.15
  minrmsrf = 0.725#allResRmsAllS[lastIter+1]#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  hsize = 960
  vsize = 350
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,hsize=hsize,vsize=vsize,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)  
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])

    #Normalize
  ncguess = normalizeM(cguess)
  nbguess = normalizeM(bguess)
  ncw = normalizeM(cw)
  nbw = normalizeM(bw)
  ndw = normalizeM(dw)
  dti = 0.004

  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndw = interpolate(nc,kc,ndw,dt,nc2,dt2,error,freq)
  nbw = interpolate(nb,kb,nbw,dt,nb2,dt2,error,freq)
  ncguess = interpolate(nc,kc,ncguess,dt,nc2,dt2,error,freq)
  nbguess = interpolate(nb,kb,nbguess,dt,nb2,dt2,error,freq)
  nc = nc2
  nb = nb2
  dt = dt2


  #Wavelets
   #Wavelets
  pngDir = directory
  #pngDir = None     
  title = "Shaping Filter (h) GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncguess],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Estimated c GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)  

  pngDir = directory
  #pngDir = None     
  title = "Estimated d GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Interpolated Guess d GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "ImpulseGuess d GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  ndgw[-kc*4] = 1
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Guess (b) GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbguess],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Estimated b GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

def goSinopecCyclic():
  #directory = None
  directory = "./thesisFigures/4pt2/"
  nb,kb = 23,-11#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/0sinopec/"
  #directory = "./slides/2sinopec/"
  #directory = "./slides/3sinopec/"
  #directory = "./slides/4sinopec/"
  #directory = "./slides/5sinopec/"
  #directory = "./slides/10sinopec/"
  #get sino images
  x0,nx = 260,100
  f,gNoNR,u = getSinoImage(x0,nx)
  du = computeBackDiff2D(u)
  SimplePlot.asPixels(u)
  SimplePlot.asPixels(du)

  #halfwidth = 0
  halfwidth = 2
  #halfwidth = 3
  #halfwidth = 4
  #halfwidth = 5
  #halfwidth = 10
  if halfwidth==0:
    print "No filtering"
    g = copy(gNoNR)
  else:
    ref = RecursiveExponentialFilter(halfwidth)
    g = zerofloat(len(gNoNR[0]),len(gNoNR))
    ref.apply2(gNoNR,g)

  nx = len(f)
  nt = len(f[0])
  ng = len(g[0])
  nu = len(u[0])

 #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  niter = 100
  ww = WaveletWarpingCBCyclic()
  ww.setMinPercentChange(0.01)
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(0.000)
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone[-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  print "beginning h"
  dump(hw)

  cguess = copy(hw)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  lastIter = ww.getLastIter()
  print "lastIter = "+str(lastIter)
  allResRmsAllS = ww.getAllResRmsAllS()
  #Estimated Wavelets
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)

  #Processing
  warp = Warper()
  sg = warp.applyS(u,g)
  bg = ww.applyC(nb,kb,bw,g)

  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  g = ww.makeRms1(225,800,g)
  gNoNR = ww.makeRms1(225,800,gNoNR)
  f = ww.makeRms1(f)
  fZoom = ww.makeRms1(100,200,f)
  sbgZoom = ww.makeRms1(100,200,sbg)
  sgZoom = ww.makeRms1(100,200,sg)
  bg = ww.makeRms1(225,800,bg)
  csbgZoom = ww.makeRms1(100,200,csbg)
  hsgZoom = ww.makeRms1(100,200,hsg)
  csbgmhsg = ww.makeRms1(sub(csbg,hsg))

  #Print
  dt = 0.004
  dx = 0.015
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)
  print "Number of time samples: "+str(nt)

  #Plotting
    #Data
  pngDir = directory
  #pngDir = None
  title= "f small" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  hsize,vsize = 480,560
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "g No Noise reduction SMall" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  hsize,vsize = 480,560
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[gNoNR],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "gNoNR gNR large" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[gNoNR,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "gNR large" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "f csbg csbg hsg cyclic" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,200*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)","Distance (km)","Distance (km)","Distance (km)"]
  hminmax = None
  hint = 1.0
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[fZoom,csbgZoom,csbgZoom,hsgZoom],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "f csbg csbg sg cyclic" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,200*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)","Distance (km)","Distance (km)","Distance (km)"]
  hminmax = None
  hint = 1.0
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[fZoom,csbgZoom,csbgZoom,sgZoom],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)


  pngDir = directory
  #pngDir = None
  title= "du large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  clipmin,clipmax = 1.3,1.7
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySideDU(st,sx,[du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clipmin=clipmin,clipmax=clipmax,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "duabc"
  dump(du[0]);

  dusmall = copy(399,len(du),101,0,1,1,du)
  SimplePlot.asPixels(dusmall)
  print "maxdusmall = "+str(max(dusmall))
  print "mindusmall = "+str(min(dusmall))

  vpvs = sub(mul(du,2.0),1.0)
  pngDir = directory
  #pngDir = None
  title= "vpvs large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  clipmin,clipmax = 1.6,2.4
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySideDU(st,sx,[vpvs],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clipmin=clipmin,clipmax=clipmax,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)  
  print "duabc"
  dump(du[0]);

    #GN Meas
  print "All RMS cyclic"
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals cyclic"
  maxrmsri = 0.736#max(allResRmsAllS)#0.15
  minrmsrf = 0.725#allResRmsAllS[lastIter+1]#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  hsize = 960
  vsize = 350
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,hsize=hsize,vsize=vsize,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)  
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])

    #Normalize
  ncguess = normalizeM(cguess)
  nbguess = normalizeM(bguess)
  ncw = normalizeM(cw)
  nbw = normalizeM(bw)
  ndw = normalizeM(dw)
  dti = 0.004

  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndw = interpolate(nc,kc,ndw,dt,nc2,dt2,error,freq)
  nbw = interpolate(nb,kb,nbw,dt,nb2,dt2,error,freq)
  ncguess = interpolate(nc,kc,ncguess,dt,nc2,dt2,error,freq)
  nbguess = interpolate(nb,kb,nbguess,dt,nb2,dt2,error,freq)
  nc = nc2
  nb = nb2
  dt = dt2


  #Wavelets
   #Wavelets
  pngDir = directory
  #pngDir = None     
  title = "Shaping Filter (h)cyclic"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncguess],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Estimated c cyclic"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)  

  pngDir = directory
  #pngDir = None     
  title = "Estimated d cyclic"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Interpolated Guess d cyclic"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "ImpulseGuess d cyclic"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  ndgw[-kc*4] = 1
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Guess (b) cyclic"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbguess],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Estimated b cyclic"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

def goGBCGN():
  directory = None
  #directory = "./thesisFigures/4pt2/"
  nb,kb = 23,-11#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/0sinopec/"
  #directory = "./slides/2sinopec/"
  #directory = "./slides/3sinopec/"
  #directory = "./slides/4sinopec/"
  #directory = "./slides/5sinopec/"
  #directory = "./slides/10sinopec/"
  #get sino images
  x0,nx = 260,100
  f,gNoNR,u = getSinoImage(x0,nx)
  du = computeBackDiff2D(u)
  SimplePlot.asPixels(u)
  SimplePlot.asPixels(du)

  #halfwidth = 0
  halfwidth = 2
  #halfwidth = 3
  #halfwidth = 4
  #halfwidth = 5
  #halfwidth = 10
  if halfwidth==0:
    print "No filtering"
    g = copy(gNoNR)
  else:
    ref = RecursiveExponentialFilter(halfwidth)
    g = zerofloat(len(gNoNR[0]),len(gNoNR))
    ref.apply2(gNoNR,g)

  nx = len(f)
  nt = len(f[0])
  ng = len(g[0])
  nu = len(u[0])

  #Wavelet estimation parameters
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  niter = 300
  ww = WaveletWarpingCBGN()
  ww.setMinPercentChange(0.01)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  #First guesses of c and b. 
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bguess = copy(bone)
  hstabfact = 0.0
  hw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
  print "beginning h"
  dump(hw)
  cguess = copy(hw)
  print "cguess";
  dump(cguess);
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
  allResRmsAllS = ww.getAllResRmsAllS()
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
  sg = warp.applyS(u,g)
  bg = ww.applyC(nb,kb,bw,g)


  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  g = ww.makeRms1(225,800,g)
  gNoNR = ww.makeRms1(225,800,gNoNR)
  f = ww.makeRms1(f)
  fZoom = ww.makeRms1(100,200,f)
  sbgZoom = ww.makeRms1(100,200,sbg)
  sgZoom = ww.makeRms1(100,200,sg)
  bg = ww.makeRms1(225,800,bg)
  csbgZoom = ww.makeRms1(100,200,csbg)
  hsgZoom = ww.makeRms1(100,200,hsg)
  csbgmhsg = ww.makeRms1(sub(csbg,hsg))

  #Print
  dt = 0.004
  dx = 0.015
  print "tmin: "+str(tmin)+" or "+str(tmin*dt)
  print "tmax: "+str(tmax)+" or "+str(tmax*dt)
  print "Number of time samples: "+str(nt)

  #Plotting
    #Data
  pngDir = directory
  #pngDir = None
  title= "f small" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  hsize,vsize = 480,560
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "g No Noise reduction SMall" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  hsize,vsize = 480,560
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[gNoNR],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  hsize=hsize,vsize=vsize,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "gNoNR gNR large" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)"]
  hminmax = None
  hint = 0.5
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[gNoNR,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "gNR large" 
  print title
  st = Sampling(len(g[0]),dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 225*dt,800*dt
  vlabel,vminmax,vint = "PS Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  pngDir = directory
  #pngDir = None
  title= "f csbg csbg hsg GN" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,200*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)","Distance (km)","Distance (km)","Distance (km)"]
  hminmax = None
  hint = 1.0
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[fZoom,csbgZoom,csbgZoom,hsgZoom],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = directory
  #pngDir = None
  title= "f csbg csbg sg GN" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,200*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)","Distance (km)","Distance (km)"]
  hminmax = None
  hint = 1.0
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[fZoom,csbgZoom,csbgZoom,sgZoom],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)


  pngDir = directory
  #pngDir = None
  title= "du large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  clipmin,clipmax = 1.3,1.7
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySideDU(st,sx,[du],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clipmin=clipmin,clipmax=clipmax,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "duabc"
  dump(du[0]);

  dusmall = copy(399,len(du),101,0,1,1,du)
  SimplePlot.asPixels(dusmall)
  print "maxdusmall = "+str(max(dusmall))
  print "mindusmall = "+str(min(dusmall))

  vpvs = sub(mul(du,2.0),1.0)
  pngDir = directory
  #pngDir = None
  title= "vpvs large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  clipmin,clipmax = 1.6,2.4
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySideDU(st,sx,[vpvs],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clipmin=clipmin,clipmax=clipmax,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)  
  print "duabc"
  dump(du[0]);

    #GN Meas
  print "All RMS "
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals GN"
  maxrmsri = 0.736#max(allResRmsAllS)#0.15
  minrmsrf = 0.725#allResRmsAllS[lastIter+1]#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  hsize = 960
  vsize = 350
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,hsize=hsize,vsize=vsize,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  paper=True,twocol=True)  
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])

    #Normalize
  ncguess = normalizeM(cguess)
  nbguess = normalizeM(bguess)
  ncw = normalizeM(cw)
  nbw = normalizeM(bw)
  ndw = normalizeM(dw)
  dti = 0.004

  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  nb2 = scale*(nb-1)+1
  dt2 = dt/scale
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndw = interpolate(nc,kc,ndw,dt,nc2,dt2,error,freq)
  nbw = interpolate(nb,kb,nbw,dt,nb2,dt2,error,freq)
  ncguess = interpolate(nc,kc,ncguess,dt,nc2,dt2,error,freq)
  nbguess = interpolate(nb,kb,nbguess,dt,nb2,dt2,error,freq)
  nc = nc2
  nb = nb2
  dt = dt2


  #Wavelets
   #Wavelets
  pngDir = directory
  #pngDir = None     
  title = "Shaping Filter (h) GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncguess],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Estimated c GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)  

  pngDir = directory
  #pngDir = None     
  title = "Estimated d GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Interpolated Guess d GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "ImpulseGuess d GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  ndgw[-kc*4] = 1
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Guess (b) GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbguess],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)

  pngDir = directory
  #pngDir = None     
  title = "Estimated b GN"
  hint = None
  hsize = 960
  vsize = 350
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbw],dotsonsticks=True,hsize=hsize,vsize=vsize,hint=hint,title=title,pngDir=pngDir,paper=True,twocol=True)



def getWavelet(fpeak,decay,nc,kc,mp=False):
  x = zerofloat(nc)
  x[-kc] = 1.0
  return synthetic.addWavelet(fpeak,decay,x,mp)

def normalizeM(f):
  nf = len(f)
  sums = 0.0

  y = zerofloat(nf)
  maxv = max(f)
  f = div(f,maxv)
  return f

def normalizeMAAWOS(f):
  minf = min(f)
  maxf = max(f)
  absminf = abs(minf) 
  absmaxf = abs(maxf) 
  if absminf<absmaxf:
    return mul(1.0/maxf,f)
  else:
    return mul(1.0/minf,f)


def interpolate(nc,kc,c,dt,nc2,dt2,error,freq):
  ci = zerofloat(nc2)
  si = SincInterpolator.fromErrorAndFrequency(error,freq)
  si.interpolate(nc,dt,-kc*dt,c,nc2,dt2,-kc*dt,ci)
  return ci

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
  gain(100,fr)
  gain(100,gr)
  return fr,gr,ur
  #gain(100,f)
  #gain(100,g)
  #return f,g,u

def getGBCImage(x0,nx):
  dataDir = "/Users/Chris/data/gbc/dat/"
  n1f,n1g,d1,f1 = 2000,852,0.002,0.0
  n2,d2,f2 =  721,0.0150,0.000
  f = readImage(dataDir+"pp.dat",n1f,n2)
  g = readImage(dataDir+"ps1.dat",n1g,n2)
  u = readImage(dataDir+"u_sag_linear.dat",n1f,n2)
  u = add(u,rampfloat(0.0,1.0,0.0,n1f,n2))
  fr = zerofloat(n1f,nx)
  gr = zerofloat(n1g,nx)
  ur = zerofloat(n1f,nx)
  for ix in range(0,nx):
    fr[ix] = f[x0+ix]
    gr[ix] = g[x0+ix]
    ur[ix] = u[x0+ix]
  gain(100,fr)
  gain(100,gr)
  return fr,gr,ur
  #gain(100,f)
  #gain(100,g)
  #return f,g,u


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














  

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

