#############################################################################
# Demo of 2 wavelet estimations from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.awt.ColorMap import *
from edu.mines.jtk.lapack import *
from wwarp import WaveletWarpingCBGN, WaveletWarpingCBCyclic, Warper
from wwarp import WaveletWarping#display only, do not use
from wwarp import ShapingFilter 
import synthetic
import plotting
from java.util import Random

############################################################################


def main(args):
  #displaySimpleCaseOfWarpingWithWavelets()
  #goWarpDifference2t()
  #goWarpDifferencelog()
  #goSimpleEstimate()
  #oneDNoNoiseHighSqueezingTwoWavelets()
  #increaseBy1NoNoise()
  #increaseBy1Sinopec()
  #increaseBy1UseOldSolutionToStartSinopec()
  #oneDNoiseHighSqueezingTwoWavelets()
  #oneDNoiseHighSqueezingTwoWaveletsCBScaling()
  #goSinopecGN()
  #goSinopecCyclic()
  #goDiffEstimate()
  #goDiffNoiseEstimate()
  #go2D()
  #chooseNoise()
  #goSinopecGNAlphaAnalysis()


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



#A separate test to figure out if the Gauss-Newton method is working.
def oneDNoNoiseHighSqueezingTwoWavelets():
  directory = None#"./slides/nonoiseHighSqueezing/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
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
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 
  na,ka = nb,kb

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

  pngDir = directory
  #pngDir = None
  title= "[f,csbg] [f,hsg]"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],None
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
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
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "[f,sg] [f,hsg]"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],None
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-10.0
  hmax1,hmax2 = 1.0,10.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,sg],[f,hsg]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing=tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  title= "[f,f] [f,f]"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],None
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  color = [Color.BLACK,Color.BLACK]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,f],[f,f]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  title= "[f,CSBg] [f,f]"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],None
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,csbg],[f,f]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None
  title= "f g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-5.0
  hmax1,hmax2 = 1.0,5.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  #pngDir = directory
  pngDir = None
  title= "f HSg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-5.0
  hmax1,hmax2 = 1.0,5.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  #pngDir = directory
  pngDir = None
  title= "f HSgZoom"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[0.0,0.928],0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-3.0
  hmax1,hmax2 = 1.0,3.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  #pngDir = directory
  pngDir = None
  title= "f HSgZoom"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[0.0,0.928],0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-3.0
  hmax1,hmax2 = 1.0,3.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  #pngDir = directory
  pngDir = None
  title= "f CSBgZoom"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[0.0,0.928],0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-3.0
  hmax1,hmax2 = 1.0,3.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  #pngDir = directory
  pngDir = None
  title= "f csbg p sbg bg q g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-8.0
  hmax1,hmax2 = 1.0,8.0
  hlabel = ["f","HSg","CSBg","du","p","SBg","Bg","q","g"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin1,hmax1],\
  [hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin2,hmax2]]
  hint = [None,None,None,None,None,None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,csbg,du,p,sbg,bg,q,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

    #GN Meas
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals"
  maxrmsri = max(allResRmsAllS)#0.15
  minrmsrf = 0.0#min(allResRmsAllS)#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter+1],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])

  pngDir = directory
  #pngDir = None
  title = "All Rms ResidualsZoom"
  maxrmsri = 0.006#max(allResRmsAllS)#0.15
  minrmsrf = 0.0#min(allResRmsAllS)#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[177.0,lastIter+1],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])


  #pngDir = directory
  pngDir = None
  title = "Total Initial RMS of Residuals (Black) Final RMS of Residuals (Blue)"
  for i in range(lastIter,niter):
    allResRmsInitS[i] = allResRmsInitS[lastIter]
    allResRmsFinaS[i] = allResRmsFinaS[lastIter]
  maxrmsri = max(allResRmsInitS)#0.15
  minrmsrf = min(allResRmsFinaS)#0.0
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
  
  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
  ndw = interpolate(nc,kc,ndw,dt,nc2,dt2,error,freq)
  nbk = interpolate(nb,kb,nbk,dt,nb2,dt2,error,freq)
  nbw = interpolate(nb,kb,nbw,dt,nb2,dt2,error,freq)
  ncguess = interpolate(nc,kc,ncguess,dt,nc2,dt2,error,freq)
  nbguess = interpolate(nb,kb,nbguess,dt,nb2,dt2,error,freq)
  nc = nc2
  nb = nb2
  dt = dt2
  #"""

    #Wavelets
  #pngDir = directory
  pngDir = None     
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated c"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  title = "Known c"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[nck],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None     
  title = "Estimated d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  title = "known d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndk],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  #pngDir = directory
  pngDir = None     
  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated b"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  title = "Known b"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbk],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


#A separate test to figure out if the Gauss-Newton method is working.
def oneDNoiseHighSqueezingTwoWavelets():
  directory = "./slides/noiseHighSqueezing/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.40
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

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
  ww.setMinPercentChange(maxpc)#units are percentage.
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

  pngDir = directory
  #pngDir = None
  title= "[f,csbg] [f,hsg]"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],None
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
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
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  title= "[f,f] [f,f]"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],None
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  color = [Color.BLACK,Color.BLACK]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,f],[f,f]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  title= "[f,CSBg] [f,f]"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,0.98
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],None
  hint1,hint2 = 0.5,1.0
  hmin1,hmin2 = -1.0,-2.5
  hmax1,hmax2 = 1.0,2.5
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  color = [Color.BLACK,Color.RED]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plot2TracesInSamePlotSideBySideWithOtherPlots(st,\
  [[f,csbg],[f,f]],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  color=color,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None
  title= "f g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-5.0
  hmax1,hmax2 = 1.0,5.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  #pngDir = directory
  pngDir = None
  title= "f HSg"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-5.0
  hmax1,hmax2 = 1.0,5.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  #pngDir = directory
  pngDir = None
  title= "f HSgZoom"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[0.0,0.928],0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-3.0
  hmax1,hmax2 = 1.0,3.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  #pngDir = directory
  pngDir = None
  title= "f CSBgZoom"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[0.0,0.928],0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-3.0
  hmax1,hmax2 = 1.0,3.0
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[f,csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  #pngDir = directory
  pngDir = None
  title= "f csbg p sbg bg q g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-8.0
  hmax1,hmax2 = 1.0,8.0
  hlabel = ["f","HSg","CSBg","du","p","SBg","Bg","q","g"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin1,hmax1],\
  [hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin2,hmax2]]
  hint = [None,None,None,None,None,None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,csbg,du,p,sbg,bg,q,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

    #GN Meas
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals"
  maxrmsri = max(allResRmsAllS)#0.15
  minrmsrf = allResRmsAllS[lastIter+1]#min(allResRmsAllS)#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter+1],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])

  pngDir = directory
  #pngDir = None
  title = "All Rms ResidualsZoom"
  maxrmsri = 0.006#max(allResRmsAllS)#0.15
  minrmsrf = 0.0#min(allResRmsAllS)#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter+1],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])


  #pngDir = directory
  pngDir = None
  title = "Total Initial RMS of Residuals (Black) Final RMS of Residuals (Blue)"
  for i in range(lastIter,niter):
    allResRmsInitS[i] = allResRmsInitS[lastIter]
    allResRmsFinaS[i] = allResRmsFinaS[lastIter]
  maxrmsri = max(allResRmsInitS)#0.15
  minrmsrf = min(allResRmsFinaS)#0.0
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
  
  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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
  dti = 0.004

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
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  ndk = interpolate(nc,kc,ndk,dt,nc2,dt2,error,freq)
  ndw = interpolate(nc,kc,ndw,dt,nc2,dt2,error,freq)
  nbk = interpolate(nb,kb,nbk,dt,nb2,dt2,error,freq)
  nbw = interpolate(nb,kb,nbw,dt,nb2,dt2,error,freq)
  ncguess = interpolate(nc,kc,ncguess,dt,nc2,dt2,error,freq)
  nbguess = interpolate(nb,kb,nbguess,dt,nb2,dt2,error,freq)
  nc = nc2
  nb = nb2
  dt = dt2
  #"""

    #Wavelets
  #pngDir = directory
  pngDir = None     
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated c"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  title = "Known c"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[nck],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None     
  title = "Estimated d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  title = "known d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndk],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  #pngDir = directory
  pngDir = None     
  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated b"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  title = "Known b"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbk],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

#A separate test to figure out if the Gauss-Newton method is working.
def oneDNoiseHighSqueezingTwoWaveletsCBScaling():
  directory = "./slides/scaling/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.20
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

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
  ww.setMinPercentChange(maxpc)#units are percentage.
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
  cwtenth = mul((0.1),cw)
  bw = cbw[1]
  bwten = mul(10,bw)

  dw = ww.getWaveletC(nb,kb,bw,nc,kc)

  #Processing
  warp = Warper()
  bg = ww.applyC(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  bgScaled = ww.applyC(nb,kb,bwten,g)
  sbgScaled = warp.applyS(u,bgScaled)
  csbgScaled = ww.applyC(nc,kc,cwtenth,sbgScaled)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))
  SimplePlot.asPoints(warp.applyS(u,g))

  dt = 0.004

  pngDir = directory
  #pngDir = None
  title= "residual scaled"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[0.0,0.928],0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-0.5
  hmax1,hmax2 = 1.0,0.5
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[sub(csbgScaled,f)],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "residual not scaled"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",[0.0,0.928],0.5
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-0.5
  hmax1,hmax2 = 1.0,0.5
  hlabel = ["Amplitude","Amplitude"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2]]
  hint = [None,None]
  hsize,vsize = 960,560
  tilespacing = None
  plotting.plotTracesSideBySide(st,[sub(csbg,f)],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)



def increaseBy1NoNoise():
  directory = "./slides/Increase1Synthetic/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.16,0.07
  freqd,decayd = 0.08,0.07
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0#0.4
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

  nhfinal = 85
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

    #Get iteration information
    lastIter = ww.getLastIter()
    allResRmsAllS = ww.getAllResRmsAllS()
    allRmsResCB[i] = allResRmsAllS[lastIter]
    warp = Warper()
    sg = warp.applyS(u,g)
    hsg = ww.applyC(nc,kc,hw,sg)
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
  title = "Increaseby1"+" nhfinal "+str(nhfinal)+" nhinitial "+str(nhinitial)+" maxiter "+str(niter)+" noise "+str(nrmsf)+" r0 "+str(r0)+" r1 "+str(r1)+" nhfinal"+str(nhfinal)+"Syn Increase nb and nh nc constant first rms ri (black) last rms rf (blue)"
  maxrms = 0.14#max([max(allRmsResCB),max(allRmsResH)])
  minrms = 0.0#min([min(allRmsResCB),min(allRmsResH)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(allRmsResH))
  print "rmsrcb = "+str(len(allRmsResCB))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of residuals",[minrms,maxrms],0.04
  hlabel,hminmax,hint = "nh=nb+80",[nhinitial,nhfinal],1.0
  plotting.plotMeasInSamePlotLargeMarks(si, [allRmsResH,allRmsResCB],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

def increaseBy1Sinopec():
  directory = "./slides/Increase1Sino/"
  x0,nx = 260,100
  f,gNoNR,u = getSinoImage(x0,nx)

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


  maxpc = 0.01
  niter = 300
  #set tmin and tmax 
  tmin = 100
  tmax = 500

  nhfinal = 91
  nhinitial = 81
  khinitial = 0
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
    ww.setLineSearchMinScale(0.0000)
  
    #First guesses of c and b. 
    bone = zerofloat(nb)
    bone [-kb] = 1.0
    bguess = copy(bone)
    hstabfact = 0.0
    hw =  ww.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)
    if (i != 0):
      cguess =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
      cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
      #Estimated Wavelets
      cw = cbw[0]
      bw = cbw[1]
      dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    #Get iteration information
    if (i != 0):
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
    nb = nb+1
    nh = nh+1

  allRmsResCB[0] = allRmsResH[0]
  pngDir = None
  pngDir = directory
  title = "Increaseby1SinohalfWidth"+str(halfwidth)
  maxrms = 0.770#max([max(allRmsResCB),max(allRmsResH)])
  minrms = 0.761#min([min(allRmsResCB),min(allRmsResH)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(allRmsResH))
  print "rmsrcb = "+str(len(allRmsResCB))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+10",[nhinitial,nhfinal],1.0
  plotting.plotMeasInSamePlotLargeMarks(si, [allRmsResH,allRmsResCB],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

def increaseBy1UseOldSolutionToStartSinopec():
  directory = "./slides/Increase1SinoUseOldSolutionToStart/"
  x0,nx = 260,100
  f,gNoNR,u = getSinoImage(x0,nx)

  #halfwidth = 0
  #halfwidth = 2
  halfwidth = 2
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


  maxpc = 0.01
  niter = 300
  #set tmin and tmax 
  tmin = 100
  tmax = 500

  nhfinal = 91
  nhinitial = 81
  khinitial = 0
  nb = 1
  kb = 0
  nc = nhinitial
  kc = khinitial
  nh = nhinitial
  kh = khinitial
  n = (nhfinal-nh)+1
  allRmsResH = zerofloat(n)
  allRmsResCB = zerofloat(n)
  bone = zerofloat(nb)
  bone [-kb] = 1.0
  bw = copy(bone)
  hstabfact = 0.0
  ww = WaveletWarpingCBGN()
  ww.setMinPercentChange(maxpc)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setLineSearchMinScale(0.0000)
  cw =  ww.getWaveletC(nc,kc,nb,kb,bone,hstabfact,u,f,g)
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
    ww.setLineSearchMinScale(0.0000)


    #First guesses of c and b. 
    #print "old qrl b"
    #dump(bw)
    #print "old crl c"
    #dump(cw)
    bw = addZeros(bw,nb)
    cw = addZeros(cw,nc)
    print "old wth zero qrl b"
    dump(bw)
    bone = zerofloat(nb)
    bone [-kb] = 1.0
    hw =  ww.getWaveletC(nh,kh,nb,kb,bone,hstabfact,u,f,g)
    if (i != 0 ):
      cbw = ww.getWaveletCInverseB(nb,kb,bw,nc,kc,cw,u,f,g,niter)
      #Estimated Wavelets
      cw = cbw[0]
      bw = cbw[1]
      print "new qrl b"
      dump(bw)
      #print "new crl c"
      #dump(cw)
      dw = ww.getWaveletC(nb,kb,bw,nc,kc)

    if (i != 0 ):
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
    nb = nb+1
    nh = nh+1

  allRmsResCB[0] = allRmsResH[0]
  pngDir = None
  pngDir = directory
  title = "Increaseby1SinohalfWidth"+str(halfwidth)
  maxrms = 0.770#max([max(allRmsResCB),max(allRmsResH)])
  minrms = 0.761#min([min(allRmsResCB),min(allRmsResH)])
  print "min = "+str(minrms)
  print "max = "+str(maxrms)
  print "n = "+str(n)
  print "rmsrh = "+str(len(allRmsResH))
  print "rmsrcb = "+str(len(allRmsResCB))
  si = Sampling(n,1.0,nc)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrms,maxrms],None
  hlabel,hminmax,hint = "nh=nb+10",[nhinitial,nhfinal],1.0
  plotting.plotMeasInSamePlotLargeMarks(si, [allRmsResH,allRmsResCB],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)




def goSinopecGN():
  directory = None
  #directory = "./slides/pc0_021_10_0sinopec100/"
  #directory = "./slides/21_10_0sinopec100/"
  #directory = "./slides/5_2_3sinopec100/"
  #nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/7_3_3sinopec100/"
  #nb,kb = 7,-3#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/9_4_3sinopec100/"
  #nb,kb = 9,-4#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/21_10_3sinopec100/"
  #nb,kb = 21,-10#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/31_15_3sinopec100/"
  #nb,kb = 31,-15#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/23_11_3sinopec100/"
  nb,kb = 23,-11#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/19_9_3sinopec100/"
  #nb,kb = 19,-9#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/17_8_3sinopec100/"
  #nb,kb = 17,-8#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/15_7_3sinopec100/"
  #nb,kb = 15,-7#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/13_6_3sinopec100/"
  #nb,kb = 13,-6#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/11_5_3sinopec100/"
  #nb,kb = 11,-5#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/25_12_3sinopec100/"
  #nb,kb = 25,-12#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/27_13_3sinopec100/"
  #nb,kb = 27,-13#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/29_14_3sinopec100/"
  #nb,kb = 29,-14#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/33_16_3sinopec100/"
  #nb,kb = 33,-16#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/35_17_3sinopec100/"
  #nb,kb = 35,-17#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/37_18_3sinopec100/"
  #nb,kb = 37,-18#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/39_19_3sinopec100/"
  #nb,kb = 39,-19#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/41_20_3sinopec100/"
  #nb,kb = 41,-20#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/43_21_3sinopec100/"
  #nb,kb = 43,-21#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/45_22_3sinopec100/"
  #nb,kb = 45,-22#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/47_23_3sinopec100/"
  #nb,kb = 47,-23#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/49_24_3sinopec100/"
  #nb,kb = 49,-24#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/51_25_3sinopec100/"
  #nb,kb = 51,-25#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/53_26_3sinopec100/"
  #nb,kb = 53,-26#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/61_30_3sinopec100/"
  #nb,kb = 61,-30#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/71_35_3sinopec100/"
  #nb,kb = 71,-35#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/81_40_3sinopec100/"
  #nb,kb = 81,-40#sampling for inverse wavelet A #Note ka<=kc 
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

  #Display only
  wwOld = WaveletWarping()
  wwOld.setTimeRange(tmin,tmax)
  lbg = wwOld.applyL(u,bg)
  ######

  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  g = ww.makeRms1(225,800,g)
  gNoNR = ww.makeRms1(225,800,gNoNR)
  f = ww.makeRms1(f)
  sbg = ww.makeRms1(sbg)
  sg = ww.makeRms1(sg)
  bg = ww.makeRms1(225,800,bg)
  lbg = ww.makeRms1(225,800,lbg)
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
  slide=True,fracWidth=0.45,fracHeight=0.8,aspectRatio=16.0/9.0)

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
  slide=True,fracWidth=0.45,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "gNoNR large" 
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
  plotting.plotImagesSideBySide(st,sx,[gNoNR],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

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
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "Sg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[sg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "Bg large" 
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
  plotting.plotImagesSideBySide(st,sx,[bg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "LBg large" 
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
  plotting.plotImagesSideBySide(st,sx,[lbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None
  title= "SBg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[sbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "CSBg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "CSBg zoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,150*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "CSBg bottomzoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 450*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)



  pngDir = directory
  #pngDir = None
  title= "HSg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "HSg zoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,150*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "HSg bottom zoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 450*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)




  pngDir = directory
  #pngDir = None
  title= "f large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "f zoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,150*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)



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
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "duabc"
  dump(du[0]);






  pngDir = directory
  #pngDir = None
  title= "f HSg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f,hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "f CSBg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f,csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)




    #GN Meas
  print "All RMS"
  dump(allResRmsAllS)
  print "alphab"
  dump(alphaBS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals"
  maxrmsri = 0.736#max(allResRmsAllS)#0.15
  minrmsrf = 0.725#allResRmsAllS[lastIter+1]#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter+1],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])

  pngDir = directory
  #pngDir = None
  title = "All Rms ResidualsZoom"
  maxrmsri = 0.006#max(allResRmsAllS)#0.15
  minrmsrf = 0.0#min(allResRmsAllS)#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter+1],None
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "RMS = "+str(allResRmsAllS[lastIter])
  print "RMS = "+str(allResRmsAllS[lastIter+1])
  print "RMS = "+str(allResRmsAllS[lastIter+2])


  #pngDir = directory
  pngDir = None
  title = "Total Initial RMS of Residuals (Black) Final RMS of Residuals (Blue)"
  for i in range(lastIter,niter):
    allResRmsInitS[i] = allResRmsInitS[lastIter]
    allResRmsFinaS[i] = allResRmsFinaS[lastIter]
  maxrmsri = max(allResRmsInitS)#0.15
  minrmsrf = min(allResRmsFinaS)#0.0
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
  
  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
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

  #pngDir = directory
  pngDir = None
  title = "Initial 2NormSq of individual and all residuals"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,
  [allRes2NormSqInitS,dataRes2NormSqInitS,bPenaltyRes2NormSqInitS,cPenaltyRes2NormSqInitS,alphaBS,alphaCS],\
  color=None,\
  vlabel=["allRes","dataRes","bPenRes","cPenRes","alphaB","alphaC"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0.0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  pngDir = directory
  #pngDir = None
  title = "alphab"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,
  [bPenaltyRes2NormSqFinaS,cPenaltyRes2NormSqFinaS,alphaBS,alphaCS],\
  color=None,\
  vlabel=["bPenRes","cPenRes","alphaB","alphaC"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0.0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  #pngDir = directory
  pngDir = None
  title = "Step Length and Two Norm of Gradient"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,[stepLengthS,gradient2NormS],\
  color=None,\
  vlabel=["StepLength","2NormGrad"],\
  vminmax=[None,None],\
  vint=[None,None],\
  hlabel="Iterations",hminmax=[0.0,lastIter],hint=None,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)

  pngDir = directory
  #pngDir = None
  title = "2NormDeltaB,2NormDeltaC,2NormB,2NormC,ConditionNumber,RMS Percent Change"
  siter = Sampling(niter,1.0,0.0)
  plotting.plotMeasOnTopOfEachOther(siter,[deltaB2NormS,deltaC2NormS,\
  b2NormS,c2NormS,conditionNumberS,rmsPercentChange],\
  color=None,\
  vlabel=["2NormDelB","2NormDelC","2NormB","2NormC","CondNum","RMSPercChange"],\
  vminmax=[None,None,None,None,None,None],\
  vint=[None,None,None,None,None,None],\
  hlabel="Iterations",hminmax=[0.0,lastIter],hint=None,\
  hsize=960,vsize=760,\
  title=title,pngDir=pngDir,\
  slide=None,fracWidth=None,fracHeight=None,\
  paper=True,onecol=True,twocol=False)


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
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated c"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Interpolated Guess d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "ImpulseGuess d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  ndgw[-kc*4] = 1
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None     
  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated b"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  title = "estimated PP and PS wavelets"
  pngDir = directory
  sh = Sampling(nc,dt,kc*dti)
  plotWaveletsPpPs(sh,ncw,ndw,0.9,0.8,16.0/9.0,title=title,pngDir=pngDir,twocol=True)

  title = "starting PP and PS wavelets"
  pngDir = directory
  plotWaveletsPpPs(sh,ncguess,ndgw,0.9,0.8,16.0/9.0,title=title,pngDir=pngDir,twocol=True)



def goSinopecCyclic():
  #directory = "./slides/cyclic23_11_3sinopec100/"
  directory = None
  nb,kb = 23,-11#sampling for inverse wavelet A #Note ka<=kc 

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

  #Display only
  wwOld = WaveletWarping()
  wwOld.setTimeRange(tmin,tmax)
  lbg = wwOld.applyL(u,bg)
  ######

  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)
  hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

  g = ww.makeRms1(225,800,g)
  gNoNR = ww.makeRms1(225,800,gNoNR)
  f = ww.makeRms1(f)
  sbg = ww.makeRms1(sbg)
  sg = ww.makeRms1(sg)
  bg = ww.makeRms1(225,800,bg)
  lbg = ww.makeRms1(225,800,lbg)
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
  slide=True,fracWidth=0.45,fracHeight=0.8,aspectRatio=16.0/9.0)

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
  slide=True,fracWidth=0.45,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "gNoNR large" 
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
  plotting.plotImagesSideBySide(st,sx,[gNoNR],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

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
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "Sg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[sg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "Bg large" 
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
  plotting.plotImagesSideBySide(st,sx,[bg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "LBg large" 
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
  plotting.plotImagesSideBySide(st,sx,[lbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None
  title= "SBg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[sbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "CSBg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "CSBg zoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,150*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "CSBg bottomzoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 450*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)



  pngDir = directory
  #pngDir = None
  title= "HSg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "HSg zoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,150*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "HSg bottom zoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 450*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)




  pngDir = directory
  #pngDir = None
  title= "f large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "f zoom large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,150*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.1
  hlabel = ["Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)



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
  slide=True,fracWidth=0.7,fracHeight=0.8,aspectRatio=16.0/9.0)
  print "duabc"
  dump(du[0]);






  pngDir = directory
  #pngDir = None
  title= "f HSg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f,hsg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None
  title= "f CSBg large" 
  print title
  st = Sampling(nt,dt,0.0)
  sx = Sampling(nx,dx,x0*dx)
  vmin,vmax = 100*dt,500*dt
  vlabel,vminmax,vint = "PP Time (s)",[vmin,vmax],0.5
  hlabel = ["Distance (km)","Distance (km)"]
  hminmax = None
  hint = 0.25
  tilespacing = None
  clip = 3.0
  plotting.plotImagesSideBySide(st,sx,[f,csbg],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  clip = clip,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)




    #GN Meas
  print "All RMS"
  dump(allResRmsAllS)
  pngDir = directory
  #pngDir = None
  title = "All Rms Residuals"
  maxrmsri = 0.736#max(allResRmsAllS)#0.15
  minrmsrf = 0.725#allResRmsAllS[lastIter]#0.0
  siter = Sampling(niter,1.0,0.0)
  color=[Color.BLACK,Color.RED]
  vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
  hlabel,hminmax,hint = "Iterations",[0.0,lastIter],5.0
  plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
  color=color,\
  vlabel=vlabel, vminmax=vminmax, vint=vint,\
  hlabel=hlabel, hminmax=hminmax, hint=hint,\
  title=title, pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


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
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated c"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Interpolated Guess d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "ImpulseGuess d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  ndgw[-kc*4] = 1
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None     
  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated b"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  title = "estimated PP and PS wavelets"
  pngDir = directory
  sh = Sampling(nc,dt,kc*dti)
  plotWaveletsPpPs(sh,ncw,ndw,0.9,0.8,16.0/9.0,title=title,pngDir=pngDir,twocol=True)

  title = "starting PP and PS wavelets"
  pngDir = directory
  plotWaveletsPpPs(sh,ncguess,ndgw,0.9,0.8,16.0/9.0,title=title,pngDir=pngDir,twocol=True)



def goWarpDifference2t():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 181,-90 # sampling for wavelet H 
  niter = 1
  r0,r1 = 2.0,2.0#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 1.00

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
  pngDir = "./report15/synthetics/"
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -1.5,1.5
  vlabel,vminmax,vint = "time (s)",[vmin,vmax],1.0
  hlabel,hminmax,hint = ["p","q","sq","wlq"],[[hmin,hmax],[hmin,hmax],[hmin,hmax],[hmin,hmax]],1.0
  hsize,vsize = 960,560
  title= "test2tpqwlqsq"
  plotting.plotTracesSideBySide(st,[p,q,sq,wlq],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=["p","q","sq","wlq"],hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)
  ##########################################################
  #plot4TracesSideBySide(st,p,q,sq,wlq,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  #paper=True,onecol=True)

def goWarpDifferencelog():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 181,-90 # sampling for wavelet H 
  niter = 100
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 1.00

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D2(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
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
  title="logpqwlqsq"
  tmin = 0
  tmax = nt-1
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  amax = 1.5
  hint = 1.0
  plot4TracesSideBySide(st,p,q,sq,wlq,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)

def goSimpleEstimate():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 1.00

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D2(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Warp g to f using new method
  warp = Warper()
  sg = warp.applyS(u,g)

  #Estimate Wavelet
  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 81,-40 # sampling for wavelet H 
  nd,kd = 81,-40 # sampling for wavelet H 
  niter = 300  
  sfac = 0.0
  tmin = tmin
  tmax = 275
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  ww.setLineSearchMinScale(0.0001);
  #First guesses of c and b. 
  bguess = zerofloat(nb)
  bguess[-kb] = 1.0
  cguess = ww.getWaveletC(nc,kc,nb,kb,bguess,u,f,g)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  dump(bw)
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  stepl = ww.getStepLength()
  rapprapp = ww.getRappRapp()
  deltamag = ww.getDeltaMag()
  condnum = ww.getCondNum()

  #Get known wavelet
  dk = getWavelet(freqd,decayd,nc,kc,mpd)
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Create processed g (CSBg)
  bg = ww.applyC(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)

  #Create shaped squeezed g (HSg)
  one = zerofloat(nb)
  one[-kb] = 1.0
  hw = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyC(nc,kc,hw,sg)

  #Normalize wavelets
  nck = normalizeMAAWOS(ck)
  ncw = normalizeMAAWOS(cw)
  nhw = normalizeMAAWOS(hw)
  SimplePlot.asPoints(sub(f,hsg))

  #Plotting
  title = "fgsg"
  dt = 0.004
  hint = 3.0
  amax = 6.0
  tmin = 0
  tmax = nt-1
  st = Sampling(nt,dt,0.0)
  plot3TracesSideBySide(st,f,g,sg,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)
  title = "simplefcsbghsg"
  plot3TracesSideBySide(st,f,csbg,hsg,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  title = "simplecwck"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)
  title = "simpleh"
  plotWavelets(Sampling(nc,dt,kc*dt),[nhw],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "simple3Meas"
  siter = Sampling(niter,1.0,0.0)
  itermin = 0
  itermax = 266
  plot3MeasSideBySide(siter,rfrf,deltamag,stepl,itermin,itermax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)

def goDiffEstimate():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.06,0.05
  mpd = True#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 1.00

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D2(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Warp g to f using new method
  warp = Warper()
  sg = warp.applyS(u,g)

  #Estimate Wavelet
  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 141,-70# sampling for wavelet H 
  nd,kd = 141,-70# sampling for wavelet H 
  niter = 500
  sfac = 1.0
  tmin = 0
  tmax = nt-1
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setLineSearchMinScale(0.0001);
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  bguess = zerofloat(nb)
  bguess[-kb] = 1.0
  cguess = ww.getWaveletC(nc,kc,nb,kb,bguess,u,f,g)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nd,kd)
  dump(bw)
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  stepl = ww.getStepLength()
  rapprapp = ww.getRappRapp()
  deltamag = ww.getDeltaMag()
  condnum = ww.getCondNum()
  SimplePlot.asPoints(rfrf)
  SimplePlot.asPoints(deltamag)
  SimplePlot.asPoints(stepl)

  #Get known wavelet
  dk = getWavelet(freqd,decayd,nd,kd,mpd)
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Create processed g (CSBg)
  bg = ww.applyC(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)

  #Create shaped squeezed g (HSg)
  one = zerofloat(nb)
  one[-kb] = 1.0
  hw = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyC(nc,kc,hw,sg)

  #Normalize wavelets
  ndk = normalizeMAAWOS(dk)
  ndw = normalizeMAAWOS(dw)
  nck = normalizeMAAWOS(ck)
  ncw = normalizeMAAWOS(cw)
  nhw = normalizeMAAWOS(hw)
  SimplePlot.asPoints(sub(f,hsg))

  #Plotting
  tmin = 0
  tmax = nt-1
  title = "fgsg"
  dt = 0.004
  hint = 3.0
  amax = 6.0
  st = Sampling(nt,dt,0.0)
  plot3TracesSideBySide(st,f,g,sg,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)
  title = "difffcsbghsg"
  plot3TracesSideBySide(st,f,csbg,hsg,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  title = "diffcwck"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)
  title = "diffdwdk"
  plotWavelets(Sampling(nd,dt,kd*dt),[ndw,ndk],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)
  title = "diffh"
  plotWavelets(Sampling(nc,dt,kc*dt),[nhw],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "diff3Meas"
  siter = Sampling(niter,1.0,0.0)
  itermin = 0
  itermax = 111
  plot3MeasSideBySide(siter,rfrf,deltamag,stepl,itermin,itermax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)


def goDiffNoiseEstimate():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.5
  nrmsg = nrmsf
  sfac = 1.00

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D2(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  #Warp g to f using new method
  warp = Warper()
  sg = warp.applyS(u,g)

  #Estimate Wavelet
  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 141,-70# sampling for wavelet H 
  nd,kd = 141,-70# sampling for wavelet H 
  niter = 300
  sfac = 1.0
  tmin = 0
  tmax = nt-1
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setLineSearchMinScale(0.0001);
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  bguess = zerofloat(nb)
  bguess[-kb] = 1.0
  cguess = ww.getWaveletC(nc,kc,nb,kb,bguess,u,f,g)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  dump(bw)
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  stepl = ww.getStepLength()
  rapprapp = ww.getRappRapp()
  deltamag = ww.getDeltaMag()
  condnum = ww.getCondNum()

  #Get known wavelet
  dk = getWavelet(freqd,decayd,nc,kc,mpd)
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Create processed g (CSBg)
  bg = ww.applyC(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)

  #Create shaped squeezed g (HSg)
  one = zerofloat(nb)
  one[-kb] = 1.0
  hw = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyC(nc,kc,hw,sg)

  #Normalize wavelets
  nck = normalizeMAAWOS(ck)
  ncw = normalizeMAAWOS(cw)
  nhw = normalizeMAAWOS(hw)
  SimplePlot.asPoints(sub(f,hsg))

  #Plotting
  tmin = 0
  tmax = nt-1
  title = "fgsg"
  dt = 0.004
  hint = 3.0
  amax = 7.0
  st = Sampling(nt,dt,0.0)
  plot3TracesSideBySide(st,f,g,sg,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)
  title = "diffnoisecsbghsg"
  plot3TracesSideBySide(st,f,csbg,hsg,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)

  title = "diffnoisecwck"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)
  title = "diffnoiseh"
  plotWavelets(Sampling(nc,dt,kc*dt),[nhw],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "diffnoise3Meas"
  siter = Sampling(niter,1.0,0.0)
  itermin = 0
  itermax = 81
  plot3MeasSideBySide(siter,rfrf,deltamag,stepl,itermin,itermax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)

def go2D():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 1.00

  #Create synthetic f and g.
  nx = 10
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn2D2(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nx,nt,ni,randomi,moreps)

  #Warp g to f using new method
  warp = Warper()
  sg = warp.applyS(u,g)

  #Estimate Wavelet
  #Wavelet estimation parameters
  nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 141,-70 # sampling for wavelet H 
  nd,kd = 141,-70 # sampling for wavelet H 
  niter = 300
  sfac = 1.0
  tmin = 25
  tmax = 275
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  bguess = zerofloat(nb)
  bguess[-kb] = 1.0
  cguess = ww.getWaveletC(nc,kc,nb,kb,bguess,u,f,g)
  cbw = ww.getWaveletCInverseB(nb,kb,bguess,nc,kc,cguess,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletC(nb,kb,bw,nc,kc)
  dump(bw)
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  stepl = ww.getStepLength()
  rapprapp = ww.getRappRapp()
  deltamag = ww.getDeltaMag()
  condnum = ww.getCondNum()

  #Get known wavelet
  dk = getWavelet(freqd,decayd,nc,kc,mpd)
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletC(nc,kc,ck,nb,kb)

  #Create processed g (CSBg)
  bg = ww.applyC(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyC(nc,kc,cw,sbg)

  #Create shaped squeezed g (HSg)
  one = zerofloat(nb)
  one[-kb] = 1.0
  hw = ww.getWaveletC(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyC(nc,kc,hw,sg)

  #Normalize wavelets
  nck = normalizeMAAWOS(ck)
  ncw = normalizeMAAWOS(cw)
  nhw = normalizeMAAWOS(hw)
  #SimplePlot.asPoints(sub(f,hsg))

  #Plotting
  title = "fgsg"
  dt = 0.004
  hint = 3.0
  amax = 6.0
  #st = Sampling(nt,dt,0.0)
  #plot3TracesSideBySide(st,f,g,sg,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  #paper=True,onecol=True)
  #title = "simplefcsbghsg"
  #plot3TracesSideBySide(st,f,csbg,hsg,tmin*dt,tmax*dt,hint,amax,title=title,pngDir=pngDir,\
  #paper=True,onecol=True)

  title = "simplecwck"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)
  title = "simpleh"
  plotWavelets(Sampling(nc,dt,kc*dt),[nhw],0.1,title=title,pngDir=pngDir,paper=True,
  onecol=True)

  title = "simple3Meas"
  siter = Sampling(niter,1.0,0.0)
  itermin = 0
  itermax = 266
  plot3MeasSideBySide(siter,rfrf,deltamag,stepl,itermin,itermax,title=title,pngDir=pngDir,\
  paper=True,onecol=True)

def chooseNoise():
  x0,nx = 260,100
  f,gNoNR,u = getSinoImage(x0,nx)
  du = computeBackDiff2D(u)
  SimplePlot.asPixels(u)
  SimplePlot.asPixels(du)

  halfwidth = [0,1,2,3,4,5,6,7,8,9,10]
  nhalfwidth = len(halfwidth)
  for i in range(0,nhalfwidth):
    if halfwidth[i]==0:
      print "No filtering"
      g = copy(gNoNR)
    else:
      ref = RecursiveExponentialFilter(halfwidth[i])
      g = zerofloat(len(gNoNR[0]),len(gNoNR))
      ref.apply2(gNoNR,g)
    sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
    sp.addPixels(g)
    sp.addTitle("half width = "+str(halfwidth[i]))

def goSinopecGNAlphaAnalysis():
  directory = None
  #directory = "./slides/pc0_021_10_0sinopec100/"
  #directory = "./slides/21_10_0sinopec100/"
  #directory = "./slides/5_2_3sinopec100/"
  #nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/7_3_3sinopec100/"
  #nb,kb = 7,-3#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/9_4_3sinopec100/"
  #nb,kb = 9,-4#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/21_10_3sinopec100/"
  #nb,kb = 21,-10#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/31_15_3sinopec100/"
  #nb,kb = 31,-15#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/23_11_3sinopec100/"
  #nb,kb = 23,-11#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/19_9_3sinopec100/"
  #nb,kb = 19,-9#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/17_8_3sinopec100/"
  #nb,kb = 17,-8#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/15_7_3sinopec100/"
  #nb,kb = 15,-7#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/13_6_3sinopec100/"
  #nb,kb = 13,-6#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/11_5_3sinopec100/"
  nb,kb = 11,-5#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/25_12_3sinopec100/"
  #nb,kb = 25,-12#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/27_13_3sinopec100/"
  #nb,kb = 27,-13#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/29_14_3sinopec100/"
  #nb,kb = 29,-14#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/33_16_3sinopec100/"
  #nb,kb = 33,-16#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/35_17_3sinopec100/"
  #nb,kb = 35,-17#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/37_18_3sinopec100/"
  #nb,kb = 37,-18#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/39_19_3sinopec100/"
  #nb,kb = 39,-19#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/41_20_3sinopec100/"
  #nb,kb = 41,-20#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/43_21_3sinopec100/"
  #nb,kb = 43,-21#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/45_22_3sinopec100/"
  #nb,kb = 45,-22#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/47_23_3sinopec100/"
  #nb,kb = 47,-23#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/49_24_3sinopec100/"
  #nb,kb = 49,-24#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/51_25_3sinopec100/"
  #nb,kb = 51,-25#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/53_26_3sinopec100/"
  #nb,kb = 53,-26#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/61_30_3sinopec100/"
  #nb,kb = 61,-30#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/71_35_3sinopec100/"
  #nb,kb = 71,-35#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/81_40_3sinopec100/"
  #nb,kb = 81,-40#sampling for inverse wavelet A #Note ka<=kc 
  #directory = "./slides/0sinopec/"
  #directory = "./slides/2sinopec/"
  #directory = "./slides/3sinopec/"
  #directory = "./slides/4sinopec/"
  #directory = "./slides/5sinopec/"
  #directory = "./slides/10sinopec/"
  #get sino images
  x0,nx = 250,5
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
  alphab = [0.0,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
  nalphab = len(alphab)
  for ialphab in range(0,nalphab):
    print "alphab = "+str(alphab[ialphab])
    niter = 300
    ww = WaveletWarpingCBGN()
    ww.setMinPercentChange(0.01)#units are percentage.
    ww.setTimeRange(tmin,tmax)
    ww.setAlphaB(alphab[ialphab])
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

    #Display only
    wwOld = WaveletWarping()
    wwOld.setTimeRange(tmin,tmax)
    lbg = wwOld.applyL(u,bg)
    ######

    sbg = warp.applyS(u,bg)
    csbg = ww.applyC(nc,kc,cw,sbg)
    hsg = ww.applyC(nc,kc,hw,warp.applyS(u,ww.applyC(nb,kb,bone,g)))

    g = ww.makeRms1(225,800,g)
    gNoNR = ww.makeRms1(225,800,gNoNR)
    f = ww.makeRms1(f)
    sbg = ww.makeRms1(sbg)
    sg = ww.makeRms1(sg)
    bg = ww.makeRms1(225,800,bg)
    lbg = ww.makeRms1(225,800,lbg)
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
      #GN Meas
    print "All RMS"
    dump(allResRmsAllS)
    print "alphab"
    dump(alphaBS)
    pngDir = directory
    #pngDir = None
    title = "All Rms Residuals"+str(alphab[ialphab])
    maxrmsri = max(allResRmsAllS)#0.15
    minrmsrf = 0.0#allResRmsAllS[lastIter+1]#0.0
    siter = Sampling(niter,1.0,0.0)
    color=[Color.BLACK,Color.RED]
    vlabel,vminmax,vint = "RMS of all residuals",[minrmsrf,maxrmsri],None
    hlabel,hminmax,hint = "Iterations",[0.0,lastIter+1],None
    plotting.plotMeasInSamePlot(siter, [allResRmsAllS],\
    color=color,\
    vlabel=vlabel, vminmax=vminmax, vint=vint,\
    hlabel=hlabel, hminmax=hminmax, hint=hint,\
    title=title, pngDir=pngDir,\
    slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
    print "RMS = "+str(allResRmsAllS[lastIter])
    print "RMS = "+str(allResRmsAllS[lastIter+1])
    print "RMS = "+str(allResRmsAllS[lastIter+2])

    """
    pngDir = directory
    #pngDir = None
    title = "alphab"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,
    [bPenaltyRes2NormSqFinaS,cPenaltyRes2NormSqFinaS,alphaBS,alphaCS],\
    color=None,\
    vlabel=["bPenRes","cPenRes","alphaB","alphaC"],\
    vminmax=[None,None,None,None,None,None],\
    vint=[None,None,None,None,None,None],\
    hlabel="Iterations",hminmax=[0.0,lastIter],hint=None,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)
    """


  """
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
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated c"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ncw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[ndw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Interpolated Guess d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "ImpulseGuess d"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  ndgw = zerofloat(nc)
  ndgw[-kc*4] = 1
  plotting.plotWavelets(st,[ndgw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  pngDir = directory
  #pngDir = None     
  title = "Guess (b)"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbguess],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  pngDir = directory
  #pngDir = None     
  title = "Estimated b"
  hint = None
  st = Sampling(nb,dt,kb*dti)
  plotting.plotWavelets(st,[nbw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

  title = "estimated PP and PS wavelets"
  pngDir = directory
  sh = Sampling(nc,dt,kc*dti)
  plotWaveletsPpPs(sh,ncw,ndw,0.9,0.8,16.0/9.0,title=title,pngDir=pngDir,twocol=True)

  title = "starting PP and PS wavelets"
  pngDir = directory
  plotWaveletsPpPs(sh,ncguess,ndgw,0.9,0.8,16.0/9.0,title=title,pngDir=pngDir,twocol=True)
  """



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
      pp.setVInterval(ix,tmark[ix])
  pp.setHLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setVisible(True)
  pf.setSize(800,900)
  if pngDir:
     pf.setFontSizeForPrint(5.0,222.0)
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
     pf.setFontSizeForPrint(5.0,222.0)
     pngDir = pngDir+title+"1wavelet"+"onecol.png"
     pf.paintToPng(720.0,3.08,pngDir)
  else:
    pp.setTitle(title)


def plotSequenceMeas(st,xs,amax=None,tmark=None,hlabel=None,pngDir=None,labels=None,title=None):
  nx = len(xs)
  pp = PlotPanel(nx,1)
  for ix,xi in enumerate(xs):
    pv = pp.addPoints(ix,0,st,xi)
    if labels:
      pp.setVLabel(ix,labels[ix])
    if amax:
      pp.setVLimits(ix,0,amax[ix])
    if tmark:
      pp.setVInterval(ix,tmark[ix])
  pp.setHLabel(hlabel)
  pf = PlotFrame(pp)
  pf.setVisible(True)
  pf.setSize(800,900)
  if pngDir:
     pf.setFontSizeForPrint(5.0,222.0)
     pngDir = pngDir+title+"1wavelet"+"onecol.png"
     pf.paintToPng(720.0,3.08,pngDir)
  else:
    pp.setTitle(title)



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


def plot3MeasSideBySide(st, f, g, h, itermin, itermax, title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):
  niter = len(f)
  pv1 = PointsView(st,f)
  pv1.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  pv2 = PointsView(st,g)
  pv2.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  pv3 = PointsView(st,h)
  pv3.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  minf = 10000000.0
  for iter in range(itermin,itermax):
    if f[iter]<minf:
      minf = f[iter]
  for iter in range(itermax,niter):
    f[iter] = minf
  maxf = int(max(f)*10.0+.999)/10.0
  minf = int(minf*10.0-.999)/10.0
  maxg = int(max(g)*10.0+.999)/10.0
  maxfd2 = (maxf-minf)/2.0
  maxgd2 = maxg/2.0
  print "minf = "+str(minf)
  print "maxf = "+str(maxf)
  
  pp = PlotPanel(3,1,PlotPanel.Orientation.X1RIGHT_X2UP,
  PlotPanel.AxesPlacement.LEFT_BOTTOM)
  pp.addTiledView(0,0,pv1)
  pp.addTiledView(1,0,pv2)
  pp.addTiledView(2,0,pv3)
  pp.setHLimits(itermin,itermax)
  pp.setVLimits(0,minf,maxf)
  pp.setVLimits(1,0,maxg)
  pp.setVLimits(2,0,1.0)
  pp.setVInterval(0,maxfd2)
  pp.setVInterval(1,maxgd2)
  pp.setVInterval(2,0.5)
  pp.setHLabel("Iterations")
  pp.setVLabel(0,"ssd")
  pp.setVLabel(1,"dmag")
  pp.setVLabel(2,"scale")
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

def testAddNoise():
  n = 100
  p = zerofloat(n)
  p[n/2] = 1
  nrms,seed = .1,5
  addNoise(nrms,seed,p)

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
  y = zerofloat(nf)
  maxv = max(f)
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

def getSinoImages(x0,nx):
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

def testWarpingandBandPassFilter():
  nt,ni = 480,1#number of time samples in p and q; number of random impulses in p and q.
  na,ka = 81,-40 #sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 181,-90 # sampling for wavelet H 
  v = 0.0#The amount of shift between p and q.
  r0,r1 = 2.0,1.0#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  freqc,decayc,= 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd,= 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.00
  nrmsg = nrmsf
  randomi = False 
  moreps = False 
  p,q,f,g,noisef,noiseg,u,tmin,tmax = createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
  ww = WaveletWarpingAEig()
  ww.setTimeRange(tmin,tmax)
  r = r0
  nu = len(u)
  nuu = int(r*(nu-1)+1)
  nq = len(q)
  nuq = int(r*(nu-1)+1)
  dtqu = 1.0/r
  r = r0
  w = 0.20/r
  lq = ww.applyL(r,u,q)
  ulq = ww.upSample(nq,1.0,0.0,lq,nuq,dtqu,0.0)
  uu = ww.upSampleLinear(nu,1.0,0.0,u,nuu,dtqu,0.0)
  uq = ww.upSample(nq,1.0,0.0,q,nuq,dtqu,0.0)
  suq  = ww.applyS(dtqu,uu,uq)
  lsuq = ww.applyL(r,uu,suq)
  dlsuq = ww.subSample(r,nu,lsuq)
  stp = Sampling(nu,1.0,0.0)
  stuq = Sampling(nuq,dtqu,0.0)
  stsuq = Sampling(nuu,dtqu,0.0)
  maxp = max(p)
  maxpd2 = maxp/2.0
  amax = [maxp,maxp,maxp,maxp,maxp]
  tmark = [maxpd2,maxpd2,maxpd2,maxpd2,maxpd2]
  plotSequences(stp,[p,q,dlsuq],amax=amax,tmark=tmark,\
  labels=["p","q","dlsuq"],title="pqdlsuq")
  amax = [maxp]
  tmark = [maxpd2]
  plotSequences(stuq,[uq],amax=amax,tmark=tmark,\
  labels=["uq"],title="uq")
  plotAmplitudeSpectrumT(stp, dlsuq, 0, nq, "amp DLSUq", amax=None)
  plotAmplitudeSpectrumT(stsuq, lsuq, 0, nuu, "amp LSUq", amax=None)
  plotAmplitudeSpectrumT(stsuq, suq, 0, nuu, "amp SUq", amax=None)
  plotAmplitudeSpectrumT(stuq, uq, 0, nuq, "amp Uq", amax=None)
  plotAmplitudeSpectrumT(stp, q, 0, nq, "amp q", amax=None)
  plotAmplitudeSpectrumT(stp, p, 0, nq, "amp p", amax=None)

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

def plotWaveletsPpPs(st,h1,h2,fracWidth,fracHeight,aspectRatio,title=None,
  pngDir=None,halfcol=None,onecol=None,twocol=None):
  wpt = 240
  pp = PlotPanel(2,1)
  h1 = mul(h1,1.0/max(abs(h1)))
  h2 = mul(h2,1.0/max(abs(h2)))
  sv1 = pp.addSequence(0,0,st,h1)
  sv2 = pp.addSequence(1,0,st,h2)
  pp.setVLimits(0,-1.05,1.05)
  pp.setVLimits(1,-1.05,1.05)
  if pngDir:
    pp.setVLabel(0,"Amplitude")
    pp.setVLabel(1,"Amplitude")
  else:
    pp.setVLabel(0,"PP wavelet")
    pp.setVLabel(1,"PS wavelet")

  pp.setHLabel("Time (s)")
  pf = PlotFrame(pp)
  pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if onecol:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"onecol.png"
      pf.paintToPng(720.0,3.08,pngDir)
    if twocol:
      print "zz = "+str(aspectRatio)
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"twocol.png"
      #pf.paintToPng(720.0,3.0,pngDir)
      pf.paintToPng(720.0,14.0,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)







  

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())



