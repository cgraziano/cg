#############################################################################
# Demo of 2 wavelet estimations from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.awt.ColorMap import *
from edu.mines.jtk.lapack import *
from wwarp import WaveletWarpingCA,Warper
from wwarp import ShapingFilter 
import synthetic
import plotting
from java.util import Random

############################################################################


def main(args):
  oneDNoNoiseHighSqueezingOneWavelets()

#A separate test to figure out if the Gauss-Newton method is working.
def oneDNoNoiseHighSqueezingOneWavelets():
  directory = None#"./slides/nonoiseHighSqueezing/"
  #Synthetic parameters
  nt,ni,randomi = 581,30,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  r0,r1 = 3.15,1.55#
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

  pngDir = directory
  #pngDir = None
  title= "[f,csag] [f,hsg]"
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
  [[f,csag],[f,hsg]],\
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

  title= "[f,CSAg] [f,f]"
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
  [[f,csag],[f,f]],\
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
  title= "f CSAgZoom"
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
  plotting.plotTracesSideBySide(st,[f,csag],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


  #pngDir = directory
  pngDir = None
  title= "f csag p sag ag q g"
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,nt*dt
  vlabel,vminmax,vint = "Time (s)",None,None
  hint1,hint2 = 0.5,4.0
  hmin1,hmin2 = -1.0,-8.0
  hmax1,hmax2 = 1.0,8.0
  hlabel = ["f","HSg","CSAg","du","p","SAg","Ag","q","g"]
  hminmax = [[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin2,hmax2],[hmin1,hmax1],\
  [hmin1,hmax1],[hmin1,hmax1],[hmin1,hmax1],[hmin2,hmax2]]
  hint = [None,None,None,None,None,None,None,None,None]
  hsize,vsize = 960,560
  tilespacing = 5
  plotting.plotTracesSideBySide(st,[f,hsg,csag,du,p,sag,ag,q,g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  tilespacing = tilespacing,\
  title=title,pngDir=pngDir,\
  paper=True,onecol=True)

    #Normalize
  nhw = normalizeM(hw)
  ncw = normalizeM(cw)
  naw = normalizeM(aw)
  nck = normalizeM(ck)
  nak = normalizeM(ak)
  nak = normalizeM(ak)
  dti = dt

  #"""
  #Wavelet interpolation
  error = 0.001
  freq = 0.49
  dt = 0.004
  scale = 4
  nc2 = scale*(nc-1)+1
  na2 = scale*(na-1)+1
  dt2 = dt/scale
  nck = interpolate(nc,kc,nck,dt,nc2,dt2,error,freq)
  ncw = interpolate(nc,kc,ncw,dt,nc2,dt2,error,freq)
  nhw = interpolate(nc,kc,nhw,dt,nc2,dt2,error,freq)
  nak = interpolate(na,ka,nak,dt,na2,dt2,error,freq)
  naw = interpolate(na,ka,naw,dt,na2,dt2,error,freq)
  nc = nc2
  na = na2
  dt = dt2
  #"""

    #Wavelets
  #pngDir = directory
  pngDir = None     
  title = "Shaping Filter (h)"
  hint = None
  st = Sampling(nc,dt,kc*dti)
  plotting.plotWavelets(st,[nhw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)

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
  title = "Estimated a"
  hint = None
  st = Sampling(na,dt,ka*dti)
  plotting.plotWavelets(st,[naw],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)
  title = "Known a"
  hint = None
  st = Sampling(na,dt,ka*dti)
  plotting.plotWavelets(st,[nak],dotsonsticks=True,hint=hint,title=title,pngDir=pngDir,slide=True,fracWidth=0.9,fracHeight=0.8,aspectRatio=16.0/9.0)


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





  #############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())


