#############################################################################
# Demo of 1 wavelet estimations from warping casued by time-to-depth conversions.

from imports import *

from wwarp import WaveletWarpingCA, Warper
from wwarp import ShapingFilter 
import synthetic, syntheticTimeDepth
import plotting

############################################################################

def main(args):
  convertTToDConstantV()

def testTimeToDepthConversion():
  #Synthetic parameters
  nt,ni,randomi = 581,2,False# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
  r0 = 2000#1st Velocity (m/s)
  r1 = r0#Final Velocity (m/s)
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.07
  freqd,decayd = freqc,decayc
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

def convertTToDConstantV():
  nt,ni,randomi = 581,2,False# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
  r0 = 2000#1st Velocity (m/s)
  r1 = r0#Final Velocity (m/s)
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.07
  freqd,decayd = freqc,decayc
  mpc = False#is wavelet in f mininmum phase?
  mpd = False#is wavelet in g mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = syntheticTimeDepth.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)

  pngDir = None
  dt = 0.004
  st = Sampling(nt,dt,0.0)
  vmin,vmax = 0*dt,481*dt
  hmin,hmax = -13.0,13.0
  vlabel,vminmax,vint = "Time (s)",[vmin,vmax],0.5
  hlabel,hminmax,hint = ["Amplitude"],[[hmin,hmax]],[5.0]
  hsize,vsize = 960,560
  title= "f one wavelet"
  plotting.plotTracesSideBySide(st,[f],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)

  pngDir = None
  nz = 481
  dz = 5
  sz = Sampling(nz,dz,0.0)
  vmin,vmax = 0*dz,481*dz
  hmin,hmax = -13.0,13.0
  vlabel,vminmax,vint = "Depth (m)",[vmin,vmax],None
  hlabel,hminmax,hint = ["Amplitude"],[[hmin,hmax]],[5.0]
  hsize,vsize = 960,560
  title= "g one wavelet"
  plotting.plotTracesSideBySide(st,[g],\
  vlabel=vlabel,vminmax=vminmax,vint=vint,\
  hlabel=hlabel,hminmax=hminmax,hint=hint,\
  title=title,pngDir=pngDir,\
  paper=True,twocol=True)




  


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

