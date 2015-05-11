#############################################################################
# Demo of 2 wavelet estimations from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.awt.ColorMap import *
from edu.mines.jtk.lapack import *
from wwarp import WaveletWarpingCBGN, WaveletWarpingCBGNOld, Warper
from wwarp import ShapingFilter 
import synthetic
import plotting
from java.util import Random

############################################################################

def main(args):
  #goSFacTest()
  #go1Dvs2DTest()
  #goNoiseTest()
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
  #goSino1000ShapingFilter()
  #goTestBrentMinFinder()
  #goWarpTest()
  #goEstimateWaveletTest()
  #goEstimateWaveletTest2D()
  goSinopec1()
  #goSinopec()
  #goVaryUPTest()
  #goSimpleTest()
  #goBuild2DSynthetics()
  #goTestUniqMeas()
  #goNaVarySyntheticTest()
  #plotNaVarySyntheticTest()
  #goTestU()
  #goTestTMaxTheory()
  #plotTestTMaxTheory()
  #goAcVsAa()
  #goInverseLengthScan()
  #goFreqPlot()
  #printSumSqDiff()jj
  #testSLG(nas[i],kas[i])
  #testAddNoise()
  #testWarpingandBandPassFilter()
 
#This test is designed to show that increasing the sfac parameter will not change the overall
#result.
def goSFacTest():
  #####1D##########
  #Synthetic parameters
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 81,-40 # sampling for wavelet H 
  niter = 500
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
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
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
  ww.setMaxPercentChange(0.000)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  bk = ww.getWaveletH(nc,kc,ck,nb,kb)

  #Normalize wavelets
  ncw = normalizeM(cw)
  nck = normalizeM(ck)
  nbw = normalizeM(bw)
  nbk = normalizeM(bk)
  ndw = normalizeM(dw)
  ndk = normalizeM(dk)
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
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyH(nc,kc,h,sg)

  #deconvolve f
  bf = ww.applyH(nb,kb,bw,f)
  
  #Create warped g
  bg = ww.applyH(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,cw,sbg)

  #plotting
  pngDir = None
  if (write):
    pngDir = "./png/estimatewavelet/"

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title1 = "c"
  title2 = "d"
  title3 = "b"
  title4 = "h"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[h],title=title4,pngDir=pngDir,onecol=True)

  if len(rfrf)>0:
    maxriri = max(riri)
    title = "rms ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    maxstepl = max(stepl)
    title = "Scale"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[stepl],amax=maxstepl,hlabel="Iterations",labels="Scale",pngDir=None,title=title)

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
  title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,csbg,hsg,g,sg,sub(hsg,f)],amax=amax,tmark=tmark,\
  labels=["f","CSGg","HSg","g","Sg","hsg-f"],pngDir=pngDir,\
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
def goTestNcNbTradeOff():
  #Synthetic parameters
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
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
    ww.setStabilityFactor(sfac)
    ww.setMaxPercentChange(0.001)
    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
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
    ww.setStabilityFactor(sfac)
    ww.setMaxPercentChange(0.001)
    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
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
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
  ww.setMaxPercentChange(0.001)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  bk = ww.getWaveletH(nc,kc,ck,nb,kb)

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
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
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
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 77,-38# sampling for wavelet H 
  niter = 500
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
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
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
  ww.setMaxPercentChange(0.001)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  bk = ww.getWaveletH(nc,kc,ck,nb,kb)

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
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)
  nhw = normalizeMAAWOS(h)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
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
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  nb,kb = 45,-22#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 37,-18# sampling for wavelet H 
  niter = 500
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
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
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
  ww.setMaxPercentChange(0.001)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  bk = ww.getWaveletH(nc,kc,ck,nb,kb)

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
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
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
    ww.setStabilityFactor(sfac)
    ww.setMaxPercentChange(0.001)
    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  ww.setStabilityFactor(sfac)
  ww.setLineSearchMinScale(0.0)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
  one = zerofloat(nb)
  one[-kb] = 1
  hw = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
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
  vmin,vmax = 0*dt,481*dt
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
    ww.setStabilityFactor(sfac)
    ww.setMaxPercentChange(0.001)
    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
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
    ww.setStabilityFactor(sfac)
    ww.setMaxPercentChange(maxpc)
    #Shaping filter
    b = zerofloat(nb)
    b[-kb] = 1.0
    h = ww.getWaveletH(nh,kh,nb,kb,b,u,f,g)
    rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,b,u,f,g))


    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  vmin,vmax = 0*dt,481*dt
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
    ww.setStabilityFactor(sfac)
    ww.setLineSearchMinScale(minsl)
    #Shaping filter
    b = zerofloat(nb)
    b[-kb] = 1.0
    h = ww.getWaveletH(nh,kh,nb,kb,b,u,f,g)
    rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,b,u,f,g))

    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
    ww.setStabilityFactor(sfac)
    ww.setMaxPercentChange(pc)
    #Shaping filter
    b = zerofloat(nb)
    b[-kb] = 1.0
    h = ww.getWaveletH(nh,kh,nb,kb,b,u,f,g)
    rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,b,u,f,g))


    """
    #First guesses of c and b. 
    if (i==0):
      b = zerofloat(nb)
      b[-kb] = 1.0
      c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
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
    c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
    title = "Increase start regular"

    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletH(nb,kb,bw,nc,kc)

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
    ww.setStabilityFactor(sfac)
    ww.setMaxPercentChange(0.001)
    #Shaping filter
    b = zerofloat(nb)
    b[-kb] = 1.0
    h = ww.getWaveletH(nh,kh,nb,kb,b,u,f,g)
    rmsrh[i] = ww.rms(ww.computeResidual(nh,kh,h,nb,kb,b,u,f,g))


    #First guesses of c and b. 
    b = zerofloat(nb)
    b[-kb] = 1.0
    c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
    cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
    cw = cbw[0]
    bw = cbw[1]
    dw = ww.getWaveletH(nb,kb,bw,nc,kc)

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

    
#In this test, a synthetic trace is created. This trace is then replicated to create a 2D
#synthetic. This 1D and 2D synthetic should have the same system of equations
def go1Dvs2DTest():
  #####1D##########
  #Synthetic parameters
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 81,-40 # sampling for wavelet H 
  niter = 500
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 0.00

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
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  bk = ww.getWaveletH(nc,kc,ck,nb,kb)

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
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyH(nc,kc,h,sg)

  #deconvolve f
  bf = ww.applyH(nb,kb,bw,f)
  
  #Create warped g
  bg = ww.applyH(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,cw,sbg)

  #plotting
  pngDir = None
  if (write):
    pngDir = "./png/estimatewavelet/"

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title1 = "c"
  title2 = "d"
  title3 = "b"
  title4 = "h"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[h],title=title4,pngDir=pngDir,onecol=True)

  if len(rfrf)>0:
    maxriri = max(riri)
    title = "rms ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    maxstepl = max(stepl)
    title = "Scale"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[stepl],amax=maxstepl,hlabel="Iterations",labels="Scale",pngDir=None,title=title)

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
  title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,csbg,hsg,g,sg,sub(hsg,f)],amax=amax,tmark=tmark,\
  labels=["f","CSGg","HSg","g","Sg","hsg-f"],pngDir=pngDir,\
  title=title)
  #plotSequences(st,[bf,sbg],amax=[1,1],tmark=[0.1,0.1],\
  #labels=["Bf","SBg"],pngDir=pngDir,\
  #title=title)

  #####2D#####
  nx = 100
  p2,q2,f2,g2,noisef2,noiseg2,u2,ntmin,ntmax = synthetic.createSyntheticLn2D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,\
  nrmsf,nrmsg,nx,nt,ni,randomi,moreps)

  #Estimate wavelet
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u2,f2,g2)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u2,f2,g2,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
  dump(bw)

  #Get iteration information
  riri2 = ww.getRiRi()
  rfrf2 = ww.getRfRf()
  vv2 = ww.getVV()
  stepl2 = ww.getStepLength()
  deltamag2 = ww.getDeltaMag()
  condnum2 = ww.getCondNum()
  SimplePlot.asPoints(sub(rfrf2,rfrf));
  #SimplePlot.asPoints(sub(f2[0],f2[1]));
  #SimplePlot.asPoints(sub(f2[0],f2[2]));
  #SimplePlot.asPoints(sub(f2[0],f2[3]));
  #SimplePlot.asPoints(sub(f2[0],f2[4]));
  #SimplePlot.asPoints(sub(g2[0],g2[1]));
  #SimplePlot.asPoints(sub(g2[0],g2[2]));
  #SimplePlot.asPoints(sub(g2[0],g2[3]));
  #SimplePlot.asPoints(sub(g2[0],g2[4]));

  #Get known wavelet
  dk = getWavelet(freqd,decayd,nc,kc,mpd)
  ck = getWavelet(freqc,decayc,nc,kc,mpc)
  bk = ww.getWaveletH(nc,kc,ck,nb,kb)

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
  sg2 = warp.applyS(u2,g2)
  sq2 = warp.applyS(u2,q2)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletH(nc,kc,nb,kb,one,u2,f2,sg2)
  hsg2 = ww.applyH(nc,kc,h,sg2)


  #deconvolve f
  bf2 = ww.applyH(nb,kb,bw,f2)
  
  #Create warped g
  bg2 = ww.applyH(nb,kb,bw,g2)
  sbg2 = warp.applyS(u2,bg2)
  csbg2 = ww.applyH(nc,kc,cw,sbg2)

  #plotting
  pngDir = None
  if (write):
    pngDir = "./png/estimatewavelet/"

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title1 = "c 2D"
  title2 = "d 2D"
  title3 = "b 2D"
  title4 = "h 2D"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[h],title=title4,pngDir=pngDir,onecol=True)

  if len(rfrf)>0:
    maxriri = max(riri2)
    title = "rms ri (black) rms rf (blue) 2D"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri2,rfrf2],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    maxstepl = max(stepl2)
    title = "Scale 2D"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[stepl2],amax=maxstepl,hlabel="Iterations",labels="Scale",pngDir=None,title=title)

    maxvv= max(vv2)
    title = "2 norm of neg. gradient square 2D"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[vv2],amax=maxvv,tmark=maxvv/100.0,hlabel="Iterations",labels="V'V",pngDir=None,title=title)


    maxdm = max(deltamag2)
    title = "Delta Mag 2D"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[deltamag2],amax=maxdm,tmark=maxdm/100.0,hlabel="Iterations",labels="Delta Mag",pngDir=None,title=title)

    maxcn = max(condnum2)
    title = "X Condition Number 2D"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[condnum2],amax=maxcn,tmark=maxcn/100.0,hlabel="Iterations",labels="X Condition Number",pngDir=None,title=title)
  #ssdzoom,fiter = zoomConverge(ssd,1.0)
  #nsubiter = len(ssdzoom)
  #ssub = Sampling(nsubiter,1.0,fiter)
  #med = getMedian(ssdzoom)
  #ssdzoom = sub(ssdzoom,med)
  #sp = SimplePlot()
  #sp.addPoints(ssub,ssdzoom)
  #sp.addTitle("Sum of square differences between f and CSBg")

  maxf = max(f2[nx-1])
  maxfd2 = maxf/2.0
  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u2[nx-1],dt)
  title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f2[nx-1],csbg2[nx-1],hsg2[nx-1],g2[nx-1],sg2[nx-1],sub(hsg2[nx-1],f2[nx-1])],\
  amax=amax,tmark=tmark,labels=["f","CSGg","HSg","g","Sg","hsg-f"],pngDir=pngDir,\
  title=title+" nx-1")
  #plotSequences(st,[bf,sbg],amax=[1,1],tmark=[0.1,0.1],\
  #labels=["Bf","SBg"],pngDir=pngDir,\
  #title=title)

#In this test I examine the effect of noise on the results of the Gauss-Newton method.
#I am only examining the 1D case.
def goNoiseTest():

  #####1D No Noise##########
  #Synthetic parameters
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
  nb,kb = 5,-2#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 81,-40 # sampling for wavelet H 
  niter = 500
  r0,r1 = 2.3,1.3#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = 0.0#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf
  sfac = 0.00

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
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)

  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
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
  bk = ww.getWaveletH(nc,kc,ck,nb,kb)

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
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,sg)
  hsg = ww.applyH(nc,kc,h,sg)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
  hsg = ww.applyH(nc,kc,h,sg)

  #deconvolve f
  bf = ww.applyH(nb,kb,bw,f)
  
  #Create warped g
  bg = ww.applyH(nb,kb,bw,g)
  sbg = warp.applyS(u,bg)
  csbg = ww.applyH(nc,kc,cw,sbg)

  #plotting
  pngDir = None
  if (write):
    pngDir = "./png/estimatewavelet/"

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  title1 = "c No Noise"
  title2 = "d No Noise"
  title3 = "b No Noise"
  title4 = "h No Noise"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[h],title=title4,pngDir=pngDir,onecol=True)

  if len(rfrf)>0:
    maxriri = max(riri)
    title = "rms ri (black) rms rf (blue) No Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    maxstepl = max(stepl)
    title = "Scale"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[stepl],amax=maxstepl,hlabel="Iterations",labels="Scale",pngDir=None,title=title)

    maxvv= max(vv)
    title = "2 norm of neg. gradient square No Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[vv],amax=maxvv,tmark=maxvv/100.0,hlabel="Iterations",labels="V'V",pngDir=None,title=title)


    maxdm = max(deltamag)
    title = "Delta Mag No Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[deltamag],amax=maxdm,tmark=maxdm/100.0,hlabel="Iterations",labels="Delta Mag",pngDir=None,title=title)

    maxcn = max(condnum)
    title = "X Condition Number No Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[condnum],hlabel="Iterations",labels="X Condition Number",pngDir=None,title=title)
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
  title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,csbg,hsg,g,sg,sub(hsg,f)],amax=amax,tmark=tmark,\
  labels=["f","CSGg","HSg","g","Sg","hsg-f"],pngDir=pngDir,\
  title=title)
  #plotSequences(st,[bf,sbg],amax=[1,1],tmark=[0.1,0.1],\
  #labels=["Bf","SBg"],pngDir=pngDir,\
  #title=title)


  #####1D Noise (no)##########
  #Synthetic parameters
  nrmsf = 0.5
  nrmsg = nrmsf

  #Create synthetic f and g.
  pno,qno,fno,gno,noisef,noiseg,uno,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
  freqd,decayd,mpd,r0,r1,v,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  duno = computeBackDiff(uno)

  #Estimate wavelet
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  #First guesses of c and b. 
  bno = zerofloat(nb)
  bno[-kb] = 1.0
  cno = ww.getWaveletH(nc,kc,nb,kb,bno,uno,fno,gno)

  cbwno = ww.getWaveletCInverseB(nb,kb,b,nc,kc,cno,uno,fno,gno,niter)
  cwno = cbwno[0]
  bwno = cbwno[1]
  dwno = ww.getWaveletH(nb,kb,bwno,nc,kc)
  dump(bwno)

  #Get iteration information
  ririno = ww.getRiRi()
  rfrfno = ww.getRfRf()
  vvno = ww.getVV()
  steplno = ww.getStepLength()
  deltamagno = ww.getDeltaMag()
  condnumno = ww.getCondNum()

  #Get known wavelet
  dkno = getWavelet(freqd,decayd,nc,kc,mpd)
  ckno = getWavelet(freqc,decayc,nc,kc,mpc)
  bkno = ww.getWaveletH(nc,kc,ckno,nb,kb)

  #Normalize wavelets
  ncwno = normalizeMAAWOS(cwno)
  nckno = normalizeMAAWOS(ckno)
  nbwno = normalizeMAAWOS(bwno)
  nbkno = normalizeMAAWOS(bkno)
  ndwno = normalizeMAAWOS(dwno)
  ndkno = normalizeMAAWOS(dkno)

  #Create simply warped g
  sgno = warp.applyS(uno,gno)
  sqno = warp.applyS(uno,qno)

  #Create shaped Sg
  one = zerofloat(nb)
  one[-kb] = 1
  hno = ww.getWaveletH(nc,kc,nb,kb,one,uno,fno,sgno)
  hsgno = ww.applyH(nc,kc,hno,sgno)

  #deconvolve f
  bfno = ww.applyH(nb,kb,bwno,fno)
  
  #Create warped g
  bgno = ww.applyH(nb,kb,bwno,gno)
  sbgno = warp.applyS(uno,bgno)
  csbgno = ww.applyH(nc,kc,cwno,sbgno)

  #plotting
  pngDir = None
  if (write):
    pngDir = "./png/estimatewavelet/"

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(uno,dt)
  title1 = "c Noise"
  title2 = "d Noise"
  title3 = "b Noise"
  title4 = "h Noise"
  plotWavelets(Sampling(nc,dt,kc*dt),[ncwno,nckno],title=title1,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[ndwno,ndkno],title=title2,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nb,dt,kb*dt),[nbwno,nbkno],title=title3,pngDir=pngDir,onecol=True)
  plotWavelets(Sampling(nc,dt,kc*dt),[hno],title=title4,pngDir=pngDir,onecol=True)

  if len(rfrf)>0:
    maxriri = max(ririno)
    title = "rms ri (black) rms rf (blue) Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[ririno,rfrfno],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    maxstepl = max(steplno)
    title = "Scale Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[steplno],amax=maxstepl,hlabel="Iterations",labels="Scale",pngDir=None,title=title)

    maxvv= max(vvno)
    title = "2 norm of neg. gradient square Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[vvno],amax=maxvv,tmark=maxvv/100.0,hlabel="Iterations",labels="V'V",pngDir=None,title=title)


    maxdm = max(deltamagno)
    title = "Delta Mag Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[deltamagno],amax=maxdm,tmark=maxdm/100.0,hlabel="Iterations",labels="Delta Mag",pngDir=None,title=title)

    maxcn = max(condnumno)
    title = "X Condition Number Noise"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[condnumno],hlabel="Iterations",labels="X Condition Number",pngDir=None,title=title)
  #ssdzoom,fiter = zoomConverge(ssd,1.0)
  #nsubiter = len(ssdzoom)
  #ssub = Sampling(nsubiter,1.0,fiter)
  #med = getMedian(ssdzoom)
  #ssdzoom = sub(ssdzoom,med)
  #sp = SimplePlot()
  #sp.addPoints(ssub,ssdzoom)
  #sp.addTitle("Sum of square differences between f and CSBg")

  maxf = max(fno)
  maxfd2 = maxf/2.0
  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  ut = mul(uno,dt)
  title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0)+"r1="+\
  str(r1)+"v="+str(v)+"nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)
  amax = [maxf,maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[fno,csbgno,hsgno,gno,sgno,sub(hsgno,fno)],amax=amax,tmark=tmark,\
  labels=["f","CSGg","HSg","g","Sg","hsg-f"],pngDir=pngDir,\
  title=title)
  #plotSequences(st,[bf,sbg],amax=[1,1],tmark=[0.1,0.1],\
  #labels=["Bf","SBg"],pngDir=pngDir,\
  #title=title)

  
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
  ww.setStabilityFactor(sfac)
  #Create shaped squeezed g (HSg)
  one = zerofloat(nb)
  one[-kb] = 1.0
  hw = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
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








def goWarpTest():
  #Synthetic parameters
  nt,ni,randomi = 481,2,False# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
  nb,kb = 81,-40#sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 81,-40 # sampling for wavelet H 
  niter = 500
  r0,r1 = [2.3],[1.3]#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = [0.0]#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = [0.0]
  nrmsg = nrmsf
  sfac = 1.00

  nv = len(v)
  nr0 = len(r0)
  nr1 = len(r1)
  nnoise = len(nrmsf)
  for ir0 in range(0,nr0):
    for ir1 in range(0,nr1):
      for iv in range(0,nv):
        for inoise in range(0,nnoise):
          #Create synthetic f and g.
          p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
          freqd,decayd,mpd,r0[ir0],r1[ir1],v[iv],nrmsf[inoise],nrmsg[inoise],\
          nt,ni,randomi,moreps)
          du = computeBackDiff(u)
          print "tmin = "+str(tmin)
          print "tmax = "+str(tmax)
          tmin = 0#tmin
          tmax = nt-1#tmax
          print "tmin = "+str(tmin)
          print "tmax = "+str(tmax)

          
          #Test warper.
          warp = Warper()
          warp.plotAmplitudeSpectrums(True)
          sq = warp.applyS(u,q)

          #p
          st = Sampling(nt,1.0,0.0)
          sp = SimplePlot()
          sp.addPoints(st,p)
          sp.setTitle("p")

          SimplePlot.asPoints(sub(p,sq))





def goTestBrentMinFinder():
  ww = WaveletWarpingCBGN()
  ww.testSimple(); 

def goEstimateWaveletTest():
  #Synthetic parameters
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  nb,kb = 2,0 #sampling for inverse wavelet B #Note ka<=kc (B is in g)
  nc,kc = 2,0 # sampling for wavelet H 
  niter = 1000
  r0,r1 = [2.3],[1.3]#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  v = [0.0]#The amount of shift between p and q.
  freqc,decayc = 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd = 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = [0.0]
  nrmsg = nrmsf
  sfac = 0.0001

  nv = len(v)
  nr0 = len(r0)
  nr1 = len(r1)
  nnoise = len(nrmsf)
  for ir0 in range(0,nr0):
    for ir1 in range(0,nr1):
      for iv in range(0,nv):
        for inoise in range(0,nnoise):
          #Create synthetic f and g.
          p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn1D(freqc,decayc,mpc,\
          freqd,decayd,mpd,r0[ir0],r1[ir1],v[iv],nrmsf[inoise],nrmsg[inoise],\
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
          ww.setStabilityFactor(sfac)
          #First guesses of c and b. 
          b = zerofloat(nb)
          b[-kb] = 1.0
          c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
    
          cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
          cw = cbw[0]
          bw = cbw[1]
          dw = ww.getWaveletH(nb,kb,bw,nc,kc)
          dump(bw)
          riri = ww.getRiRi()
          rfrf = ww.getRfRf()
          vv = ww.getVV()
          stepl = ww.getStepLength()
          deltamagb = ww.getDeltaMagB()
          deltamagc = ww.getDeltaMagC()
          condnum = ww.getCondNum()
          twonormb = ww.getTwoNormb()
          twonormc = ww.getTwoNormc()

          #Get known wavelet
          dk = getWavelet(freqd,decayd,nc,kc,mpd)
          ck = getWavelet(freqc,decayc,nc,kc,mpc)
          bk = ww.getWaveletH(nc,kc,ck,nb,kb)

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
          h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,sg)
          hsg = ww.applyH(nc,kc,h,sg)

          #Create shaped Sg
          one = zerofloat(nb)
          one[-kb] = 1
          h = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
          hsg = ww.applyH(nc,kc,h,sg)

          #deconvolve f
          bf = ww.applyH(nb,kb,bw,f)
          
          #Create warped g
          bg = ww.applyH(nb,kb,bw,g)
          sbg = warp.applyS(u,bg)
          csbg = ww.applyH(nc,kc,cw,sbg)
          #print "The sum of the square differences is "+str(ww.getSumSqDiff(na))+", which "+"is equivalent to the smallest eigenvalue"+str(ww.getEigVal(0))
          #print "The sum of square differences if an impulse is "+\
          #"assumed to the be the inverse wavelet"+str(ww.getSumSqDiffNoWaveletEst(kb))
          #getWarpingPartsSpectrums(r0[ir0],tmin,tmax,f,g,u)
          #plotAmplitudeSpectrumT(Sampling(nb,1.0,kb*1.0), nbw, 0, nb, "amp b", amax=None)
          #plotAmplitudeSpectrumT(Sampling(nc,1.0,kc*1.0), ncw, 0, nc, "amp c", amax=None)
          #plotAmplitudeSpectrumT(Sampling(nt,1.0,0.0), bf, 0, nt, "amp bf", amax=None)
          #plotAmplitudeSpectrumT(Sampling(nt,1.0,0.0), sbg, 0, nt, "amp sbg", amax=None)
          #plotAmplitudeSpectrumT(Sampling(nt,1.0,0.0), noisef, 0, nt, "amp noisef", amax=None)
 
          #plotting
          pngDir = None
          #pngDir = "./png/estimatewavelet/"

          dt,ft = 0.004,0.000
          st = Sampling(nt,dt,ft)
          ut = mul(u,dt)
          #title = "nt="+str(nt)+"ni="+str(ni)+\
          #"r0="+str(r0[ir0])+"r1="+str(r1[ir1])+"v="+str(v[iv])+\
          #"nrmsf="+str(nrmsf[inoise])+"nrmsg="+str(nrmsg[inoise])
          title1 = "c"
          title2 = "d"
          title3 = "b"
          title4 = "h"
          plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
          plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
          plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
          plotWavelets(Sampling(nc,dt,kc*dt),[h],title=title4,pngDir=pngDir,onecol=True)

          if len(rfrf)>0:

            lastiter = ww.getLastIter()
            maxriri = max(riri)
            title = "rms ri (black) rms rf (blue)"
            si = Sampling(niter,1.0,0.0)
            plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

            title = "Meas1"
            siter = Sampling(niter,1.0,0.0)
            plotting.plotMeasOnTopOfEachOther(siter,[stepl,vv],\
            color=None,\
            vlabel=["StepLength","2NormGradSq"],\
            vminmax=[None,None],\
            vint=[None,None],\
            hlabel="Iterations",hminmax=[0,lastiter],hint=None,\
            title=title,pngDir=pngDir,\
            slide=None,fracWidth=None,fracHeight=None,\
            paper=True,onecol=True,twocol=False)

            title = "Meas2"
            siter = Sampling(niter,1.0,0.0)
            plotting.plotMeasOnTopOfEachOther(siter,[deltamagb,deltamagc,twonormb,twonormc,condnum],\
            color=None,\
            vlabel=["DeltMagb","DeltaMagc","2Normb","2Normc","CondNum"],\
            vminmax=[None,None,None,None,None],\
            vint=[None,None,None,None,None],\
            hlabel="Iterations",hminmax=[0,lastiter],hint=None,\
            title=title,pngDir=pngDir,\
            slide=None,fracWidth=None,fracHeight=None,\
            paper=True,onecol=True,twocol=False)

          maxf = max(f)
          maxfd2 = maxf/2.0
          dt,ft = 0.004,0.000
          st = Sampling(nt,dt,ft)
          ut = mul(u,dt)
          title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0[ir0])+"r1="+\
          str(r1[ir1])+"v="+str(v[iv])+"nrmsf="+str(nrmsf[inoise])+"nrmsg="+str(nrmsg[inoise])
          amax = [maxf,maxf,maxf,maxf,maxf,maxf]
          tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
          plotSequences(st,[f,csbg,hsg,g,sg,sub(hsg,f)],amax=amax,tmark=tmark,\
          labels=["f","CSGg","HSg","g","Sg","hsg-f"],pngDir=pngDir,\
          title=title)
          #plotSequences(st,[bf,sbg],amax=[1,1],tmark=[0.1,0.1],\
          #labels=["Bf","SBg"],pngDir=pngDir,\
          #title=title)

def goEstimateWaveletTest2D():
  #Synthetic parameters
  nt,nx,ni,randomi = 481,100,2,False# number of time samples in p and q; number of random impulses in p and q.
  moreps = False 
  v = [0.0]#The amount of shift between p and q.
  r0,r1 = [2.3],[1.3]#constant u'(t)
  freqc,decayc,= 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  freqd,decayd,= 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = nrmsf

  nv = len(v)
  nr0 = len(r0)
  nr1 = len(r1)
  for ir0 in range(0,nr0):
    for ir1 in range(0,nr1):
      for iv in range(0,nv):
        #Create synthetic f and g.
        p,q,f,g,noisef,noiseg,u,tmin,tmax = synthetic.createSyntheticLn2D(freqc,decayc,mpc,\
        freqd,decayd,mpd,r0[ir0],r1[ir1],v[iv],\
        nrmsf,nrmsg,nx,nt,ni,randomi,moreps)
        tmin = tmin
        tmax = tmax+50
        print "tmin = "+str(tmin)
        print "tmax = "+str(tmax)


        #Estimate wavelet
        nb,kb = 5,-2#sampling for inverse wavelet A #Note ka<=kc 
        nc,kc = 101,-50 # sampling for wavelet C
        nd,kd = 101,-50 # sampling for wavelet D 
        niter = 500
        sfac = 0.0
        warp = Warper()
        ww = WaveletWarpingCBGN()
        ww.setParallel(True)
        ww.setMaxPercentChange(0.001)#units are percentage.
        ww.setTimeRange(tmin,tmax)
        ww.setStabilityFactor(sfac)
        ww.setLineSearchMinScale(0.0001)
        #First guesses of c and b. 
        b = zerofloat(nb)
        b[-kb] = 1.0
        c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
        cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
        cw = cbw[0]
        bw = cbw[1]
        dw = ww.getWaveletH(nb,kb,bw,nc,kc)
        one = zerofloat(nb)
        one[-kb] = 1
        hw = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
        riri = ww.getRiRi()
        rfrf = ww.getRfRf()
        vv = ww.getVV()
        stepl = ww.getStepLength()
        deltamag = ww.getDeltaMag()
        condnum = ww.getCondNum()

        #Get known wavelet
        dk = getWavelet(freqd,decayd,nc,kc,mpd)
        ck = getWavelet(freqc,decayc,nc,kc,mpc)
        bk = ww.getWaveletH(nc,kc,ck,nb,kb)



        #Normalize wavelets
        nhw = normalizeMAAWOS(hw)
        ncw = normalizeMAAWOS(cw)
        nck = normalizeMAAWOS(ck)
        nbw = normalizeMAAWOS(bw)
        ndw = normalizeMAAWOS(dw)
        nbk = normalizeMAAWOS(bk)
        ndk = normalizeMAAWOS(dk)

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

        #rms check
        """
        hsgmf = sub(hsg,f)
        sgmf = sub(sg,f)
        rmsf = ww.rms(f)
        rmssg = ww.rms(sg)
        rmshsg = ww.rms(hsg)
        rmshsgmf = ww.rms(hsgmf)
        rmssgmf = ww.rms(sgmf)
        print "rmsf = "+str(rmsf)
        print "rmssg = "+str(rmssg)
        print "rmshsg = "+str(rmshsg)
        print "rmshsgmf = "+str(rmshsgmf)
        print "rmssgmf = "+str(rmssgmf)
        sp1 = SimplePlot()
        sp2 = SimplePlot()
        sp3 = SimplePlot()
        sp4 = SimplePlot()
        sp5 = SimplePlot()
        sp1.setTitle("f")
        sp2.setTitle("sg")
        sp3.setTitle("hsg")
        sp4.setTitle("hsg-f")
        sp5.setTitle("sg-f")
        sp1.setHLimits(tmin,tmax)
        sp2.setHLimits(tmin,tmax)
        sp3.setHLimits(tmin,tmax)
        sp4.setHLimits(tmin,tmax)
        sp5.setHLimits(tmin,tmax)
        sp1.addPoints(f)
        sp2.addPoints(sg)
        sp3.addPoints(hsg)
        sp4.addPoints(hsgmf)
        sp5.addPoints(sgmf)
        """

        #plotting
        pngDir = None
        if (write):
          pngDir = "./png/sino/"

        dt = 0.004
        ft = 0.000
        st = Sampling(nt,dt,ft)
        ut = mul(u,dt)
        title1 = "c"
        title2 = "d"
        title3 = "b"
        title4 = "h"
        plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
        plotWavelets(Sampling(nc,dt,kc*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
        plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title3,pngDir=pngDir,onecol=True)
        plotWavelets(Sampling(nc,dt,kc*dt),[hw],title=title4,pngDir=pngDir,onecol=True)
        du = computeBackDiff2D(u)

        getWarpingPartsSpectrums(4.0,tmin,tmax,f[0],g[0],u[0])
        plotAmplitudeSpectrumT(Sampling(nb,1.0,kb*1.0), nbw, 0, nb, "amp b", amax=None)
        plotAmplitudeSpectrumT(Sampling(nc,1.0,kc*1.0), ncw, 0, nc, "amp c", amax=None)

        if len(rfrf)>0:
          maxriri = max(riri)
          title = "rms ri (black) rms rf (blue)"
          si = Sampling(niter,1.0,0.0)
          plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

          maxstepl = max(stepl)
          title = "Scale"
          si = Sampling(niter,1.0,0.0)
          plotSequenceMeas(si,[stepl],amax=maxstepl,hlabel="Iterations",labels="Scale",pngDir=None,title=title)

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
          
        
        maxf = max(f)
        maxfd2 = maxf/2.0
        ntg = len(g[0])

        dt,ft = 0.004,0.000
        st = Sampling(nt,dt,ft)
        stg = Sampling(ntg,dt,ft)
        #suu = Sampling(nuu,dtgu*dt,0.0)
        #sgu = Sampling(ngu,dtgu*dt,0.0)
        #ssug = Sampling(nuu,dtgu*dt,0.0)
        ut = mul(u,dt)
        amax = [4.0,4.0,4.0,4.0,max(du)]
        plotSequences(st,[f[0],csbg[0],hsg[0],sg[0],du[0]],amax=amax,labels=["f","CSBg","HSg","Sg","du"])
        plotSequences(stg,[g[0]],amax=amax,labels=["g"])

def goSinopec1():
  #get sino images
  x0 = 250
  f,g,u = getSinoTrace(x0)
  ng = len(g)
  nu = len(u)

  #Wavelet estimation parameters
  nb,kb = 2,0#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 2,0# sampling for wavelet H 
  nd,kd = nc,kc# sampling for wavelet H 

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  niter = 200
  sfac = 0.001
  #sfac = 0.1
  warp = Warper()
  ww = WaveletWarpingCBGNOld()
  ww.setParallel(True)
  ww.setMaxPercentChange(0.01)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  ww.setLineSearchMinScale(0.0000)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
  one = zerofloat(nb)
  one[-kb] = 1.0
  hw = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  vv = ww.getVV()
  stepl = ww.getStepLength()
  deltamagb = ww.getDeltaMagB()
  deltamagc = ww.getDeltaMagC()
  condnum = ww.getCondNum()
  twonormb = ww.getTwoNormb()
  twonormc = ww.getTwoNormc()

  #Normalize wavelets
  nhw = normalizeM(hw)
  ncw = normalizeM(cw)
  nbw = normalizeM(bw)
  nbguess = normalizeM(b)
  ndw = normalizeM(dw)

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


  #rms check
  """
  hsgmf = sub(hsg,f)
  sgmf = sub(sg,f)
  rmsf = ww.rms(f)
  rmssg = ww.rms(sg)
  rmshsg = ww.rms(hsg)
  rmshsgmf = ww.rms(hsgmf)
  rmssgmf = ww.rms(sgmf)
  print "rmsf = "+str(rmsf)
  print "rmssg = "+str(rmssg)
  print "rmshsg = "+str(rmshsg)
  print "rmshsgmf = "+str(rmshsgmf)
  print "rmssgmf = "+str(rmssgmf)
  sp1 = SimplePlot()
  sp2 = SimplePlot()
  sp3 = SimplePlot()
  sp4 = SimplePlot()
  sp5 = SimplePlot()
  sp1.setTitle("f")
  sp2.setTitle("sg")
  sp3.setTitle("hsg")
  sp4.setTitle("hsg-f")
  sp5.setTitle("sg-f")
  sp1.setHLimits(tmin,tmax)
  sp2.setHLimits(tmin,tmax)
  sp3.setHLimits(tmin,tmax)
  sp4.setHLimits(tmin,tmax)
  sp5.setHLimits(tmin,tmax)
  sp1.addPoints(f)
  sp2.addPoints(sg)
  sp3.addPoints(hsg)
  sp4.addPoints(hsgmf)
  sp5.addPoints(sgmf)
  """

  #plotting
  pngDir = None
  #pngDir = "./wavelets/"

  dt = 0.004
  ft = 0.000
  nt = len(f)
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  plotUnknownWavelets(Sampling(nc,dt,kc*dt),[ncw],title="Estimated c iteration=200",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(nb,dt,kb*dt),[nbw],title="Estimated b iteration=200",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(nb,dt,kb*dt),[nbguess],title="1st Guess b (unit impulse)",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(nd,dt,kd*dt),[nhw],title="1st Guess c (shaping filter)",pngDir=pngDir,onecol=True)
  du = computeBackDiff(u)

  #getWarpingPartsSpectrums(4.0,tmin,tmax,f,g,u)
  #plotAmplitudeSpectrumT(Sampling(nb,1.0,kb*1.0), nbw, 0, nb, "amp b", amax=None)
  #plotAmplitudeSpectrumT(Sampling(nc,1.0,kc*1.0), ncw, 0, nc, "amp c", amax=None)

  if len(rfrf)>0:
    maxriri = max(riri)
    title = "rms ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    lastiter = ww.getLastIter()
    maxriri = max(riri)
    title = "rms ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    title = "Meas1"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,[stepl,vv],\
    color=None,\
    vlabel=["StepLength","2NormGradSq"],\
    vminmax=[None,None],\
    vint=[None,None],\
    hlabel="Iterations",hminmax=[0,lastiter],hint=None,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)

    title = "Meas2"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,[deltamagb,deltamagc,twonormb,twonormc,condnum],\
    color=None,\
    vlabel=["DeltMagb","DeltaMagc","2Normb","2Normc","CondNum"],\
    vminmax=[None,None,None,None,None],\
    vint=[None,None,None,None,None],\
    hlabel="Iterations",hminmax=[0,lastiter],hint=None,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)
    
  
  maxf = max(f)
  maxfd2 = maxf/2.0
  ntg = len(g)

  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  stg = Sampling(ntg,dt,ft)
  #suu = Sampling(nuu,dtgu*dt,0.0)
  #sgu = Sampling(ngu,dtgu*dt,0.0)
  #ssug = Sampling(nuu,dtgu*dt,0.0)
  su = Sampling(ntg,dt,ft)
  ut = mul(u,dt)

  maxsg = max(sg)
  amax = [maxsg,maxsg,maxsg,maxsg,max(du)]
  #plotSequences(st,[f,csbg,hsg,sg,du],amax=amax,labels=["f","CSBg","HSg","Sg","du"])
  #plotSequences(stg,[g],amax=amax,labels=["g"])

def goSinopec():
  #get sino images
  x0,nx = 250,100
  f,g,u = getSinoImage(x0,nx)
  nx = len(f)
  ng = len(g[0])
  nu = len(u[0])

  #Wavelet estimation parameters
  nb,kb = 3,-1#sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 21,-10#sampling for wavelet H 
  nd,kd = 21,-10# sampling for wavelet H 

  #set tmin and tmax 
  tmin = 100
  tmax = 500

  #Estimate wavelet
  niter = 1000
  sfac = 0.0001
  maxpc = 0.000000
  warp = Warper()
  ww = WaveletWarpingCBGN()
  ww.setParallel(False)
  ww.setMaxPercentChange(maxpc)#units are percentage.
  ww.setTimeRange(tmin,tmax)
  ww.setStabilityFactor(sfac)
  ww.setLineSearchMinScale(0.0)
  #First guesses of c and b. 
  b = zerofloat(nb)
  b[-kb] = 1.0
  c = ww.getWaveletH(nc,kc,nb,kb,b,u,f,g)
  cbw = ww.getWaveletCInverseB(nb,kb,b,nc,kc,c,u,f,g,niter)
  cw = cbw[0]
  bw = cbw[1]
  dw = ww.getWaveletH(nb,kb,bw,nc,kc)
  one = zerofloat(nb)
  one[-kb] = 1
  hw = ww.getWaveletH(nc,kc,nb,kb,one,u,f,g)
  riri = ww.getRiRi()
  rfrf = ww.getRfRf()
  vv = ww.getVV()
  stepl = ww.getStepLength()
  deltamagb = ww.getDeltaMagB()
  deltamagc = ww.getDeltaMagC()
  condnum = ww.getCondNum()
  twonormb = ww.getTwoNormb()
  twonormc = ww.getTwoNormc()

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
  #pngDir = "./png/sino/"

  
  dt = 0.004
  ft = 0.000
  nt = len(f[0])
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  plotUnknownWavelets(Sampling(nc,dt,kc*dt),[ncw],title="c",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(nb,dt,kb*dt),[nbw],title="b",pngDir=pngDir,onecol=True)
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
  plotImage(st,sx,f,tmin=0.4,tmax=1.6,title="f")
  plotImage(su,sx,g,tmin=0.9,tmax=2.8,title="smallg")
  plotImage(su,sx,g,title="smallg")
  plotImage(su,sx,bg,tmin=0.4,tmax=1.6,title="Bg")
  plotImage(st,sx,sbg,tmin=0.4,tmax=1.6,title="SBg")
  plotImage(st,sx,csbg,tmin=0.4,tmax=1.6,title="CSBg")
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
    title = "rms ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    lastiter = ww.getLastIter()
    maxriri = max(riri)
    title = "rms ri (black) rms rf (blue)"
    si = Sampling(niter,1.0,0.0)
    plotSequenceMeas(si,[riri,rfrf],amax=maxriri,hlabel="Iterations",labels="Sum of square diff",pngDir=None,title=title)

    title = "Meas1"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,[stepl,vv],\
    color=None,\
    vlabel=["StepLength","2NormGradSq"],\
    vminmax=[None,None],\
    vint=[None,None],\
    hlabel="Iterations",hminmax=[0,lastiter],hint=None,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)

    title = "Meas2"
    siter = Sampling(niter,1.0,0.0)
    plotting.plotMeasOnTopOfEachOther(siter,[deltamagb,deltamagc,twonormb,twonormc,condnum],\
    color=None,\
    vlabel=["DeltMagb","DeltaMagc","2Normb","2Normc","CondNum"],\
    vminmax=[None,None,None,None,None],\
    vint=[None,None,None,None,None],\
    hlabel="Iterations",hminmax=[0,lastiter],hint=None,\
    title=title,pngDir=pngDir,\
    slide=None,fracWidth=None,fracHeight=None,\
    paper=True,onecol=True,twocol=False)

def goSimpleTest():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  freq,decay,= 0.08,0.05
  na,ka = 9,0 #sampling for inverse wavelet A 
  nc,kc = 181,-90 # sampling for wavelet H 
  mp = False#is wavelet in f mininmum phase?
  c = 0.0#The amount of shift between p and q.
  r0,r1 = 2.0,1.00#The first and last warping amounts that will be applied to trace p. 
                 #A line is used to calculate the warping amount for a particular time.
  noise = False
  nrmsf = 0.01
  nrmsg = 0.01

  #akf = getInverseWavelet(freqf,decayf,ncf,kcf,naf,kaf,mpf)
  #akg = getInverseWavelet(freqg,decayg,ncg,kcg,nag,kag,mpg)
  #Estimation and sampling parameters.
  sfac = 1.000
  p,q,f,g,u,tmin,tmax = createSyntheticLn1D(freq,decay,mp,r0,r1,c,noise,nrmsf,nrmsg,\
  nt,ni,randomi,moreps)
  #p,q,f,g,u,tmin,tmax = createSyntheticCos(freq,decay,mp,r0,r1,c,noise,nrmsf,nrmsg,nt,ni)

  
  #main piece of wavelet estimation code
  ww = WaveletWarpingAEig()
  ww.setTimeRange(tmin,tmax)
  #ww.setGaussWeights(1.0,0.0,50,nt-1,nt)
  bw = ww.getInverseA(na,ka,u,f,g)
  w = ww.getWeights1D()
  print "2 norm squared = "+str(norm2sq(bw))
  #a = createSingleArray(naf+nag,na,akf,na,akg)
  #aunit = createUnitVector(a)
  #print "2 norm squared = "+str(norm2sq(aunit))

  #dump(aunit)
  
  hw = ww.getWaveletH(na,ka,bw,nc,kc)
  hk = getWavelet(freq,decay,nc,kc,mp) # known wavelet
  ak = ww.getWaveletH(nc,kc,hk,na,ka)

  #Eigenvalue and eigenvector analysis
  #for i in range(1):
    #print "eigenvector "+str(i)
    #dump(ww.getEigVector(i))

  print "(lambda1-lambda0)/lambdan/eps = "+str(ww.getUniqMeasure1())
  print "eigenvalues"
  #dump(ww.getEigVals())

  
  #Sum of squared differences between Fa and SGa
  print "Sum of squared differences between Fa and SGa = "+\
  str(ww.getSumSqDiff(na,ka,u,f,g))


  du = computeBackDiff(u)

  #QC wavelet estimation and warping
  sq = ww.applyS(u,q)
  slq = ww.applyW(u,q)

  slg = ww.applyW(u,g)

  af = ww.applyH(na,ka,bw,f)
  ag = ww.applyH(na,ka,bw,g)
  slag = ww.applyW(u,ag) 
  hslag = ww.applyH(nc,kc,hw,slag)
  fsubhslag = sub(f,hslag)

  akf = ww.applyH(na,ka,ak,f)
  akg = ww.applyH(na,ka,ak,g)
  slakg = ww.applyW(u,akg) 
  hslakg = ww.applyH(nc,kc,hk,slakg)

  amax = 5
  
  #Normalize wavelets
  ncw = normalizeMAAWOS(hw)
  nck = normalizeMAAWOS(hk)
  nbw = normalizeMAAWOS(bw)
  nak = normalizeMAAWOS(ak)
  print "***********************************"
  print "nak"
  #dump(nak)
  print "nbw = "
  #dump(nbw)
  print "***********************************"
  print "naw"
  #dump(naw)
  print "nak"
  #dump(nak)
  print "tmax = "+str(tmax)

  #plotting
  dt,ft = 0.004,0.000
  st = Sampling(nt,dt,ft)
  title = "na="+str(na)+" ka="+str(ka)+" nc="+str(nc)+" kc="+str(kc)+\
  " r0="+str(r0)+"r1="+str(r1)+"c="+str(c)
  amax = [max(f),max(f),nt*dt,max(du),max(slg),max(f),max(f)]
  tmark = [1,1,1.0,max(du),max(slg),1,1]
  plotSequences(st,[f,g,u,du,slg,sub(f,slg),hslag],labels=["f","g","u","du/dt","Sg","f-Sg","HSAg","f-HSAg"],pngDir=pngDir,title="fgslgsubhslagsub"+title)
  plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title="a"+title,pngDir=pngDir,onecol=True)
  plotSequences(st,[p,q,slq,sub(p,slq)],labels=["p","q","slq","p-slq"],pngDir=pngDir,title="fgslgsubhslagsub"+title)

def goBuild2DSynthetics():
#Synthetic parameters
  nt,nx,ni = 481,800,55# number of time samples in p and q; number of random impulses in p and q.
  v = 0.0#The amount of shift between p and q.
  r0,r1 = 2.0,1.0#constant u'(t)
  freq,decay,= 0.08,0.05
  mp = False#is wavelet in f mininmum phase?
  nrmsf = 0.0
  nrmsg = 0.0

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = createSyntheticLn2D(freq,decay,mp,\
  r0,r1,v,nrmsf,nrmsg,nx,nt,ni,randomi,moreps)
  ww = WaveletWarpingAEig()
  dlsuq = ww.applyW(u,q)
  dlsug = ww.applyW(u,g)

  #plotting 
  st = Sampling(nt,0.004,0.0)
  sx = Sampling(nx,0.025,0.0)
  maxf = max(f)
  maxfd2 = maxf/2.0
  amax = [maxf,maxf,maxf,maxf,maxf,maxf,r0]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,maxfd2,0.5]
  title = "2D test"
  plotSequences(st,[p[0],q[0],dlsuq[0],f[0],g[0],dlsug[0]],amax=amax,tmark=tmark,\
  labels=["p","q","dlsuq","f","g","dlsug"],title=title)

def goNaVarySyntheticTest():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  freq,decay,= 0.08,0.05
  nc,kc = 181,-90 # sampling for wavelet H 
  mp = True#is wavelet in f minimum phase?
  r0,r1,c = 2.0,1.0,0.0
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

  #measures = zerofloat(nc,nr1,nr0)
  measures = zerofloat(nna)
  ka = 0
  for na in range(na0,na1,dna):
    print "na = "+str(na)
    #Estimation and sampling parameters.
    sfac = 1.000
    p,q,f,g,u,tmin,tmax = createSyntheticLn1D(freq,decay,mp,r0,r1,c,\
    noise,nrmsf,nrmsg,nt,ni,randomi,moreps)
    du = computeBackDiff(u)
    nt = len(f)

    #main piece of wavelet estimation code
    ww = WaveletWarpingAEig()
    ww.setTimeRange(tmin,tmax)
    aw = ww.getInverseA(na,ka,u,f,g)
    hw = ww.getWaveletH(na,ka,aw,nc,kc)
    hk = getWavelet(freq,decay,nc,kc,mp) # known wavelet
    ak = getMinPhaseInverseWavelet(na,freq,decay)
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
  folder = "natestsr02r110"
  title = "r0"+str(r0)+"r1"+str(r1)+"na0"+str(na0)+"na1"+str(na1)+"ka"+str(ka)+"ni"+str(ni)
  aos1 = ArrayOutputStream("./"+folder+"/"+title+"measures1")
  aos1.writeFloats(measures1)
  aos1.close()
  aos2 = ArrayOutputStream("./"+folder+"/"+title+"measures2")
  aos2.writeFloats(measures2)
  aos2.close()
  aos3 = ArrayOutputStream("./"+folder+"/"+title+"naks")
  aos3.writeFloats(naks)
  aos3.close()
  aos4 = ArrayOutputStream("./"+folder+"/"+title+"naws")
  aos4.writeFloats(naws)
  aos4.close()
  aos5 = ArrayOutputStream("./"+folder+"/"+title+"l0ns")
  aos5.writeFloats(l0ns)
  aos5.close()
  aos6 = ArrayOutputStream("./"+folder+"/"+title+"l1ns")
  aos6.writeFloats(l1ns)
  aos6.close()
  aos7 = ArrayOutputStream("./"+folder+"/"+title+"ncks")
  aos7.writeFloats(ncks)
  aos7.close()
  aos8 = ArrayOutputStream("./"+folder+"/"+title+"ncws")
  aos8.writeFloats(ncws)
  aos8.close()

def plotNaVarySyntheticTest():
  #Synthetic parameters

  title = "r02.0r11.0na03na1480ka0ni54"
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  freq,decay,= 0.08,0.05
  nc,kc = 181,-90 # sampling for wavelet H 
  mp = True#is wavelet in f minimum phase?
  r0,r1,c = 2.0,1.0,0.0
  noise = False
  nrmsf = 0.01
  nrmsg = 0.01
  na0 = 3
  na1 = 480
  dna = 10
  ka = 0

  #title = "r02.0r11.99999na03na1480ka0ni54"
  #nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  #freq,decay,= 0.08,0.05
  #nc,kc = 181,-90 # sampling for wavelet H 
  #mp = True#is wavelet in f minimum phase?
  #r0,r1,c = 2.0,1.99999,0.0
  #noise = False
  #nrmsf = 0.01
  #nrmsg = 0.01
  #na0 = 3
  #na1 = 480
  #dna = 10
  #ka = 0

  #title = "r02.0r11.99999na03na1960ka0ni108"
  #nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  #freq,decay,= 0.08,0.05
  #nc,kc = 181,-90 # sampling for wavelet H 
  #mp = True#is wavelet in f minimum phase?
  #r0,r1,c = 2.0,1.99999,0.0
  #noise = False
  #nrmsf = 0.01
  #nrmsg = 0.01
  #na0 = 3
  #na1 = (480*2)
  #dna = (10*2)
  #ka = 0

  p,q,f,g,u,tmin,tmax = createSyntheticLn1D(freq,decay,mp,r0,r1,c,\
  noise,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)

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
  print str(tmax)

  dt = 0.004
  #Read and plot
  #title = "r0"+str(r0)+"r1"+str(r1)+"na0"+str(na0)+"na1"+str(na1)+"ka"+str(ka)+"ni"+str(ni)
  folder = "natestsr02r110"
  ais1 = ArrayInputStream("./"+folder+"/"+title+"measures1")
  ais1.readFloats(measures1)
  ais1.close()
  ais2 = ArrayInputStream("./"+folder+"/"+title+"measures2")
  ais2.readFloats(measures2)
  ais2.close()
  ais3 = ArrayInputStream("./"+folder+"/"+title+"naks")
  ais3.readFloats(naks)
  ais3.close()
  ais4 = ArrayInputStream("./"+folder+"/"+title+"naws")
  ais4.readFloats(naws)
  ais4.close()
  ais5 = ArrayInputStream("./"+folder+"/"+title+"l0ns")
  ais5.readFloats(l0ns)
  ais5.close()
  ais6 = ArrayInputStream("./"+folder+"/"+title+"l1ns")
  ais6.readFloats(l1ns)
  ais6.close()
  ais7 = ArrayInputStream("./"+folder+"/"+title+"ncks")
  ais7.readFloats(ncks)
  ais7.close()
  ais8 = ArrayInputStream("./"+folder+"/"+title+"ncws")
  ais8.readFloats(ncws)
  ais8.close()


  snaks = zerofloat(nna,maxna)
  snaws = zerofloat(nna,maxna)
  for i in range(0,nna):
    for j in range(0,maxna):
      snaks[j][i] = naks[i][j]
      snaws[j][i] = naws[i][j]
  #for i in range(0,nna):
  #  plotWavelets(Sampling(nc,1.0,kc*dt),[ncws[i],ncks[i]],title="na = "+str(i*dna+na0)+title,pngDir=pngDir,onecol=True)
  #for i in range(0,13):
  #  plotWavelets(Sampling(nna,dna,na0),[snaws[i],snaks[i]],title="naws"+str(i)+" and naks"+str(i)+title,pngDir=pngDir,onecol=True)
  plotSequenceMeas(Sampling((nna),dna,na0),[measures1],amax=[10.0],tmark=[5.0],labels=["m"],pngDir=pngDir,title=title+"zoommeasure1")
  plotSequenceMeas(Sampling((nna),dna,na0),[measures2],amax=[10.0],tmark=[5.0],labels=["m2"],pngDir=pngDir,title=title+"zoommeasure2")
  plotSequenceMeas(Sampling((nna),dna,na0),[measures1],labels=["m"],pngDir=pngDir,title=title+"measure1")
  plotSequenceMeas(Sampling((nna),dna,na0),[measures2],labels=["m2"],pngDir=pngDir,title=title+"measure2")
  plotSequenceMeas(Sampling((nna),dna,na0),[l0ns],labels=["l0ns"],pngDir=pngDir,title=title+"l0ns")
  plotSequenceMeas(Sampling((nna),dna,na0),[l1ns],labels=["l1ns"],pngDir=pngDir,title=title+"l1ns")
  st = Sampling(nt,0.004,0.00)
  ut = mul(u,0.004)
  amax = [3.2,3.2,3.2,3.2,nt*0.004,r0]
  tmark = [3.2/2.0,3.2/2.0,3.2/2.0,3.2/2.0,nt*0.004/2.0,r0/2.0]
  plotSequences(st,[p,q,f,g,ut,du],amax=amax,tmark=tmark,labels=["p","q","f","g","u","du"],pngDir=pngDir,title="pqfgudu"+title)



def goVarySyntheticTest():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
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
            hw = ww.getWaveletH(na,ka,aw,nc,kc)
            hk = getWavelet(freq,decay,nc,kc,mp) # known wavelet
            #ak = ww.getWaveletH(nc,kc,hk,na,ka)
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
def goAcVsAa():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
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

#Uses a range of r0, r1, and c to test different u(t) expressions
def goTestU():
  nt = 481
  ni = 54
  r00,r01,dr0 = 2.0,2.0,0.2
  nr0 = int((r01-r00)/dr0)+1
  r10,r11,dr1 = 1.0,2.0,0.2
  nr1 = int((r11-r10)/dr1)+1
  c0,c1,dc = 200.0,200.0,10.0
  nc = int((c1-c0)/dc)+1

  print "nr0 = "+str(nr0)
  print "nr1 = "+str(nr1)
  print "nc = "+str(nc)
  umax = nt-1 
  for ir0 in range(nr0):
    for ir1 in range(nr1):
      for ic in range(nc):
        r0 = r00+dr0*ir0  
        r1 = r10+dr1*ir1  
        c = c0+dc*ic
        #u,p,q,itmin,itmax = cosupq(r0,r1,c,nt,ni)
        #u,p,q,itmin,itmax = lnupq(r0,r1,c,nt,ni)
        u,p,q,itmin,itmax = sqrtupq(r0,r1,c,nt,ni)
        du = computeBackDiff(u)

        #plotting
        nt,dt,ft = nt,0.004,0.000
        st = Sampling(nt,dt,ft)
        title = "r0 = "+str(r0)+" r1 = "+str(r1)+" c = "+str(c)
        amax = [umax*dt,max(du)]
        tmark = [1.0,1.0]
        print "itmax = "+str(itmax)
        plotSequences(st,[mul(u,dt),du],labels=["u","du"],amax=amax,tmark=tmark,title=title,pngDir=None)
        #plotSequences(st,[u,duf],labels=["u","duf"],amax=amax,tmark=tmark,title=title,pngDir=None)

def goTestTMaxTheory():
  #Synthetic parameters
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  freq,decay,= 0.08,0.05
  nc,kc = 181,-90 # sampling for wavelet H 
  mp = True#is wavelet in f minimum phase?
  r0,r1,c = 2.0,1.0,0.0
  noise = False
  nrmsf = 0.01
  nrmsg = 0.01

  #Estimation and sampling parameters.
  sfac = 1.000
  p,q,f,g,u,tmin,tmax = createSyntheticLn1D(freq,decay,mp,r0,r1,c,\
  noise,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)
  nt = len(f)

  na0 = 3
  na1 = 10
  dna = 1
  nna = (na1-na0)/dna+1
  maxna = ((na1-na0)/dna)*dna+na0-1
  z = zerofloat(maxna,maxna,nna)
  measures = zerofloat(nna)
  ncks = zerofloat(nc,nna)
  ncws = zerofloat(nc,nna)
  lambdas = zerofloat(maxna,nna)
  print str(maxna)
  print str(nna)


  ka = 0
  for na in range(na0,na1,dna):
    print "na = "+str(na)

    #main piece of wavelet estimation code
    ww = WaveletWarpingAEig()
    ww.setTimeRange(tmin,tmax)
    aw = ww.getInverseA(na,ka,u,f,g)
    hw = ww.getWaveletH(na,ka,aw,nc,kc)
    hk = getWavelet(freq,decay,nc,kc,mp) # known wavelet
    ak = getMinPhaseInverseWavelet(na,freq,decay)
    ncw = normalizeMAAWOS(hw)
    nck = normalizeMAAWOS(hk)
    measures[(na-na0)/dna] = ww.getMeasure()
    print "meaure = "+str(ww.getMeasure())
    eigvals = ww.getEigVals() 
    print "length of eigvals = "+str(len(eigvals))
    for i in range(0,na):
      lambdas[(na-na0)/dna][i] = eigvals[i]
      print "eigvals = "+str(eigvals[i])
    ww.getZ(z[(na-na0)/dna])
    for i in range(0,nc):
      ncks[(na-na0)/dna][i] = nck[i]
      ncws[(na-na0)/dna][i] = ncw[i]


  #Write out information to file for later
  folder = "tmaxtheorytests"
  title = "r0"+str(r0)+"r1"+str(r1)+"na0"+str(na0)+"na1"+str(na1)+"ka"+str(ka)+"ni"+str(ni)
  aos1 = ArrayOutputStream("./"+folder+"/"+title+"measures")
  aos1.writeFloats(measures)
  aos1.close()
  aos2 = ArrayOutputStream("./"+folder+"/"+title+"lambdas")
  aos2.writeFloats(lambdas)
  aos2.close()
  aos3 = ArrayOutputStream("./"+folder+"/"+title+"ncks")
  aos3.writeFloats(ncks)
  aos3.close()
  aos4 = ArrayOutputStream("./"+folder+"/"+title+"ncws")
  aos4.writeFloats(ncws)
  aos4.close()
  aos5 = ArrayOutputStream("./"+folder+"/"+title+"z")
  aos5.writeFloats(z)
  aos5.close()

def plotTestTMaxTheory():
  title = "r02.0r11.0na0341na1348ka0ni54"
  nt,ni,randomi = 481,2,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = False
  freq,decay,= 0.08,0.05
  nc,kc = 181,-90 # sampling for wavelet H 
  mp = True#is wavelet in f minimum phase?
  r0,r1,c = 2.0,1.0,0.0
  noise = False
  nrmsf = 0.01
  nrmsg = 0.01
  p,q,f,g,u,tmin,tmax = createSyntheticLn1D(freq,decay,mp,r0,r1,c,\
  noise,nrmsf,nrmsg,nt,ni,randomi,moreps)
  du = computeBackDiff(u)
  SimplePlot.asPoints(f)
  SimplePlot.asPoints(g)

  na0 = 341
  na1 = 348
  dna = 1
  ka = 0




  nna = (na1-na0)/dna+1
  maxna = ((na1-na0)/dna)*dna+na0
  z = zerofloat(maxna,maxna,nna)
  measures = zerofloat(nna)
  ncks = zerofloat(nc,nna)
  ncws = zerofloat(nc,nna)
  lambdas = zerofloat(maxna,nna)

  print "maxna = "+str(maxna)
  print "nna = "+str(nna)
  print "tmax = "+str(tmax)

  dt = 0.004
  #Read and plot
  folder = "tmaxtheorytests"
  ais1 = ArrayInputStream("./"+folder+"/"+title+"measures")
  ais1.readFloats(measures)
  ais1.close()
  ais2 = ArrayInputStream("./"+folder+"/"+title+"z")
  ais2.readFloats(z)
  ais2.close()
  ais3 = ArrayInputStream("./"+folder+"/"+title+"lambdas")
  ais3.readFloats(lambdas)
  ais3.close()
  ais4 = ArrayInputStream("./"+folder+"/"+title+"ncks")
  ais4.readFloats(ncks)
  ais4.close()
  ais5 = ArrayInputStream("./"+folder+"/"+title+"ncws")
  ais5.readFloats(ncws)
  ais5.close()

  for i in range(0,nna-1):
    na = i*dna+na0
    lambdan = lambdas[i][na-1]
    plotWavelets(Sampling(nc,1.0,kc*dt),[ncws[i],ncks[i]],title="na = "+str(na)+title,pngDir=pngDir,onecol=True)
    
    plotMatrix(z[i],title=str(na)+"by"+str(na))
    plotSequenceMeas(Sampling(maxna,1.0,0.0),[div(lambdas[i],lambdan)],labels=["lambdas scaled by lambdans na= "+str(na)],pngDir=pngDir,title=title)
  plotSequenceMeas(Sampling((nna),dna,na0),[measures],labels=["Measure na = "+str(na)],pngDir=pngDir,title=title)





def goTestUniqMeas():
  WaveletWarpingAEig.ulp1()
  



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
    a = ww.getWaveletH(nc,kc,hk,na,ka) # known inverse wavelet
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
