#############################################################################
# Demo of 2 wavelet estimations from warping.

from imports import *

from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.sgl import *
from edu.mines.jtk.awt.ColorMap import *
from edu.mines.jtk.lapack import *
from wwarp import WaveletWarping, WaveletWarpingABEig, Warper
from wwarp import ShapingFilter
from java.util import Random

############################################################################

write = False
#write = True

def main(args):
  goEstimateWaveletTest()
  #goEstimateWaveletTest2D()
  #goSinopec1()
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
  #printSumSqDiff()
  #testSLG(nas[i],kas[i])
  #testAddNoise()
  #testWarpingandBandPassFilter()

def goEstimateWaveletTest():
  #Synthetic parameters
  nt,ni,randomi = 481,55,True# number of time samples in p and q; number of random impulses in p and q.
  moreps = True 
  na,ka = 81,-40 #sampling for inverse wavelet a #Note ka<=kc 
  nb,kb = 81,-40 #sampling for inverse wavelet b #Note ka<=kc 
  nc,kc = 181,-90 # sampling for wavelet c
  nd,kd = 181,-90 # sampling for wavelet d 
  v = [0.0]#The amount of shift between p and q.
  r0,r1 = [5.0],[1.0]#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  freqc,decayc,= 0.08,0.05
  mpc = True#is wavelet in f mininmum phase?
  finitec = False#does the wavelet in f have a finite length?
  freqd,decayd,= 0.08,0.05
  mpd = True#is wavelet in f mininmum phase?
  finited = False#does the wavelet in f have a finite length?
  nrmsf = [0.00]
  nrmsg = nrmsf
  nab = na+nb

  nv = len(v)
  nr0 = len(r0)
  nr1 = len(r1)
  nnoise = len(nrmsf)
  for ir0 in range(0,nr0):
    for ir1 in range(0,nr1):
      for iv in range(0,nv):
        for inoise in range(0,nnoise):
          #Create synthetic f and g.
          p,q,f,g,noisef,noiseg,u,tmin,tmax = createSyntheticLn1D(freqc,decayc,mpc,finitec,\
          freqd,decayd,mpd,finited,r0[ir0],r1[ir1],v[iv],nrmsf[inoise],nrmsg[inoise],\
          nt,ni,randomi,moreps)
          du = computeBackDiff(u)
          if (abs(ka) > abs(kb)):
            tmin = tmin+-ka
            tmax = tmax+ka
          else:
            tmin = tmin+-kb
            tmax = tmax+kb

          print "tmin = "+str(tmin)
          print "tmax = "+str(tmax)

          #Estimate wavelet
          warp = Warper()
          ww = WaveletWarpingABEig()
          ww.setTimeRange(tmin,tmax)
          abw = ww.getInverseAB(na,ka,nb,kb,u,f,g)
          print "length of abw = "+str(len(abw))
          aw = ww.getA(na,abw)
          bw = ww.getB(nb,abw)
          cw = ww.getWaveletH(na,ka,aw,nc,kc)
          dw = ww.getWaveletH(nb,kb,bw,nd,kd)
          print "aw"
          dump(aw)
          print "bw"
          dump(bw)

          #Get eigenvalue and eigenvector
          #dump(ww.getEigVector(0))
          #dump(ww.getEigVals())
          #print ww.printZ()
          #z = ww.getZ()
          #SimplePlot.asPixels(z)

          #Get known wavelet
          ck = getWavelet(freqc,decayc,nc,kc,mpc,finitec)
          dk = getWavelet(freqd,decayd,nd,kd,mpd,finited)
          ak = ww.getWaveletH(nc,kc,ck,na,ka)
          bk = ww.getWaveletH(nd,kd,dk,nb,kb)

          #Normalize wavelets
          ncw = normalizeMAAWOS(cw)
          ndw = normalizeMAAWOS(dw)
          naw = normalizeMAAWOS(aw)
          nbw = normalizeMAAWOS(bw)
          nck = normalizeMAAWOS(ck)
          ndk = normalizeMAAWOS(dk)
          nak = normalizeMAAWOS(ak)
          nbk = normalizeMAAWOS(bk)

          #Create simply warped g
          dlsug = warp.applyS(u,g)
          dlsuq = warp.applyS(u,q)
          
          #Create warped g
          bg = ww.applyA(nb,kb,bw,g)
          dlsubg = warp.applyS(u,bg)
          cdlsubg = ww.applyH(nc,kc,cw,dlsubg)
          print "(lambda1-lambda0)/(lambdan*eps) = "+str(ww.getWDMeasure())
          print "The sum of the square differences is "+str(ww.getSumSqDiff(nab))+", which "+"is equivalent to the smallest eigenvalue"+str(ww.getEigVal(0))
          print "The sum of square differences if an impulse is "+\
          "assumed to the be the inverse wavelet Check!!!"+str(ww.getSumSqDiffNoWaveletEst(na,ka,kb))
          getWarpingPartsSpectrums(r0[ir0],f,g,u)
 
          #plotting
          pngDir = None
          if (write):
            pngDir = "./png/estimatewavelet/"

          dt,ft = 0.004,0.000
          st = Sampling(nt,dt,ft)
          ut = mul(u,dt)
          #title = "nt="+str(nt)+"ni="+str(ni)+\
          #"r0="+str(r0[ir0])+"r1="+str(r1[ir1])+"v="+str(v[iv])+\
          #"nrmsf="+str(nrmsf[inoise])+"nrmsg="+str(nrmsg[inoise])
          title1 = "c"
          title2 = "d"
          title3 = "a"
          title4 = "b"
          plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title1,pngDir=pngDir,onecol=True)
          print "hi"
          plotWavelets(Sampling(nd,dt,kd*dt),[ndw,ndk],title=title2,pngDir=pngDir,onecol=True)
          print "hi1"
          plotWavelets(Sampling(na,dt,ka*dt),[naw,nak],title=title3,pngDir=pngDir,onecol=True)
          plotWavelets(Sampling(nb,dt,kb*dt),[nbw,nbk],title=title4,pngDir=pngDir,onecol=True)

          maxf = max(f)
          maxfd2 = maxf/2.0

          dt,ft = 0.004,0.000
          st = Sampling(nt,dt,ft)
          ut = mul(u,dt)
          title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0[ir0])+"r1="+str(r1[ir1])+"v="+str(v[iv])+\
          "nrmsf="+str(nrmsf[inoise])+"nrmsg="+str(nrmsg[inoise])

          amax = [maxf,maxf,maxf,maxf,r0[ir0]]
          tmark = [maxfd2,maxfd2,maxfd2,maxfd2,0.5]
          plotSequences(st,[f,g,dlsug,cdlsubg,du],amax=amax,tmark=tmark,\
          labels=["f","g","dlsug","cdlsubg","du"],pngDir=pngDir,\
          title=title)
          SimplePlot.asPoints(dlsuq)

def goEstimateWaveletTest2D():
  #Synthetic parameters
  nt,nx,ni = 481,30,55# number of time samples in p and q; number of random impulses in p and q.
  na,ka = 81,-40 #sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 181,-90 # sampling for wavelet H 
  v = [0.0]#The amount of shift between p and q.
  r0,r1 = [2.0],[1.0]#constant u'(t)
  freqc,decayc,= 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  finitec = False#does the wavelet in f have a finite length?
  freqd,decayd,= 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  finited = False#does the wavelet in f have a finite length?
  nrmsf = 0.00
  nrmsg = nrmsf

  nv = len(v)
  nr0 = len(r0)
  nr1 = len(r1)
  for ir0 in range(0,nr0):
    for ir1 in range(0,nr1):
      for iv in range(0,nv):
        #Create synthetic f and g.
        p,q,f,g,noisef,noiseg,u,tmin,tmax = createSyntheticLn2D(freqc,decayc,mpc,finitec,\
        freqd,decayd,mpd,finited,r0[ir0],r1[ir1],v[iv],nrmsf,nrmsg,nx,nt,ni,randomi,moreps)
        tmin = tmin+-ka
        tmax = tmax+ka
        print "tmin = "+str(tmin)
        print "tmax = "+str(tmax)

        #Estimate wavelet
        ww = WaveletWarpingAEig()
        ww.setTimeRange(tmin,tmax)
        aw = ww.getInverseA(na,ka,u,f,g)
        cw = ww.getWaveletH(na,ka,aw,nc,kc)

        #Get eigenvalue and eigenvector
        #dump(ww.getEigVector(0))
        #dump(ww.getEigVals())
        #print ww.printZ()
        #z = ww.getZ()
        #SimplePlot.asPixels(z)

        #Get known wavelet
        ck = getWavelet(freqc,decayc,nc,kc,mpc,finitec)
        dk = getWavelet(freqd,decayd,nc,kc,mpd,finited)
        ak = ww.getWaveletH(nc,kc,ck,na,ka)

        #Normalize wavelets
        ncw = normalizeMAAWOS(cw)
        nck = normalizeMAAWOS(ck)
        ndk = normalizeMAAWOS(dk)
        naw = normalizeMAAWOS(aw)
        nak = normalizeMAAWOS(ak)
        
        #Create simply warped q
        dlsuq = ww.applyW(u,q)

        #Create simply warped g
        dlsug = ww.applyW(u,g)
        
        #Create warped g
        ag = ww.applyA(na,ka,aw,g)
        dlsuag = ww.applyW(u,ag)
        hdlsuag = ww.applyH(nc,kc,cw,dlsuag)
        print "(lambda1-lambda0)/(lambdan*eps) = "+str(ww.getWDMeasure())
        print "The sum of the square differences is "+str(ww.getSumSqDiff(na))+", which "+"is equivalent to the smallest eigenvalue"+str(ww.getEigVal(0))
        print "The sum of square differences if an impulse is "+\
        "assumed to the be the inverse wavelet"+str(ww.getSumSqDiffNoWaveletEst(ka))
        print "aw"
        dump(aw);

        #plotting
        pngDir = None
        if (write):
          pngDir = "./png/estimatewavelet/"

        dt,ft = 0.004,0.000
        st = Sampling(nt,dt,ft)
        ut = mul(u,dt)
        title = "c"
        plotWavelets(Sampling(nc,dt,kc*dt),[ncw,nck],title=title,pngDir=pngDir,onecol=True)
        title = "d"
        plotWavelets(Sampling(nc,dt,kc*dt),[ncw,ndk],title=title,pngDir=pngDir,onecol=True)

        maxf = max(f)
        maxfd2 = maxf/2.0

        dt,ft = 0.004,0.000
        st = Sampling(nt,dt,ft)
        sx = Sampling(nx,0.025,0.0)
        ut = mul(u,dt)
        plotImage(st,sx,f,title="f")
        plotImage(st,sx,g,title="g")
        plotImage(st,sx,dlsug,title="DLSUg")
        plotImage(st,sx,hdlsuag,title="HDLSUAg")
        amax = [maxf,maxf,maxf,maxf,r0[ir0]]
        tmark = [maxfd2,maxfd2,maxfd2,maxfd2,0.5]
        plotSequences(st,[f[0],g[0],dlsug[0],hdlsuag[0]],amax=amax,tmark=tmark,\
        labels=["f","g","DLSUg","HDLSUAg"],pngDir=pngDir,title=title)
        """
        title = "nt="+str(nt)+"ni="+str(ni)+"r0="+str(r0[ir0])+"r1="+str(r1[ir1])+"v="+str(v[iv])+\
        "nrmsf="+str(nrmsf)+"nrmsg="+str(nrmsg)+\

        amax = [maxf,maxf,maxf,maxf,maxf]
        tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
        plotSequences(st,[f,g,dlsug,sub(f,dlsug)],amax=amax,tmark=tmark,\
        labels=["f","g","HWg","f-HWg"],pngDir=pngDir,\
        title=title)
        """

def goSinopec1():
  #get sino images
  x0 = 333
  f,g,u = getSinoTrace(x0)
  ng = len(g)
  nu = len(u)

  #Wavelet estimation parameters
  na,ka = 21,-10 #sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 181,-90 # sampling for wavelet H 

  #set tmin and tmax 
  tmin = 100+ka
  tmax = 400+-ka
  tmin = tmin+-ka
  tmax = tmax+ka

  #Estimate wavelet
  sfac = 1.0
  ww = WaveletWarpingABEig()
  ww.setStabilityFactor(sfac)
  ww.setTimeRange(tmin,tmax)
  aw = ww.getInverseA(na,ka,u,f,g)
  cw = ww.getWaveletH(na,ka,aw,nc,kc)

  #Get eigenvalue and eigenvector
  #dump(ww.getEigVector(0))
  #dump(ww.getEigVals())
  #print ww.printZ()
  #z = ww.getZ()
  #SimplePlot.asPixels(z)

  #Normalize wavelets
  ncw = normalizeMAAWOS(cw)
  naw = normalizeMAAWOS(aw)

  #Create simply warped g
  dlsug = ww.applyW(u,g)

  #Create very simply warped g
  sg = ww.applyS(u,g)

  
  #Create warped g
  ag = ww.applyA(na,ka,aw,g)
  dlsuag = ww.applyW(u,ag)
  hdlsuag = ww.applyH(nc,kc,cw,dlsuag)
  print "(lambda1-lambda0)/(lambdan*eps) = "+str(ww.getWDMeasure())
  print "sum of the square differences = "+str(ww.getSumSqDiff(na));
  print "The sum of the square differences is "+str(ww.getSumSqDiff(na))+", which "+"is equivalent to the smallest eigenvalue"+str(ww.getEigVal(0))
  print "The sum of square differences if an impulse is "+\
  "assumed to the be the inverse wavelet"+str(ww.getSumSqDiffNoWaveletEst(ka))
  print "aw: "
  dump(aw)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)

  #The steps of creating DLSU
  #r = 4.0
  #ngu = int(r*(ng-1)+1)
  #nuu = int(r*(nu-1)+1)
  #dtgu = 1.0/r
  #uu = ww.upSampleLinear(nu,1.0,0.0,u,nuu,dtgu,0.0)
  #ug = ww.reSampleSinc(ng,1.0,0.0,g,ngu,dtgu,0.0)
  #sug = ww.applyS(dtgu,uu,ug)
  #lsug = ww.applyL(dtgu,uu,sug)
  #nlsug = len(lsug[0])
  #dlsug1 = ww.reSampleSinc(nlsug,dtgu,0.0,lsug,nu,1.0,0.0)

  #plotting
  pngDir = None
  if (write):
    pngDir = "./png/sino/"

  dt = 0.004
  ft = 0.000
  nt = len(f)
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  plotUnknownWavelets(Sampling(nc,dt,kc*dt),[ncw],title="c 1 trace",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(na,dt,ka*dt),[naw],title="a 1 trace",pngDir=pngDir,onecol=True)
  du = computeBackDiff(u)

  plotAmplitudeSpectrumT(st, dlsug, 0, nt, "amp DLSUg", amax=None)
  plotAmplitudeSpectrumT(st, f, 0, nt, "amp f", amax=None)
  getWarpingPartsSpectrums(4.0,f,g,u)

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

  amax = [maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2]
  plotSequences(st,[f,dlsug,hdlsuag],amax=amax,tmark=tmark,\
  labels=["f","DLSUg","HDLSUAg"])
  plotSequences(stg,[g],amax=amax,tmark=tmark,\
  labels=["g"])


def goSinopec():
  #get sino images
  x0,nx = 331,50
  f,g,u = getSinoImages(x0,nx)
  nx = len(f)
  ng = len(g[0])
  nu = len(u[0])

  #Wavelet estimation parameters
  na,ka = 21,-10 #sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 181,-90 # sampling for wavelet H 

  #set tmin and tmax 
  tmin = 100+ka
  tmax = 400+-ka
  tmin = tmin+-ka
  tmax = tmax+ka

  #Estimate wavelet
  sfac = 1.0
  ww = WaveletWarpingAEig()
  ww.setStabilityFactor(sfac)
  ww.setTimeRange(tmin,tmax)
  aw = ww.getInverseA(na,ka,u,f,g)
  cw = ww.getWaveletH(na,ka,aw,nc,kc)

  #Get eigenvalue and eigenvector
  #dump(ww.getEigVector(0))
  #dump(ww.getEigVals())
  #print ww.printZ()
  #z = ww.getZ()
  #SimplePlot.asPixels(z)

  #Normalize wavelets
  ncw = normalizeMAAWOS(cw)
  naw = normalizeMAAWOS(aw)

  #Create simply warped g
  dlsug = ww.applyW(u,g)

  #Create very simply warped g
  sg = ww.applyS(u,g)

  
  #Create warped g
  ag = ww.applyA(na,ka,aw,g)
  dlsuag = ww.applyW(u,ag)
  hdlsuag = ww.applyH(nc,kc,cw,dlsuag)
  print "(lambda1-lambda0)/(lambdan*eps) = "+str(ww.getWDMeasure())
  print "sum of the square differences = "+str(ww.getSumSqDiff(na));
  print "The sum of the square differences is "+str(ww.getSumSqDiff(na))+", which "+"is equivalent to the smallest eigenvalue"+str(ww.getEigVal(0))
  print "The sum of square differences if an impulse is "+\
  "assumed to the be the inverse wavelet  "+str(ww.getSumSqDiffNoWaveletEst(ka))
  print "aw: "
  dump(aw)
  print "tmin = "+str(tmin)
  print "tmax = "+str(tmax)

  #The steps of creating DLSU
  #r = 4.0
  #ngu = int(r*(ng-1)+1)
  #nuu = int(r*(nu-1)+1)
  #dtgu = 1.0/r
  #uu = ww.upSampleLinear(nu,1.0,0.0,u,nuu,dtgu,0.0)
  #ug = ww.reSampleSinc(ng,1.0,0.0,g,ngu,dtgu,0.0)
  #sug = ww.applyS(dtgu,uu,ug)
  #lsug = ww.applyL(dtgu,uu,sug)
  #nlsug = len(lsug[0])
  #dlsug1 = ww.reSampleSinc(nlsug,dtgu,0.0,lsug,nu,1.0,0.0)

  #plotting
  pngDir = None
  if (write):
    pngDir = "./png/sino/"

  
  dt = 0.004
  ft = 0.000
  nt = len(f[0])
  st = Sampling(nt,dt,ft)
  ut = mul(u,dt)
  plotUnknownWavelets(Sampling(nc,dt,kc*dt),[ncw],title="c",pngDir=pngDir,onecol=True)
  plotUnknownWavelets(Sampling(na,dt,ka*dt),[naw],title="a",pngDir=pngDir,onecol=True)
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
  plotImage(st,sx,f,title="f")
  plotImage(su,sx,g,title="g")
  plotImage(su,sx,ag,title="Ag")
  plotImage(st,sx,dlsuag,title="DLSUAg")
  plotImage(st,sx,hdlsuag,title="HDLSUAg")
  #plotImage(st,sx,sg,title="Sg")
  #plotImage(st,sx,du,title="du")
  #plotImage(st,sx,dlsug,title="DLSUg")
  #plotImage(st,sx,dlsug1,title="DLSUg1")

  #amax = [maxf,maxf,maxf,maxf]
  #tmark = [maxfd2,maxfd2,maxfd2,maxfd2]
  #plotSequences(st,[f[0],dlsug[0],sg[0]],amax=amax,tmark=tmark,\
  #labels=["f","HWg","Sg"])

def goBuild2DSynthetics():
#Synthetic parameters
  nt,nx,ni = 481,800,55# number of time samples in p and q; number of random impulses in p and q.
  v = 0.0#The amount of shift between p and q.
  r0,r1 = 2.0,1.0#constant u'(t)
  freq,decay,= 0.08,0.05
  mp = False#is wavelet in f mininmum phase?
  finite = False#does the wavelet in f have a finite length?
  nrmsf = 0.0
  nrmsg = 0.0

  #Create synthetic f and g.
  p,q,f,g,noisef,noiseg,u,tmin,tmax = createSyntheticLn2D(freq,decay,mp,finite,\
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

def createSyntheticCos(freq,decay,mp,finite,
  r0,r1,v,noise,nrmsf,nrmsg,nt,ni):
  u,p,q,itmin,itmax = cosupq(r0,r1,v,nt,ni)
  f = addWavelet(freq,decay,p,mp,finite)
  g = addWavelet(freq,decay,q,mp,finite)
  if noise:
    f = addNoise(nrmsf,42,f)
    g = addNoise(nrmsg,43,g)
  return p,q,f,g,u,itmin,itmax

def createSyntheticLn1D(freqc,decayc,mpc,finitec,freqd,decayd,mpd,finited,
  r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps):
  u,p,q,itmin,itmax = lnupq(r0,r1,v,nt,ni,randomi,moreps)
  f = addWavelet(freqc,decayc,p,mpc,finitec)
  g = addWavelet(freqd,decayd,q,mpd,finited)
  f,noisef = addNoise(nrmsf,42,f)
  g,noiseg = addNoise(nrmsg,43,g)
  return p,q,f,g,noisef,noiseg,u,itmin,itmax

def createSyntheticLn2D(freqc,decayc,mpc,finitec,freqd,decayd,mpd,finited,
  r0,r1,v,nrmsf,nrmsg,nx,nt,ni,randomi,moreps):
  u1,p1,q1,itmin,itmax = lnupq(r0,r1,v,nt,ni,randomi,moreps)
  u = replicateTrace(nx,u1)
  p = replicateTrace(nx,p1)
  q = replicateTrace(nx,q1)
  f = addWavelet2D(freqc,decayc,p,mpc,finitec)
  g = addWavelet2D(freqd,decayd,q,mpd,finited)
  f,noisef = addNoise2D(nrmsf,42,f)
  g,noiseg = addNoise2D(nrmsg,43,g)
  #dump(f)
  print "itmin = "+str(itmin)
  print "itmax = "+str(itmax)
  return p,q,f,g,noisef,noiseg,u,itmin,itmax

#Copies the same trace to create a 2D image
def replicateTrace(nx,f1):
  nt = len(f1)
  f = zerofloat(nt,nx)
  for i in range(0,nx):
    f[i] = f1
  return f
  

def zeroOnAndAfter(f,itzero):
  nf = len(f)
  for i in range(itzero,nf):
    f[i] = 0.0
  return f

def cosupq(r0,r1,c,nt,ni):
  umax = nt-1
  p = zerofloat(nt)
  q = zerofloat(nt)
  u = zerofloat(nt)
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.45)
  if (r0==r1):
    for n in range(nt):
      u[n] = r0*n+c
    dq = umax/float(ni+1)#float cast ensures float division.
    ts = rampfloat(dq,dq,ni)
    ran = Random(55)
    for ji in range(ni):
      rj = 2.0*ran.nextFloat()-1.0
      tq = ts[ji] #time u
      tp = (tq-c)/r0#time t(u)
      si.accumulate(tp,rj,nt,1.0,0.0,p)
      si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = 0
    itmax = int((umax-c)/r0+0.5)
    return u,p,q,itmin,itmax
  else:
    b = r0
    a = b*sqrt(1-(r1/r0)*(r1/r0))/(umax-c)
    for n in range(nt):
      u[n] = b*sin(a*n)/a+c
    dq = umax/float(ni+1)#float cast ensures float division.
    ts = rampfloat(dq,dq,ni)
    ran = Random(55)
    for ji in range(ni):
      rj = 2.0*ran.nextFloat()-1.0
      tq = ts[ji] #time u
      tp = asin(a*(tq-c)/b)/a
      si.accumulate(tp,rj,nt,1.0,0.0,p)
      si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = 0
    itmax = int((asin(a*(umax-c)/b)/a)+0.5)
    return u,p,q,itmin,itmax

def lnupq(r0,r1,v,nt,ni,randomi,moreps):
  umax = nt-1
  p = zerofloat(nt)
  q = zerofloat(nt)
  u = zerofloat(nt)
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.40)
  print "max length = "+str(si.getMaximumLength())
  if (r0==r1):
    print "r0==r1"
    for n in range(nt):
      u[n] = r0*n+v
    dq = (umax-v)/float(ni+1)
    ts = rampfloat(dq+v+20,dq,ni)  
    ran = Random(55)
    #rj=1.0
    for ji in range(ni):
      rj = 1.0
      if randomi:
        rj = 2.0*ran.nextFloat()-1.0
      tq = ts[ji]#time u
      tp = (tq-v)/r0#time t(u)
      #print "tq = "+str(tq)
      #print "tp = "+str(tp)
      si.accumulate(tp,rj,nt,1.0,0.0,p)
      si.accumulate(tq,rj,nt,1.0,0.0,q)
    #for impulses seen in p, but not in q
    if moreps:
      tsni = ts[ni-1]
      ulastp = r0*nt+v#the u that corresponds to the last impulse in p
      nexi = int((ulastp-tsni)/dq)#extra i to fill in the rest of p
      exts = rampfloat(tsni+dq,dq,nexi)
      for ji in range(nexi):
        rj = 1.0
        if randomi:
          rj = 2.0*ran.nextFloat()-1.0
        tq = exts[ji]#time u
        tp = (tq-v)/r0#time t(u)
        #print "tq = "+str(tq)
        #print "tp = "+str(tp)
        si.accumulate(tp,rj,nt,1.0,0.0,p)
        si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = 0
    itmax = int((umax-v)/r0+0.5)
    return u,p,q,itmin,itmax
  else: 
    print "r0>r1"
    a = (umax-v)/log(r0/r1)
    b = r0*log(r0/r1)/(umax-v)
    for n in range(nt):
      u[n] = a*log(1.0+b*n)+v
    dq = (umax)/float(ni+1)
    ts = rampfloat(dq+v+20,dq,ni)  
    ran = Random(55)
    #rj = 1.0
    for ji in range(ni):
      rj = 1.0
      if randomi:
        rj = 2.0*ran.nextFloat()-1.0
      tq = ts[ji]#time u
      tp = (exp((tq-v)/a)-1.0)/b#time t(u)
      #print "tq = "+str(tq)
      #print "tp = "+str(tp)
      si.accumulate(tp,rj,nt,1.0,0.0,p)
      si.accumulate(tq,rj,nt,1.0,0.0,q)
    #for impulses seen in p, but not in q
    if moreps:
      tsni = ts[ni-1]
      ulastp = a*log(b*nt+1)+v#the u that corresponds to the last impulse in p
      nexi = int((ulastp-tsni)/dq)#extra i to fill in the rest of p
      exts = rampfloat(tsni+dq,dq,nexi)
      for ji in range(nexi):
        rj = 1.0
        if randomi:
          rj = 2.0*ran.nextFloat()-1.0
        tq = exts[ji]#time u
        tp = (exp((tq-v)/a)-1.0)/b#time t(u)
        #print "tq = "+str(tq)
        #print "tp = "+str(tp)
        si.accumulate(tp,rj,nt,1.0,0.0,p)
        si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = 0
    itmax = int(((exp((umax-v)/a)-1.0)/b)+0.5)
    return u,p,q,itmin,itmax

def sqrtupq(r0,r1,v,nt,ni):
  umax = nt-1
  p = zerofloat(nt)
  q = zerofloat(nt)
  u = zerofloat(nt)
  si = SincInterpolator.fromErrorAndFrequency(0.01,0.45)
  if (r0==r1):
    for n in range(nt):
      u[n] = r0*n+v
    dq = umax/float(ni+1)
    ts = rampfloat(dq,dq,ni)  
    ran = Random(55)
    for ji in range(ni):
      rj = 2.0*ran.nextFloat()-1.0
      tq = ts[ji]#time u
      tp = (tq-v)/r0#time t(u)
      si.accumulate(tp,rj,nt,1.0,0.0,p)
      si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = 0
    itmax = int((tq-v)/r0+0.5)
    return u,p,q,itmin,itmax
  else:
    a = 2.0*(pow(r1,3.0)-pow(r0,3.0))/(3.0*(umax-v))
    b = r0*r0
    for n in range(nt):
      u[n] = (2.0/(3.0*a))*(pow((a*n+b),(3.0/2.0))-pow(b,(3.0/2.0)))+v
    dq = umax/float(ni+1)
    ts = rampfloat(dq,dq,ni)  
    si = SincInterpolator.fromErrorAndFrequency(0.01,0.45)
    ran = Random(55)
    for ji in range(ni):
      rj = 2.0*ran.nextFloat()-1.0
      tq = ts[ji]#time u
      tp = (1.0/a)*(pow(3.0*a*(tq-v)/2.0+pow(b,3.0/2.0),(2.0/3.0))-b)#time t(u)
      si.accumulate(tp,rj,nt,1.0,0.0,p)
      si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = 0
    itmax = int(((1.0/a)*(pow(3.0*a*(umax-v)/2.0+pow(b,3.0/2.0),(2.0/3.0))-b))+0.5)
    return u,p,q,itmin,itmax

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


def addWavelet(fpeak,decay,p,mp=False,finite=False):
  w = 2.0*PI*fpeak
  #print "w = "+str(w)
  if not mp:
    decay *= 2.0
    w -= 2.0*PI*0.04
    #print "decay = "+str(decay)
    #print "w= "+str(w)
  r = exp(-decay)
  #print "r= "+str(r)
  a1,a2 = -2.0*r*cos(w),r*r
  #print "a =",[1,a1,a2]
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  #print "pole1 = "+poles[0].toString()
  #print "pole2 = "+poles[1].toString()
  zeros = []
  #print "zero1 = "+zeros[0].toString()
  gain = 1.0
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.applyForward(p,t)
  if not mp:
    w = 2.0*PI*(fpeak+0.04)
    a1r,a2r = -2.0*r*cos(w),r*r
    #print "areverse =",[1,a1r,a2r]
    #print "wreverse= "+str(w)
    poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
    #print "pole1reverse = "+poles[0].toString()
    #print "pole2reverse = "+poles[1].toString()
    zeros = []
    gain = 1.0
    rcf = RecursiveCascadeFilter(poles,zeros,gain)
    rcf.applyReverse(t,x)
  #if finite:
  #  w = 2.0*PI*fpeak
  #  a1r,a2r = -2.0*r*cos(w),r*r
  #  print "areverse =",[1,a1r,a2r]
  #  print "wreverse= "+str(w)
  #  poles = []
  #  zeros = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  #  print "zero1reverse = "+zeros[0].toString()
  #  print "zero2reverse = "+zeros[1].toString()
  #  gain = 1.0
  #  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  #  rcf.applyReverse(t,x)

  else:
    copy(t,x)
  if not mp:
    conv(2,0,[1.0,-0.95],len(x),0,copy(x),len(x),0,x) # attenuate DC
  return x

def addWavelet2D(fpeak,decay,p,mp,finite):
  nx = len(p)
  nt = len(p[0])
  f = zerofloat(nt,nx)
  for ix in range(0,nx):
    f[ix] = addWavelet(fpeak,decay,p[ix],mp,finite)
  return f

def getSimpleInverseWavelet(fpeak,decay):
  w = 2.0*PI*fpeak
  #print "w = "+str(w)
  r = exp(-decay)
  #print "r= "+str(r)
  a1,a2 = -2.0*r*cos(w),r*r
  a = [1,a1,a2]
  return a

def getInverseWavelet(fpeak,decay,nc,kc,na,ka,mp,finite):
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
    
def getWavelet(fpeak,decay,nc,kc,mp=False,finite=False):
  x = zerofloat(nc)
  x[-kc] = 1.0
  return addWavelet(fpeak,decay,x,mp,finite)

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

def plotImages(st,sx,xs,amax=None,tmark=None,pngDir=None,labels=None,title=None):
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
      pp.setVLimits(ix,-amax[ix],amax[ix])
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
      #pv.setLineStyle(PointsView.Line.SOLID)
      pv.setLineStyle(PointsView.Line.NONE)
      #pv.setLineColor(Color.BLUE)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(3.4)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
    elif hs[ih]:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.SOLID)
      pv.setLineColor(Color.RED)
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
      pngDir = title+"onecol.png"
      sp.paintToPng(720.0,3.08,pngDir)
    if twocol:
      sp.setFontSizeForPrint(8.0,469.0)
      pngDir = title+"twocol.png"
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
      #pv.setLineColor(Color.BLUE)
      pv.setMarkStyle(PointsView.Mark.FILLED_CIRCLE)
      pv.setMarkSize(3.4)
      hsmax = max(hsmax,abs(max(hs[ih])),abs(min(hs[ih])))
    elif hs[ih]:
      pv = sp.addPoints(st,hs[ih])
      pv.setLineStyle(PointsView.Line.SOLID)
      pv.setLineColor(Color.RED)
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
      pngDir = title+"onecol.png"
      sp.paintToPng(720.0,3.08,pngDir)
    if twocol:
      sp.setFontSizeForPrint(8.0,469.0)
      pngDir = title+"twocol.png"
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
  maxf = max(f)
  absminf = abs(minf) 
  absmaxf = abs(maxf) 
  if absminf<absmaxf:
    return mul(1.0/maxf,f)
  else:
    return mul(1.0/minf,f)


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
  for i in range(1,nu2):
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

def plotImage(st,sx,f,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.setSize(350,450)
  pv = sp.addPixels(st,sx,f)
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

def getWarpingPartsSpectrums(r,f,g,u):
  nu = len(u)
  nuu = int(r*(nu-1)+1)
  ng = len(g)
  nug = int(r*(ng-1)+1)
  nf = len(f)
  dtgu = 1.0/r
  w = 0.20/r
  warp = Warper()

  uu = warp.upSampleLinear(nu,1.0,0.0,u,nuu,dtgu,0.0)
  ug = warp.upSample(ng,1.0,0.0,g,nug,dtgu,0.0)
  sug  = warp.applyS(dtgu,uu,ug)
  lsug = warp.applyL(r,uu,sug)
  dlsug = warp.subSample(r,nu,lsug)
  stf = Sampling(nf,1.0,0.0)
  stg = Sampling(ng,1.0,0.0)
  stug = Sampling(nug,dtgu,0.0)
  stsug = Sampling(nuu,dtgu,0.0)
  maxf = max(f)
  maxfd2 = maxf/2.0
  amax = [maxf,maxf,maxf,maxf,maxf]
  tmark = [maxfd2,maxfd2,maxfd2,maxfd2,maxfd2]
  print len(f)
  print nf
  plotSequences(Sampling(len(f),1.0,0.0),[f,dlsug],amax=amax,tmark=tmark,\
  labels=["f","dlsug"],title="fdlsug")
  plotSequences(stg,[g],amax=amax,tmark=tmark,\
  labels=["g"],title="g")
  amax = [maxf]
  tmark = [maxfd2]
  plotSequences(stug,[ug],amax=amax,tmark=tmark,\
  labels=["ug"],title="ug")
  plotAmplitudeSpectrumT(stf, dlsug, 0, nf, "amp DLSUg", amax=None)
  plotAmplitudeSpectrumT(stsug, lsug, 0, nuu, "amp LSUg", amax=None)
  plotAmplitudeSpectrumT(stsug, sug, 0, nuu, "amp SUg", amax=None)
  plotAmplitudeSpectrumT(stug, ug, 0, nug, "amp Ug", amax=None)
  plotAmplitudeSpectrumT(stg, g, 0, ng, "amp g", amax=None)


def testWarpingandBandPassFilter():
  nt,ni = 480,1#number of time samples in p and q; number of random impulses in p and q.
  na,ka = 81,-40 #sampling for inverse wavelet A #Note ka<=kc 
  nc,kc = 181,-90 # sampling for wavelet H 
  v = 0.0#The amount of shift between p and q.
  r0,r1 = 2.0,1.0#constant u'(t) $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$2.0
  freqc,decayc,= 0.08,0.05
  mpc = False#is wavelet in f mininmum phase?
  finitec = False#does the wavelet in f have a finite length?
  freqd,decayd,= 0.08,0.05
  mpd = False#is wavelet in f mininmum phase?
  finited = False#does the wavelet in f have a finite length?
  nrmsf = 0.00
  nrmsg = nrmsf
  randomi = False 
  moreps = False 
  p,q,f,g,noisef,noiseg,u,tmin,tmax = createSyntheticLn1D(freqc,decayc,mpc,finitec,\
  freqd,decayd,mpd,finited,r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps)
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


  

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())


