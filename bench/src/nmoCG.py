#############################################################################
# Demo nmo stretch

from imports import *
from linalgebra import ToeplitzRecursion 
from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.interp.CubicInterpolator import *
from edu.mines.jtk.util.ArrayMath import *
from warp import WaveletNmo,WarpedWavelet,ShapingFilter
from data import SUDataGrabber
from dspCG import Correlator

#############################################################################

#pngDir = "./png"
pngDir = None

def main(args):
  #Define gather file locations and vnmo and tnmo values
  ################################################################
  cmp300 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=300.strip"
  cmp301 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=301.strip"
  cmp1298 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=1298.strip"
  cmp1299 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=1299.strip"
  cmp1300 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=1300.strip"
  cmp1301 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=1301.strip"
  cmp1302 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=1302.strip"
  cmp1700 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=1700.strip"
  cmp1701 = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/m1cdp=1701.strip"
  tnmo300 = [0.0570429,0.456343,1.35952,2.09157,2.7856,3.73631,4.90569]
  vnmo300 = [1.51219,1.51597,1.848,2.1465,2.36416,2.58181,2.67509]
  tnmo1300 =     [0.0380286,0.475358,0.950715,1.25,1.75882,2.33876,3.54617,4.35428,5.68528]
  vnmo1300 =     [1.50597  ,1.52462 ,1.79203 ,1.87,2.0    ,2.28953,2.45122,2.65644,2.69997]
  tnmo1700 = [0.0475358,0.456343,0.988744,1.71129,2.54792,3.30849,4.36378,5.50464]
  vnmo1700 = [1.50597,1.52462,1.77959,2.09675,2.36416,2.48853,2.56938,2.61291]
#ts=0.010,0.646,0.846,1.150,1.511,1.825,2.490,2.985
#  vs=1.510,1.575,1.687,1.817,1.938,1.980,2.446,2.735
  ################################################################
  na,ka = 11, 0 # sampling for inverse wavelet a
  nh,kh = 201,-50 # sampling for wavelet h
  goPredictiveDeconvCompareSynthetic()
  #goPredictiveDeconvCompare(cmp1300,tnmo1300,vnmo1300,na,ka,nh,kh)
  #goDrHaleCorrTest()
  #test2Corr(cmp1300)
  #test1Corr(5,5)
  #goEstimateWaveletFromVikingGrabenW(cmp1300,tnmo1300,vnmo1300,na,ka,nh,kh)
  #goEstimateWaveletFromVikingGrabenNoW(cmp1300,tnmo1300,vnmo1300,na,ka,nh,kh)
  #goPlotWaveletFor5AdjacentCMPs(cmp1298,cmp1299,cmp1300,cmp1301,cmp1302,tnmo1300,vnmo1300,na,ka,nh,kh)
  #goPlotWaveletFor2DistantCMPs(cmp300,cmp1700,tnmo300,tnmo1700,vnmo300,vnmo1700,na,ka,nh,kh)
  #goPlotEnhancedNMOFor5AdjacentCMPs(cmp1298,cmp1299,cmp1300,cmp1301,cmp1302,tnmo1300,vnmo1300,na,ka,nh,kh)
  #compareStacksDiffNMOCor()
  #goCmpGatherWithKnownWavelet()
  #goEstimateWaveletForOneOffsetVikingGraben()
  #goEstimateWaveletForOneOffset()

def goPredictiveDeconvCompare(cmp1,tnmo1,vnmoRaw1,na,ka,nh,kh):
  #Input data
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first
  hp = SUDataGrabber.grab2DFile(cmp1, nx, nt)
  #apply gain due to compensate for geometric spreading
  hp = tpow(3.0,st,hp)

  #Vnmo values for viking-graben dataset
  vc = 1.5
  #vnmoRaw = [vc,vc,vc,vc,vc,vc,vc,vc,]
  vnmoInterp = interpolateVel(tnmo1,st,vnmoRaw1)

  #################### Enhanced NMO Estimation #################
  wn = WaveletNmo(st,sx,vnmoInterp)
  aenha = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  henha = wn.getWaveletH(na,ka,aenha,nh,kh); # estimate wavelet
  g = wn.applyHNmoA(na,ka,aenha,nh,kh,henha,hp) #convolve inverse, apply NMO, convolve wavelet
  #apply NMO to original data
  #e = wn.applyNmo(hp)
  
  #tmin,tmax,perc = 0.0,2.5,98.0
  #plotSequence(Sampling(na,st.delta,kh*st.delta),aenha,title="inverse wavelet")
  #plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(henha),title="estimated wavelet")
  #plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  #plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  #plotGather(st,sx,hp,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  #plotSequence(Sampling(na,st.delta,ka*st.delta),normalize(a),title="inverse")
 # plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(h),
               #title="estimated wavelet")

 
  #######################################################################
  #################### Predictive Deonvolution Estimation ##########
  napred = na+1 #preddecon reduces filter length from na to na-1, +1 compensates
  #hp = tpow(3.0,st,hp)
  z = zerofloat(napred)
  zsum = zerofloat(napred)
  khp = 0 
  kz = 0 
  for i in range(0,nx): #sum up autocorrelations of gather
    xcor(nt,khp,hp[i],nt,khp,hp[i],napred,kz,z)
    add(zsum,z,zsum)
  acorfunc = zerofloat(na)
  g = zerofloat(napred)
  for i in range(0,na):
    g[i] = zsum[i+1]
    acorfunc[i] = zsum[i]
  #g[0] = 1

  tr = ToeplitzRecursion(acorfunc,g)
  apred = tr.solve()
 
  one = [1.0]
  hpred = ShapingFilter.design(nh,kh,na,ka,normalize(apred),1,0,one); # estimate wavelet
  
  #######################################################################
  #plotSequence(Sampling(na,st.delta,kh*st.delta),normalize(ak),xmax=1.2,title="known inverse")
  #plotSequence(Sampling(na,st.delta,kh*st.delta),normalize(apred),xmax=1.2,title="Pred Decon Inverse")
  #plotSequence(Sampling(na,st.delta,kh*st.delta),normalize(aenha),xmax=1.2,title="enhanced inverse wavelet")
  #plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(hknwn),xmax=1.2,title="known wavelet")
  #plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(hpred),xmax=1.2,title="Pred Decon Wavelet")
  #plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(henha),xmax=1.2,title="enhanced estimated wavelet")

  samph = Sampling(nh,st.delta,kh*st.delta)
  sampa = Sampling(na,st.delta,ka*st.delta)

  plot2Sequences(samph,samph,normalize(hpred),normalize(henha),xmax=1.2,title="Estimated Wavelet, r pred, g enchanced")
  plot2Sequences(sampa,sampa,normalize(apred),normalize(aenha),xmax=1.5,title="Estimated Inverse, r pred, g enchanced")

def goPredictiveDeconvCompareSynthetic():
  #################### Build wavelet and CMP gather with const v##########
  st = Sampling(501,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(201,0.010,0.0); nx,dx,fx = sx.count,sx.delta,sx.first
  nref,vnmo = 1,2.0 # number of reflectors and NMO velocity
  freq,decay = 30.0,0.1 # peak frequency and decay for wavelet
  p = makeCmpReflections(vnmo,nref,st,sx) # cmp gather without wavelet
  hp = addArWavelet(freq,decay,st,sx,p) # cmp gather with wavelet
  na,ka = 3,0 # sampling for inverse wavelet a
  ak = zerofloat(na) # array for the known inverse wavelet a
  r,w = exp(-decay),2.0*PI*freq*st.delta # radius and frequency of poles
  a1,a2 = -2.0*r*cos(w),r*r # coefficients for inverse wavelet
  ak[0-ka] = 1.0
  ak[1-ka] = a1
  ak[2-ka] = a2
  nh,kh = 100,-20 # sampling for wavelet h
  one = [1.0]
  hknwn = ShapingFilter.design(nh,kh,na,ka,ak,1,0,one); # estimate wavelet
  #######################################################################
  #################### Predictive Deonvolution Estimation ##########
  napred = na+1 #preddecon reduces filter length from na to na-1, +1 compensates
  #hp = tpow(3.0,st,hp)
  z = zerofloat(napred)
  zsum = zerofloat(napred)
  khp = 0 
  kz = 0 
  for i in range(0,nx): #sum up autocorrelations of gather
    xcor(nt,khp,hp[i],nt,khp,hp[i],napred,kz,z)
    add(zsum,z,zsum)
  acorfunc = zerofloat(na)
  g = zerofloat(napred)
  for i in range(0,na):
    g[i] = zsum[i+1]
    acorfunc[i] = zsum[i]
  print acorfunc
  #g[0] = 1

  tr = ToeplitzRecursion(acorfunc,g)
  apred = tr.solve()
 
  one = [1.0]
  hpred = ShapingFilter.design(nh,kh,na,ka,normalize(apred),1,0,one); # estimate wavelet
  #######################################################################
  #################### Enhanced NMO Estimation #########################
  wn = WaveletNmo(st,sx,vnmo)
  aenha = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  henha = wn.getWaveletH(na,ka,aenha,nh,kh); # estimate wavelet
  #######################################################################
  #plotSequence(Sampling(na,st.delta,kh*st.delta),normalize(ak),xmax=1.2,title="known inverse")
  #plotSequence(Sampling(na,st.delta,kh*st.delta),normalize(apred),xmax=1.2,title="Pred Decon Inverse")
  #plotSequence(Sampling(na,st.delta,kh*st.delta),normalize(aenha),xmax=1.2,title="enhanced inverse wavelet")
  #plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(hknwn),xmax=1.2,title="known wavelet")
  #plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(hpred),xmax=1.2,title="Pred Decon Wavelet")
  #plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(henha),xmax=1.2,title="enhanced estimated wavelet")

  samph = Sampling(nh,st.delta,kh*st.delta)
  sampa = Sampling(na,st.delta,ka*st.delta)

  plot3Sequences(samph,samph,samph,normalize(hknwn),normalize(hpred),normalize(henha),xmax=1.2,title="Estimated Wavelet, r known, g pred, b enchanced")
  plot3Sequences(sampa,sampa,sampa,normalize(ak),normalize(apred),normalize(aenha),xmax=1.5,title="Estimated Inverse, r known, g pred, b enchanced")

def goEstimateWaveletFromVikingGrabenW(cmp1,tnmo1,vnmoRaw1,na,ka,nh,kh):
  #Input data
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first
  hp = SUDataGrabber.grab2DFile(cmp1, nx, nt)
  #apply gain due to compensate for geometric spreading
  hp = tpow(3.0,st,hp)

  #Vnmo values for viking-graben dataset
  vc = 1.5
  #vnmoRaw = [vc,vc,vc,vc,vc,vc,vc,vc,]
  vnmoInterp = interpolateVel(tnmo1,st,vnmoRaw1)

  wn = WaveletNmo(st,sx,vnmoInterp)
  a = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
  g = wn.applyHNmoA(na,ka,a,nh,kh,h,hp) #convolve inverse, apply NMO, convolve wavelet
  #apply NMO to original data
  e = wn.applyNmo(hp)
  
  d = sub(g,e)

  tmin,tmax,perc = 0.0,2.5,98.0
  plotSequence(Sampling(na,st.delta,kh*st.delta),a,title="inverse wavelet")
  plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(h),title="estimated wavelet")
  #plotGather(st,sx,d,tmin=tmin,tmax=tmax,perc=perc,title="Difference")
  #plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  #plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  #plotGather(st,sx,hp,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  #plotSequence(Sampling(na,st.delta,ka*st.delta),normalize(a),title="inverse")
 # plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(h),
               #title="estimated wavelet")

def goEstimateWaveletFromVikingGrabenNoW(cmp1,tnmo1,vnmoRaw1,na,ka,nh,kh):
  #Input data
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first
  hp = SUDataGrabber.grab2DFile(cmp1, nx, nt)
  #apply gain due to compensate for geometric spreading
  hp = tpow(3.0,st,hp)

  #Vnmo values for viking-graben dataset
  vc = 1.5
  #vnmoRaw = [vc,vc,vc,vc,vc,vc,vc,vc,]
  vnmoInterp = interpolateVel(tnmo1,st,vnmoRaw1)

  wn = WaveletNmo(st,sx,vnmoInterp)
  a = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
  g = wn.applyHNmoA(na,ka,a,nh,kh,h,hp) #convolve inverse, apply NMO, convolve wavelet
  #apply NMO to original data
  e = wn.applyNmo(hp)
  
  d = sub(g,e)

  tmin,tmax,perc = 0.0,2.5,98.0
  #plotGather(st,sx,d,tmin=tmin,tmax=tmax,perc=perc,title="Difference")
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  #plotGather(st,sx,hp,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  #plotSequence(Sampling(na,st.delta,ka*st.delta),normalize(a),title="inverse")
 # plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(h),
               #title="estimated wavelet")

def goPlotEnhancedNMOFor5AdjacentCMPs(cmp1,cmp2,cmp3,cmp4,cmp5,tnmo,vnmoRaw,na,ka,nh,kh):
  #Input data
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first

  vnmoInterp = interpolateVel(tnmo,st,vnmoRaw)
  wn = WaveletNmo(st,sx,vnmoInterp)

  hp1 = SUDataGrabber.grab2DFile(cmp1, nx, nt)
  g1 = goEnhancedNMOFromVikingGraben(hp1,st,sx,wn,na,ka,nh,kh) 
  hp2 = SUDataGrabber.grab2DFile(cmp2, nx, nt)
  g2 = goEnhancedNMOFromVikingGraben(hp2,st,sx,wn,na,ka,nh,kh) 
  hp3 = SUDataGrabber.grab2DFile(cmp3, nx, nt)
  g3 = goEnhancedNMOFromVikingGraben(hp3,st,sx,wn,na,ka,nh,kh) 
  hp4 = SUDataGrabber.grab2DFile(cmp4, nx, nt)
  g4 = goEnhancedNMOFromVikingGraben(hp4,st,sx,wn,na,ka,nh,kh) 
  hp5 = SUDataGrabber.grab2DFile(cmp5, nx, nt)
  g5 = goEnhancedNMOFromVikingGraben(hp5,st,sx,wn,na,ka,nh,kh) 
  tmin,tmax,perc = 0.0,2.5,98.0
  plotGather(st,sx,g1,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1298 improved NMO")
  plotGather(st,sx,cmp1,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1298 input gather")
  plotGather(st,sx,g2,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1299 improved NMO")
  plotGather(st,sx,cmp2,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1299 input gather")
  plotGather(st,sx,g3,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1300 improved NMO")
  plotGather(st,sx,cmp3,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1300 input gather")
  plotGather(st,sx,g4,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1301 improved NMO")
  plotGather(st,sx,cmp4,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1301 input gather")
  plotGather(st,sx,g5,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1302 improved NMO")
  plotGather(st,sx,cmp5,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1302 input gather")

def goPlotWaveletFor5AdjacentCMPs(cmp1,cmp2,cmp3,cmp4,cmp5,tnmo,vnmoRaw,na,ka,nh,kh):
  #Input data
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first

  vnmoInterp = interpolateVel(tnmo,st,vnmoRaw)
  wn = WaveletNmo(st,sx,vnmoInterp)

  hp1 = SUDataGrabber.grab2DFile(cmp1, nx, nt)
  w1 = goEstimateWaveletFromVikingGraben(hp1,st,sx,wn,na,ka,nh,kh) 
  hp2 = SUDataGrabber.grab2DFile(cmp2, nx, nt)
  w2 = goEstimateWaveletFromVikingGraben(hp2,st,sx,wn,na,ka,nh,kh) 
  hp3 = SUDataGrabber.grab2DFile(cmp3, nx, nt)
  w3 = goEstimateWaveletFromVikingGraben(hp3,st,sx,wn,na,ka,nh,kh) 
  hp4 = SUDataGrabber.grab2DFile(cmp4, nx, nt)
  w4 = goEstimateWaveletFromVikingGraben(hp4,st,sx,wn,na,ka,nh,kh) 
  hp5 = SUDataGrabber.grab2DFile(cmp5, nx, nt)
  w5 = goEstimateWaveletFromVikingGraben(hp5,st,sx,wn,na,ka,nh,kh) 
  sp = SimplePlot()
  sp.addPoints(Sampling(nh,st.delta,kh*st.delta),w1)
  sp.addPoints(Sampling(nh,st.delta,kh*st.delta),w2)
  sp.addPoints(Sampling(nh,st.delta,kh*st.delta),w3)
  sp.addPoints(Sampling(nh,st.delta,kh*st.delta),w4)
  sp.addPoints(Sampling(nh,st.delta,kh*st.delta),w5)
  sp.addTitle("5 Adjacent CMP wavelets")

def goPlotWaveletFor2DistantCMPs(cmp1,cmp2,tnmo1,tnmo2,vnmoRaw1,vnmoRaw2,na,ka,nh,kh):
  #Input data
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first
  vnmoInterp1 = interpolateVel(tnmo1,st,vnmoRaw1)
  wn1 = WaveletNmo(st,sx,vnmoInterp1)

  vnmoInterp2 = interpolateVel(tnmo2,st,vnmoRaw2)
  wn2 = WaveletNmo(st,sx,vnmoInterp2)

  hp1 = SUDataGrabber.grab2DFile(cmp1, nx, nt)
  g1 = goEnhancedNMOFromVikingGraben(hp1,st,sx,wn1,na,ka,nh,kh)
  w1 = goEstimateWaveletFromVikingGraben(hp1,st,sx,wn1,na,ka,nh,kh) 
  hp2 = SUDataGrabber.grab2DFile(cmp2, nx, nt)
  g2 = goEnhancedNMOFromVikingGraben(hp2,st,sx,wn2,na,ka,nh,kh)
  w2 = goEstimateWaveletFromVikingGraben(hp2,st,sx,wn2,na,ka,nh,kh) 
  tmin,tmax,perc = 0.0,6.0,98.0
  sp = SimplePlot()
  sp.addPoints(Sampling(nh,st.delta,kh*st.delta),w1)
  sp.addPoints(Sampling(nh,st.delta,kh*st.delta),w2)
  sp.addTitle("Two Wavelets from Distant CMP Gathers")
  plotGather(st,sx,g1,tmin=tmin,tmax=tmax,perc=perc,title="CMP 300 improved NMO")
  plotGather(st,sx,g2,tmin=tmin,tmax=tmax,perc=perc,title="CMP 1700 improved NMO")

def goEnhancedNMOFromVikingGraben(hp,st,sx,wn,na,ka,nh,kh):
  #Input data
  nt,dt,ft = st.count,st.delta,st.first
  nx,dx,fx = sx.count,sx.delta,sx.first
  #apply gain due to compensate for geometric spreading
  hp = tpow(3.0,st,hp)
  
  a = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
  g = wn.applyHNmoA(na,ka,a,nh,kh,h,hp) #convolve inverse, apply NMO, convolve wavelet
  return g;

def goEstimateWaveletFromVikingGraben(hp,st,sx,wn,na,ka,nh,kh):
  #Input data
  nt,dt,ft = st.count,st.delta,st.first
  nx,dx,fx = sx.count,sx.delta,sx.first
  #hp = SUDataGrabber.grab2DFile(fileLoc, nx, nt)
  #apply gain due to compensate for geometric spreading
  hp = tpow(3.0,st,hp)
  
  a = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
  return h;

def compareStacksDiffNMOCor():
  fileLoc = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/mcdp=1300.strip"
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first
  hp = SUDataGrabber.grab2DFile(fileLoc, nx, nt)
  #apply gain due to compensate for geometric spreading
  hp = tpow(3.0,st,hp)

  #Vnmo values for viking-graben dataset
  tnmo = [0.0380286,0.475358,0.950715,1.75882,2.33876,3.54617,4.35428,5.68528]
  t = []
  for i in range(0,nt):
    t.append(i*dt)
  vnmoRaw = [1.50597,1.52462,1.79203,2.04078,2.28953,2.45122,2.65644,2.69997]
  #Build squared nmo values
  nvnmo = len(vnmoRaw)
  vnmoRawSq = []
  for i in range(0,nvnmo):
    vnmoRawSq.append(vnmoRaw[i]*vnmoRaw[i])

  ci = CubicInterpolator(Method.MONOTONIC, tnmo, vnmoRawSq)
  vnmoInterpSq = ci.interpolate(t)

  #Square Root vnmo squared values
  vnmoInterp = []
  for i in range(0,nt):
    vnmoInterp.append(sqrt(vnmoInterpSq[i]))

  na,ka = 11, 0 # sampling for inverse wavelet a
  nh,kh = 201,-50 # sampling for wavelet h
  wn = WaveletNmo(st,sx,vnmoInterpSq)
  a = wn.getInverseA(na,ka,hp) # estimate inverse wavelet
  print "a ="; dump(a);
  h = wn.getWaveletH(na,ka,a,nh,kh); # estimate wavelet
  g = wn.applyHNmoA(na,ka,a,nh,kh,h,hp) #convolve inverse, apply NMO, convolve wavelet
  #apply NMO to original data
  e = wn.applyNmo(hp)
  sg = wn.stack(g)
  se = wn.stack(e)

  d = sub(sg,se)

  tmin,tmax,perc = 0.0,6.0,98.0
  plotSequence(st,se,None,title="conventional NMO corrected stack")
  plotSequence(st,sg,None,title="improved NMO corrected stack")

def goEstimateWaveletForOneOffset():
  """ Estimates wavelet from a non-zero-offset and zero-offset trace """
  fileLoc = "C:/Users/Chris/Documents/CWP/Research/research/vikinggrabenCMP/mcdp=1300.strip"
  st = Sampling(1500,0.004,0.0); nt,dt,ft = st.count,st.delta,st.first
  sx = Sampling(60,0.0496,.262); nx,dx,fx = sx.count,sx.delta,sx.first
  hp = SUDataGrabber.grab2DFile(fileLoc, nx, nt)
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  hp = tpow(3.0,st,hp)

 #Vnmo values for viking-graben dataset
  tnmo = [0.0380286,0.475358,0.950715,1.75882,2.33876,3.54617,4.35428,5.68528]
  t = []
  for i in range(0,nt):
    t.append(i*dt)
  vnmoI = [1.50597,1.52462,1.79203,2.04078,2.28953,2.45122,2.65644,2.69997]
  ci = CubicInterpolator(Method.MONOTONIC, tnmo, vnmoI)
  vnmoF = ci.interpolate(t)

  ixx,ixy = 59,0 # indices for non-zero-offset and zero-offset traces
  offset = sx.getValue(ixx) # the non-zero offset
  x,y = hp[ixx],hp[ixy] # x is non-zero-offset trace; y is zero-offset trace
  na,ka = 11,0 # number of samples and index of 1st sample for inverse
  nh,kh = 201,-50 # sampling for wavelet h
  ww = WarpedWavelet(WarpedWavelet.Nmo(st,offset,vnmoF))
  ae = ww.estimateInverse(na,ka,x,y) # the estimated wavelet inverse
  print "a = ",; dump(ae)
  h = ww.estimateWavelet(na,ka,ae,nh,kh); # estimate wavelet
  wn = WaveletNmo(st,sx,vnmoF)
  g = wn.applyHNmoA(na,ka,ae,nh,kh,h,hp)

  #apply NMO to original data
  e = wn.applyNmo(hp)
  tmin,tmax,perc = 0.0,6.0,98.0
  plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(h),title="estimated wavelet")
  plotGather(st,sx,hp,tmin=tmin,tmax=tmax,perc=perc,title="input gather")
  plotGather(st,sx,e,tmin=tmin,tmax=tmax,perc=perc,title="conventional NMO")
  plotGather(st,sx,g,tmin=tmin,tmax=tmax,perc=perc,title="improved NMO")
  #for a in [ak,ae]:
  #  if a is ak:
  #    title = "known wavelet"
  #  else:
  #    title = "estimated wavelet"
  #  print title
  #  ax = zerofloat(nt)
  #  ay = zerofloat(nt)
  #  conv(na,ka,a,nt,0,x,nt,0,ax)
  #  conv(na,ka,a,nt,0,y,nt,0,ay)
  #  sax = applyNmo1(offset,vnmo,st,ax)
  #  rms = sqrt(sum(pow(sub(ay,sax),2))/nt)
  #  print "a = ",; dump(a)
  #  print "rms error =",rms
  #  nh,kh = 100,-20 # number of samples, index of 1st sample for wavelet h
  #  h = ww.estimateWavelet(na,ka,a,nh,kh);
  #  plotSequence(Sampling(nh,st.delta,kh*st.delta),normalize(h),title=title)
  #return
  #hsax = zerofloat(nt)
  #conv(nh,kh,h,nt,0,sax,nt,0,hsax)
  #sx = applyNmo1(offset,vnmo,st,x)
  #plotSequence(st,x,3.0,"x")
  #plotSequence(st,sx,3.0,"sx")
  #plotSequence(st,hsax,3.0,"hsax")
  #plotSequence(st,y,3.0,"y")

def goCmpGatherWithKnownWavelet():
  n = 501
  st = Sampling(501,0.004,0.0)
  sx = Sampling(201,0.010,0.0)
  nref,vnmo = 100,2.0
  freq,decay = 30.0,0.1
  p = makeCmpReflections(vnmo,nref,st,sx)
  hp = addArWavelet(freq,decay,st,sx,p)
  shp = applyNmo(vnmo,st,sx,hp)
  sp = applyNmo(vnmo,st,sx,p)
  hsp = addArWavelet(freq,decay,st,sx,sp)
  #plotGather(st,sx,p)
  plotGather(st,sx,hp,"CMP gather")
  plotGather(st,sx,shp,"NMO with wavelet distortion")
  plotGather(st,sx,hsp,"Reduced wavelet distortion")


#take tnmo and vnmo and interpolates vnmo^2 and turns it back
#to vnmo
def interpolateVel(tnmo, st, vnmoRaw):
  nt = st.getCount();
  dt = st.getDelta();
  t = []
  for i in range(0,nt):
    t.append(i*dt)
  #Build squared nmo values
  nvnmo = len(vnmoRaw)
  vnmoRawSq = []
  for i in range(0,nvnmo):
    vnmoRawSq.append(vnmoRaw[i]*vnmoRaw[i])

  ci = CubicInterpolator(Method.MONOTONIC, tnmo, vnmoRawSq)
  vnmoInterpSq = ci.interpolate(t)

  #Square Root vnmo squared values
  vnmoInterp = []
  for i in range(0,nt):
    vnmoInterp.append(sqrt(vnmoInterpSq[i]))
  return vnmoInterp


def tpow(power,st,f):
  """Applies t^power gain."""
  nt,dt,ft = st.count,st.delta,st.first
  nx = len(f)
  tp = pow(rampfloat(ft,dt,nt),power) # sampled times raised to power
  g = zerofloat(nt,nx)
  for ix in range(nx):
    mul(tp,f[ix],g[ix])
  return g

def normalize(h):
  return div(h,max(h))

def estimateInverseWavelet(na,ka,offset,vnmo,st,x,y):
  nmo = WarpedWavelet.Nmo(st,offset,vnmo)
  ww = WarpedWavelet(nmo)
  a = ww.estimateInverse(na,ka,x,y)
  return a

def applyNmo1(offset,vnmo,st,p):
  q = zerofloat(st.count)
  nmo = WarpedWavelet.Nmo(st,offset,vnmo)
  nmo.apply(p,q)
  return q

def plotGather(st,sx,p,tmin=None,tmax=None,perc=None,title=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if title:
    sp.setTitle(title)
  sp.setHLabel("Offset (km)")
  sp.setVLabel("Time (s)")
  sp.setSize(400,750)
  sp.setVLimits(tmin,tmax)
  pv = sp.addPixels(st,sx,p)
  if perc:
    pv.setPercentiles(100-perc,perc)

def plot2Sequences(sr,sg,r,g,xmax=None,title=None):
  pvr = PointsView(sr,r)
  pvr.setLineColor(Color.RED)
  pvg = PointsView(sg,g)
  pvg.setLineColor(Color.GREEN)

  pp = PlotPanel()
  pp.addTiledView(pvr)
  pp.addTiledView(pvg)
  if xmax:
    pp.setVLimits(-xmax,xmax)
  if title:
    pp.setTitle(title)

  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pf.setVisible(True)

def plot3Sequences(sr,sg,sb,r,g,b,xmax=None,title=None):
  pvr = PointsView(sr,r)
  pvr.setLineColor(Color.RED)
  pvg = PointsView(sg,g)
  pvg.setLineColor(Color.GREEN)
  pvb = PointsView(sb,b)
  pvb.setLineColor(Color.BLUE)

  pp = PlotPanel()
  pp.addTiledView(pvr)
  pp.addTiledView(pvg)
  pp.addTiledView(pvb)
  if xmax:
    pp.setVLimits(-xmax,xmax)
  if title:
    pp.setTitle(title)

  pf = PlotFrame(pp)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pf.setVisible(True)




def plotSequence(s,x,xmax=None,title=None):
  sp = SimplePlot.asPoints(s,x)
  if xmax:
    sp.setVLimits(-xmax,xmax)
  if title:
    sp.setTitle(title)

def applyNmo(vel,st,sx,p):
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  t0 = add(0.01*dt,rampfloat(ft,dt,nt))
  q = copy(p)
  si = SincInterp.fromErrorAndFrequency(0.01,0.4)
  for jx in range(nx):
    xj = sx.getValue(jx)
    cj = (xj*xj)/(vel*vel)
    ti = sqrt(add(mul(t0,t0),cj))
    si.interpolate(nt,dt,ft,p[jx],nt,ti,q[jx])
    q[jx] = mul(q[jx],div(t0,ti))
  return q

def makeCmpReflections(vel,nref,st,sx):
  nt,nx = st.count,sx.count
  dt,dx = st.delta,sx.delta
  ft,fx = st.first,sx.first
  p = zerofloat(nt,nx)
  ts = add(ft,mul((nt-1)*dt,randfloat(nref)))
  rs = sub(mul(2.0,randfloat(nref)),1.0)
  si = SincInterp.fromErrorAndFrequency(0.01,0.45)
  for jx in range(nx):
    xj = sx.getValue(jx)
    cj = (xj*xj)/(vel*vel)
    for jr in range(nref):
      tj = ts[jr]
      rj = rs[jr]
      tj = sqrt(tj*tj+cj)
      si.accumulate(tj,rj,nt,dt,ft,p[jx])
  return p

def addArWavelet(fpeak,decay,st,sx,p):
  r = exp(-decay)
  w = 2.0*PI*fpeak*st.delta
  a1,a2 = -2.0*r*cos(w),r*r
  #print "a1 =",a1," a2 =",a2
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  zeros = []
  gain = 1.0
  x = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  rcf.apply1Forward(p,x)
  return x

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),9.0)



#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
