
from imports import *
from edu.mines.jtk.dsp.Conv import *
from java.util import Random

def createSyntheticLn1D(freqc,decayc,mpc,freqd,decayd,mpd,
  r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps):
  u,p,q,itmin,itmax = lnupq(r0,r1,v,nt,ni,randomi,moreps)
  f = addWavelet(freqc,decayc,p,mpc)
  g = addWavelet(freqd,decayd,q,mpd)
  f,noisef = addNoise(nrmsf,42,f)
  g,noiseg = addNoise(nrmsg,43,g)
  return p,q,f,g,noisef,noiseg,u,itmin,itmax

#Sets 2 impulses in p to always be at the same time.
def createSyntheticLn1DSimple(freqc,decayc,mpc,freqd,decayd,mpd,
  r0,r1,v,nrmsf,nrmsg,nt,randomi,moreps):
  u,p,q,itmin,itmax = lnupqsimple(r0,r1,v,nt,randomi,moreps)
  f = addWavelet(freqc,decayc,p,mpc)
  g = addWavelet(freqd,decayd,q,mpd)
  f,noisef = addNoise(nrmsf,42,f)
  g,noiseg = addNoise(nrmsg,43,g)
  return p,q,f,g,noisef,noiseg,u,itmin,itmax

def createSyntheticLn1DFilteredNoise(knoiseupper,width,freqc,decayc,mpc,freqd,decayd,mpd,
  r0,r1,v,nrmsf,nrmsg,nt,ni,randomi,moreps):
  u,p,q,itmin,itmax = lnupq(r0,r1,v,nt,ni,randomi,moreps)
  f = addWavelet(freqc,decayc,p,mpc)
  g = addWavelet(freqd,decayd,q,mpd)
  f,noisef = addNoiseFiltered(knoiseupper,width,nrmsf,42,f)
  g,noiseg = addNoiseFiltered(knoiseupper,width,nrmsg,43,g)
  return p,q,f,g,noisef,noiseg,u,itmin,itmax


def createSyntheticLn2D(freqc,decayc,mpc,freqd,decayd,mpd,
  r0,r1,v,nrmsf,nrmsg,nx,nt,ni,randomi,moreps):
  u1,p1,q1,itmin,itmax = lnupq(r0,r1,v,nt,ni,randomi,moreps)
  u = replicateTrace(nx,u1)
  p = replicateTrace(nx,p1)
  q = replicateTrace(nx,q1)
  f = addWavelet2D(freqc,decayc,p,mpc)
  g = addWavelet2D(freqd,decayd,q,mpd)
  f,noisef = addNoise2D(nrmsf,42,f)
  g,noiseg = addNoise2D(nrmsg,43,g)
  #dump(f)
  print "itmin = "+str(itmin)
  print "itmax = "+str(itmax)
  return p,q,f,g,noisef,noiseg,u,itmin,itmax

def createSyntheticLn2DSimple(freqc,decayc,mpc,freqd,decayd,mpd,
  r0,r1,v,nrmsf,nrmsg,nx,nt,randomi,moreps):
  u1,p1,q1,itmin,itmax = lnupqsimple(r0,r1,v,nt,randomi,moreps)
  u = replicateTrace(nx,u1)
  p = replicateTrace(nx,p1)
  q = replicateTrace(nx,q1)
  f = addWavelet2D(freqc,decayc,p,mpc)
  g = addWavelet2D(freqd,decayd,q,mpd)
  f,noisef = addNoise2D(nrmsf,42,f)
  g,noiseg = addNoise2D(nrmsg,43,g)
  #dump(f)
  print "itmin = "+str(itmin)
  print "itmax = "+str(itmax)
  return p,q,f,g,noisef,noiseg,u,itmin,itmax

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
    dq = (umax-v-50)/float(ni+1)
    ts = rampfloat(dq+v+31,dq,ni)  
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
      if (nexi>0):
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
          #si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = 0#lock
    itmax = int((ts[ni-1]-v)/r0+0.5)
    return u,p,q,itmin,itmax
  else: 
    print "r0>r1" 
    a = (umax-v)/log(r0/r1)
    b = r0*log(r0/r1)/(umax-v)
    for n in range(nt):
      u[n] = a*log(1.0+b*n)+v
    print "a = "+str(a)
    print "b = "+str(b)
    dq = (umax-v-50)/float(ni+1)
    print "dq = "+str(dq)
    ts = rampfloat(dq+v+31,dq,ni)  
    #ts = rampfloat(dq+v+100,dq,ni)  
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
        si.accumulate(tp,rj,nt,1.0,0.0,p)
        if tq<umax:
          si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = 0#lock
    itmax = int(((exp((ts[ni-1]-v)/a)-1.0)/b)+0.5)#lock
    print "tsnim1 = "+str(ts[ni-1])
    return u,p,q,itmin,itmax

#Only meant to be used for 2 impulses
#Will ensure that the impulses in p will be in the same location for all 
#r0 and r1s.
def lnupqsimple(r0,r1,v,nt,randomi,moreps):
  ni = 2
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
    dq = (umax-v-50)/float(ni+1)
    ts = rampfloat(dq+v+31,dq,ni)  
    ran = Random(55)
    #1st impulse
    rj = 1.0
    tp = 75
    tq = tp*r0+v
    print "tq = "+str(tq)
    print "tp = "+str(tp)
    #2nd impulse
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)
    #2nd impulse
    rj = -1.0
    tp = 150
    tq = tp*r0+v
    print "tq = "+str(tq)
    print "tp = "+str(tp)
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)
    itmin = int((ts[0]-v)/r0+0.5)
    itmax = int((ts[ni-1]-v)/r0+0.5)
    return u,p,q,itmin,itmax
  else: 
    print "r0>r1" 
    a = (umax-v)/log(r0/r1)
    b = r0*log(r0/r1)/(umax-v)
    for n in range(nt):
      u[n] = a*log(1.0+b*n)+v
    print "a = "+str(a)
    print "b = "+str(b)
    dq = (umax-v-50)/float(ni+1)
    print "dq = "+str(dq)
    ts = rampfloat(dq+v+31,dq,ni)  
    #ts = rampfloat(dq+v+100,dq,ni)  
    ran = Random(55)
    #1st impulse
    rj = 1.0
    tp = 75
    tq = a*log(1.0+b*tp)+v
    print "tq = "+str(tq)
    print "tp = "+str(tp)
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)
    #2nd impulse
    rj = -1.0
    tp = 150
    tq = a*log(1.0+b*tp)+v
    print "tq = "+str(tq)
    print "tp = "+str(tp)
    si.accumulate(tp,rj,nt,1.0,0.0,p)
    si.accumulate(tq,rj,nt,1.0,0.0,q)

    itmin = 0#lock
    itmax = int(((exp((ts[ni-1]-v)/a)-1.0)/b)+0.5)#lock
    print "tsnim1 = "+str(ts[ni-1])
    return u,p,q,itmin,itmax


def addWavelet(fpeak,decay,p,mp):
  print "############### Add Wavelet ###############"
  w = 2.0*PI*fpeak
  if not mp:
    decay *= 2.0
    w -= 2.0*PI*0.04
  r = exp(-decay)
  print "w = "+str(w)
  print "r = "+str(r)
  a1,a2 = -2.0*r*cos(w),r*r
  print "a1 =",[1,a1,a2]
  #poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w),Cdouble.polar(r,w1),Cdouble.polar(r,-w1)]
  poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
  #poles = [Cdouble.polar(0.8,0)]
  zeros = []
  gain = 1.0
  x = copy(p)
  t = copy(p)
  rcf = RecursiveCascadeFilter(poles,zeros,gain)
  #rcf.applyReverse(p,t)
  rcf.applyForward(p,t)
  SimplePlot.asPoints(t)
  if not mp:
    w = 2.0*PI*(fpeak+0.04)
    print "w = "+str(w)
    print "r = "+str(r)
    a1,a2 = -2.0*r*cos(w),r*r
    print "a2 =",[1,a1,a2]
    poles = [Cdouble.polar(r,w),Cdouble.polar(r,-w)]
    #poles = [Cdouble.polar(0.5,0)]
    zeros = []
    gain = 1.0
    rcf = RecursiveCascadeFilter(poles,zeros,gain)
    #rcf.applyForward(t,x)
    rcf.applyReverse(t,x)
    SimplePlot.asPoints(x)
  else:
    copy(t,x)
  #conv(2,0,[1.0,-0.95],len(x),0,copy(x),len(x),0,x) # attenuate DC
  print "############################################"
  return x

def addWavelet2D(fpeak,decay,p,mp):
  nx = len(p)
  nt = len(p[0])
  f = zerofloat(nt,nx)
  for ix in range(0,nx):
    f[ix] = addWavelet(fpeak,decay,p[ix],mp)
  return f

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
  n = len(f[0])
  r = Random(seed)
  g = mul(2.0,sub(randfloat(r,n),0.5))
  g2 = zerofloat(n,nx)
  rgf = RecursiveGaussianFilter(1.0)
  rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  #print "nrms = "+str(nrms)
  #print "rmsnoise/rmssignal = "+str(rms(g)/rms(f))
  for i in range(0,nx):
    f[i] = add(f[i],g)
    g2[i] = g
  return f,g2

#Returns the signal added with noise and the noise that 
#was added.
def addNoiseFiltered(knoiseupper,width,nrms, seed, f):
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
  gfiltered = zerofloat(n)
  bpf = BandPassFilter(0.0,knoiseupper,width,0.01)
  bpf.apply(g,gfiltered)
  return add(f,gfiltered),gfiltered

#Copies the same trace to create a 2D image
def replicateTrace(nx,f1):
  nt = len(f1)
  f = zerofloat(nt,nx)
  for i in range(0,nx):
    f[i] = f1
  return f

def createSyntheticCos(freq,decay,mp,
  r0,r1,v,noise,nrmsf,nrmsg,nt,ni):
  u,p,q,itmin,itmax = cosupq(r0,r1,v,nt,ni)
  f = addWavelet(freq,decay,p,mp)
  g = addWavelet(freq,decay,q,mp)
  if noise:
    f = addNoise(nrmsf,42,f)
    g = addNoise(nrmsg,43,g)
  return p,q,f,g,u,itmin,itmax

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



