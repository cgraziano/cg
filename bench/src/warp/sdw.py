from imports import *
from warp import *
from dnp import *
from viewer import *
from wwarp import *
import tensors
import plotUtils

dataDir = "/Users/Chris/data/gbc/dat/"
sigma1,sigma2,sigma3 = 8.0,8.0,8.0
smoothingDegree = 2.0
n1,d1,f1 = 2000,0.002,0.0
n2,d2,f2 =  150,0.033531,0.0
n3,d3,f3 =  145,0.033531,0.0
n1ps = 1501#subset of PS data to use
s1 = Sampling(n1,d1,f1)
s2 = Sampling(n2,d2,f2)
s3 = Sampling(n3,d3,f3)
s1ps = Sampling(n1ps,s1.getDelta(),s1.getFirst())



# Contstraints for time shift slopes, which are physically related to interval
# Vp/Vs ratios. Compute slope parameters from vpvsMin/Max.
vpvsAvgGuess = 1.9
vpvsMin,vpvsMax = 1.5,2.5
r1min = (vpvsMin-1.0)/2
r1max = (vpvsMax-1.0)/2

# Constraints for horizontal directions. A better match between the PS and PP
# images can be achieved by loosening these constraints, but I'm not sure that
# we want to correct for all lateral variations seen in the PS image but not
# in the PP image. These differences may not be related to Vp/Vs ratios.
r2min,r2max,dg2 = -0.15,0.15,15
r3min,r3max,dg3 = -0.15,0.15,15

#Plotting parameters
bwr = ColorMap.HUE
gry = ColorMap.GRAY
jet = ColorMap.JET
iMap = gry

k1,k2,k3 = 345,83,48 # default indices for plotting
slices = [k1,k2,k3]
he0 = 332 # height of time slice panel in plots
iClips = [-1.5,1.5]
uClips = [0.0,0.75]
vClips = [1.5,2.2]


def main(args):
  #structureSmoothing("pp",s1,s2,s3,smoothingDegree,sigma1,sigma2,sigma3,"pp_smooth")
  #structureSmoothing("ps1",s1,s2,s3,smoothingDegree,sigma1,sigma2,sigma3,"ps1_smooth")
  smoothDynamicWarping("pp","ps1",s1,s2,s3)
  #warpPSToPP("u_sag_linear","ps1_smooth")
  sg = getData("sg",s1,s2,s3)
  u = getData("u_sag",s1,s2,s3)
  g = getData("ps1",s1,s2,s3)
  f = getData("pp",s1,s2,s3)
  print "min sg = "+str(min(sg))
  print "max sg = "+str(max(sg))
  plotUtils.plotPP3(sg,title="sg",s1=s1,s2=s2,s3=s3,clips1=[min(sg),max(sg)],
          cmap1=iMap,label1="PP time (s)",cbw=100,slices=slices,
          he0=he0)
  plotUtils.plotPP3(u,title="u",s1=s1,s2=s2,s3=s3,clips1=[min(u),max(u)],
          cmap1=iMap,label1="PP time (s)",cbw=100,slices=slices,
          he0=he0)
  plotUtils.plotPP3(g,title="g",s1=s1,s2=s2,s3=s3,clips1=[min(g),max(g)],
          cmap1=iMap,label1="PS time (s)",cbw=100,slices=slices,
          he0=he0)
  plotUtils.plotPP3(f,title="f",s1=s1,s2=s2,s3=s3,clips1=[min(f),max(f)],
          cmap1=iMap,label1="PP time (s)",cbw=100,slices=slices,
          he0=he0)


#Don't include extensions in file names.
def structureSmoothing(inputFileName,s1,s2,s3,smoothingDegree,sigma1,sigma2,sigma3,outputFileName):
  data = getData(inputFileName,s1,s2,s3)
  dataSmooth = tensors.createTensorsAndSmooth(s1,s2,s3,data,inputFileName,dataDir,smoothingDegree,sigma1,sigma2,sigma3)
  writeImage(outputFileName,dataSmooth)

#Do not include smooth ending in fileNames.
def smoothDynamicWarping(ppFileName,psFileName,s1,s2,s3):
  print "begin smooth dynamic warping"
  # subset for PS data
  pp = getData(ppFileName+"_smooth",s1,s2,s3)
  ps1 = getData(psFileName+"_smooth",s1,s2,s3)
  dw = DynamicWarpingC.fromVpVs(s1,s1ps,vpvsAvgGuess,0.0,s2,s3)
  dw.setStrainLimits(r1min,r1max,r2min,r2max)
  setPlotVars(dw)

  # Warping using a structure aligned grid.
  # dg1=50, ng=13 - choose 13 strongest reflectors, that must be a minimum of
  # 50 samples apart. The first and last time samples are 2 of the 13 layers.
  # Fasle to read shifts from disk, True to find shifts.
  us = goShiftsSAG(s1,pp,s1ps,ps1,dw,50,13,True)



def getData(filename,s1,s2,s3):
  n1 = s1.getCount()
  n2 = s2.getCount()
  n3 = s3.getCount()
  f = readImage(filename,n1,n2,n3)
  zm = ZeroMask(f)
  gain(100,f)
  WarpUtils.normalize(f,-1.5,1.5)
  zm.apply(0.0,f)
  return f

def readImage(fileName,n1,n2,n3):
  print dataDir+fileName+".dat"
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(dataDir+fileName+".dat")
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileName,x):
  print dataDir+fileName+".dat"
  aos = ArrayOutputStream(dataDir+fileName+".dat")
  aos.writeFloats(x)
  aos.close()


#Don't include extensions in file name.
def view(fileName,n1,n2,n3):
  f = readImage(fileName,n1,n2,n3)
  sf = SimpleFrame.asImagePanels(f)

def gain(hw,f):
  """ normalize RMS amplitude within overlapping windows, half-width hw """
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  return div(f,sqrt(g),f)

def setPlotVars(dw):
  global ne1,nel,se,su,el,ul,xl,yl
  se = dw.getSampling1() # error sampling
  su = dw.getSamplingU() # shift sampling
  ne1 = se.getCount()
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting 
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting
  yl = [s3.getFirst(),s3.getLast()] # frame limits for plotting

def goShiftsSAG(sf,f,sg,g,dw,dg1,ng,find):
  desc = " structure aligned grid"
  label = "sag"
  printConstraints(desc,dg1)
  # Can do flattening separately, comment these steps out, and read the results
  # from disk.
  print "find slopes"
  goSlopes(f,dw) # compute slopes for flattening
  print "Flatten"
  goFlat(f,dw) # flatten
  #
  x1m = readImage("x1",ne1,n2,n3) # flattening mappings
  ff = readImage("ff",ne1,n2,n3) # flattened image
  g1,g1Flat = makeAutomaticGrid(x1m,ff,dw,dg1,ng) # structure aligned grid
  g1 = getG1Floats(x1m,g1Flat,sf.getDelta())
  g2 = getG2(dg2)
  g3 = getG3(dg3)
  print "g2:"; dump(g2)
  print "g3:"; dump(g3)
  return goShifts(sf,f,sg,g,g1,g2,g3,dw,desc,label,find)

def printConstraints(desc,dg1):
  k1min = int(ceil( r1min*dg1))
  k1max = int(floor(r1max*dg1))
  k2min = int(ceil( r2min*dg2))
  k2max = int(floor(r2max*dg2))
  k3min = int(ceil( r3min*dg3))
  k3max = int(floor(r3max*dg3))
  info = """Constraints: %s
  r1min=%g, r1max=%g, dg1=%g, k1min=%g, k1max=%g
  r2min=%g, r2max=%g, dg2=%g, k2min=%g, k2max=%g
  r3min=%g, r3max=%g, dg3=%g, k3min=%g, k3max=%g"""%(desc,r1min,r1max,dg1,k1min,
                                                     k1max,r2min,r2max,dg2,
                                                     k2min,k2max,r3min,r3max,
                                                     dg3,k3min,k3max)
  print info

def goSlopes(f,dw):
  fs = copy(ne1,n2,n3,f) # f with only ne1 samples
  zm = ZeroMask(0.1,1.0,1.0,1.0,fs)
  sigma1 = 8.0
  sigma2 = 8.0
  sigma3 = 8.0
  print "LSF: sigma1=%g, sigma2=%g, sigma3=%g"%(sigma1,sigma2,sigma3)
  pmax = 2.0
  lsf = LocalSlopeFinder(sigma1,sigma2,sigma3,pmax)
  p2 = zerofloat(ne1,n2,n3)
  p3 = zerofloat(ne1,n2,n3)
  ep = zerofloat(ne1,n2,n3)
  lsf.findSlopes(fs,p2,p3,ep)
  zero = 0.00;
  tiny = 0.01;
  zm.apply(zero,p2);
  zm.apply(zero,p3);
  zm.apply(tiny,ep);
  writeImage("p2",p2)
  writeImage("p3",p3)
  writeImage("ep",ep)
  s1f = Sampling(ne1,d1,f1)
  plotUtils.plotPP3(p2,title="P2",s1=s1f,s2=s2,s3=s3,label1="PP time (s)",cbar="Slope",
          cmap1=bwr,slices=slices,he0=he0)
  plotUtils.plotPP3(p3,title="P3",s1=s1f,s2=s2,s3=s3,label1="PP time (s)",cbar="Slope",
        cmap1=bwr,slices=slices,he0=he0)
  plotUtils.plotPP3(ep,title="Planarity",s1=s1f,s2=s2,s3=s3,label1="PP time (s)",
        cbar="Planarity",cmap1=jet,clips1=[0.0,1.0],slices=slices,he0=he0)

def goFlat(f,dw):
  s1f = Sampling(ne1)
  s2f = Sampling(n2)
  s3f = Sampling(n3)
  print "hi"
  p2 = readImage("p2",ne1,n2,n3)
  print "hibye"
  p3 = readImage("p3",ne1,n2,n3)
  ep = readImage("ep",ne1,n2,n3)
  print "hey listen"
  ep = pow(ep,2)
  fl = Flattener3()
  fl.setWeight1(1.0)
  fl.setIterations(0.1,1000)
  fm = fl.getMappingsFromSlopes(s1f,s2f,s3f,p2,p3,ep)
  ff = fm.flatten(f)
  h = fm.unflatten(ff)
  s = fm.getShiftsS()
  x1 = fm.x1
  y1 = fm.u1
  writeImage("x1",x1)
  writeImage("y1",y1)
  writeImage("ff",ff)
  writeImage("fs",s)
  sfs = Sampling(ne1,d1,f1)
  plotUtils.plotPP3(ff,title="PP flat",s2=s2,s3=s3,label1="Tau",slices=slices,he0=he0)
  plotUtils.plotPP3(h,title="PP unflattened",s1=sfs,s2=s2,s3=s3,label1="PP time (s)",
          slices=slices,he0=he0)
  plotUtils.plotPP3(s,title="PP flattening shifts",s1=sfs,s2=s2,s3=s3,slices=slices,
          label1="PP time (s)",cbar="Shift (samples)",cmap1=jet,clips1=None,
          he0=he0)
  print "average shift =",sum(s)/(n1*n2),"samples"

def makeAutomaticGrid(x1m,ff,dw,dg1,ng):
  es = getEnvelopeSum(ff)
  g1Flat = Subsample.subsample(es,dg1,ng)
  ng = len(g1Flat)
  g1 = zeroint(ng,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(ng):
        x = x1m[i3][i2][g1Flat[i1]]
        g1[i3][i2][i1] = int(x+0.5)
  print "g1Flat:"; dump(g1Flat)
  return g1,g1Flat

def getEnvelopeSum(ff):
  htf = HilbertTransformFilter()
  n3 = len(ff)
  n2 = len(ff[0])
  n1 = len(ff[0][0])
  es = zerofloat(n1)
  for i3 in range(n3):
    for i2 in range(n2):
      x = ff[i3][i2]
      y = copy(x)
      htf.apply(n1,x,y)
      for i1 in range(n1):
        es[i1] = es[i1] + sqrt(x[i1]*x[i1]+y[i1]*y[i1])
  return es

def goShifts(sf,f,sg,g,g1,g2,g3,dw,desc,label,find):
  print "findShifts"
  if find:
    u = dw.findShifts(sf,f,sg,g,g1,g2,g3)
    writeImage("u_"+label,u)
  u = readImage("u_"+label,n1,n2,n3)
  checkShifts3(u)
  coordMap = Viewer3P.getSparseCoordsMap(s1,g1,s2,g2,s3,g3)
  plotGrid(f,coordMap,desc)
  plotWarped(sf,f,sg,g,u,desc)
  plotVpvs(f,u,desc)
  plotShifts(f,u,desc)
  return u
 
def checkShifts3(u):
  n3 = len(u)
  n2 = len(u[0])
  results = []
  monotonic = True
  for i3 in range(n3):
    for i2 in range(n2):
      if not isMonotonic(u[i3][i2]):
        monotonic = False
        results.append("i2,i3=%d,%d: Interpolated shifts are not monotonic!"%(
          i2,i3))
  if monotonic:
    results.append("Interpolated shifts are monotonic")
  return results

def plotGrid(pp,coordMap,desc):
  # Plot PP with grid points
  s1 = Sampling(len(pp[0][0]),d1,f1)
  cm = [coordMap,"rO",6.0]
  plotUtils.plotPP3(pp,coordMap=cm,title="PP",s1=s1,s2=s2,s3=s3,clips1=iClips,
          label1="PP time (s)",vInterval1=0.2,cmap1=iMap,cbw=100,slices=slices,
          limits1=el,limits2=xl,limits3=yl,he0=he0)

def plotWarped(sf,f,sg,g,u,desc):
  # Plot Warped PS1 image and the difference between the warped and input.
  h = WarpUtils.applyShifts(sf,u,sg,g)
  plotUtils.plotPP3(h,title="PS1 warped"+desc,s1=s1,s2=s2,s3=s3,clips1=iClips,
          cmap1=iMap,label1="PP time (s)",cbw=100,limits1=el,slices=slices,
          he0=he0)
  writeImage(h,"sg")
  nrms = WarpUtils.computeNrms(ne1,f,h)
  print desc+": nrms=%g"%nrms
  # d = sub(h,f)
  # plotUtils.plotPP3(d,title="PP-PS1 warped"+desc,s1=s1,s2=s2,s3=s3,clips1=iClips,
  #         cmap1=iMap,label1="PP time (s)",cbw=100,limits1=el,slices=slices)

def plotVpvs(f,u,desc,coordMap=None):
  # Plot Vp/Vs ratios
  zm = ZeroMask(f)
  vpvs = WarpUtils.vpvs(su,u)
  print "vpvs: min=%g, max=%g"%(min(vpvs),max(vpvs))
  zm.apply(0.00,vpvs)
  cm = None
  if coordMap:
    cm = [coordMap,"rO",6.0]
  plotUtils.plotPP3(vpvs,coordMap=cm,title="Vp/Vs"+desc,s1=s1,s2=s2,s3=s3,
          clips1=vClips,label1="PP time (s)",cbar="Vp/Vs ratio",cmap1=jet,
          cbw=100,limits1=el,limits2=xl,limits3=yl,slices=slices,he0=he0)

def plotShifts(f,u,desc):
  # Plot 2D shifts
  zm = ZeroMask(f)
  zm.apply(0.00,copy(u))
  plotUtils.plotPP3(u,title="Vertical shifts"+desc,s1=s1,s2=s2,s3=s3,clips1=uClips,
          label1="PS time (s)",cbar="Vertical shift (s)",cmap1=jet,cbw=100,
          limits1=el,slices=slices,he0=he0)

def getG1Floats(x1m,g1Flat,d1):
  ng = len(g1Flat)
  g1 = zerofloat(ng,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(ng):
        g1[i3][i2][i1] = x1m[i3][i2][g1Flat[i1]]*d1
  return g1

def getG2(dg2):
  g2 = Subsample.subsample(n2,dg2); # simple regular interval grid
  ng2 = len(g2)  
  g2f = zerofloat(ng2)
  for i2 in range(ng2):
    g2f[i2] = s2.getValue(g2[i2])
  return g2f

def getG3(dg3):
  g3 = Subsample.subsample(n3,dg3); # simple regular interval grid
  ng3 = len(g3)  
  g3f = zerofloat(ng3)
  for i3 in range(ng3):
    g3f[i3] = s3.getValue(g3[i3])
  return g3f

def toFloats(f):
  n3 = len(f)
  n2 = len(f[0])
  n1 = len(f[0][0])
  ff = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      for i1 in range(n1):
        ff[i3][i2][i1] = float(f[i3][i2][i1])
  return ff

def warpPSToPP(uFileName,gFileName):
  u = getData(uFileName,s1,s2,s3)
  g = getData(gFileName,s1,s2,s3)
  warper = Warper()
  dump(u[100][50])
  print "min u Before = "+str(min(u))
  print "max u Before = "+str(max(u))
  u = div(u,0.002)
  print "min u After = "+str(min(u))
  print "max u After = "+str(max(u))
  #u = add(u,rampfloat(0.0,1.0,0.0,0.0,s1.getCount(),s2.getCount(),s3.getCount()))
  print "min u After2 = "+str(min(u))
  print "max u After2 = "+str(max(u))
  print "u1 = "+str(u[100][50][100])
  sg = warper.applyS(u,g)
  print "min sg = "+str(min(sg))
  print "max sg = "+str(max(sg))
  writeImage("sg",sg)
  
  

















#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

