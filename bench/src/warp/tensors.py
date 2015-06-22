###############################################################################
# Structure Tensors and Structure Oriented Smoothing for GBC

#from gbcUtils import *
from imports import *
from warp import *
import plotUtils
###############################################################################

#baseDir = getBaseDir()
#s1,s2,s3 = getSamplings()
#iClips = [-1.0,1.0]

#sigma1,sigma2,sigma3 = 8.0,8.0,8.0
#c = 4.0
###############################################################################

#def main(args):
  # goPP()
  # goPS1()
  #goPS2()

def createTensorsAndSmooth(s1,s2,s3,data,fileName,dataDir,smoothingDegree,sigma1,sigma2,sigma3):
  zm = ZeroMask(data)
  print "Begin Tensor Construction"
  goTensors(data,fileName+"_tensors",dataDir,sigma1,sigma2,sigma3)
  print "Finished Tensor Construction"
  display(s1,s2,s3,data,fileName+"_tensors",dataDir)
  print "Begin Smoothing"
  dataSmooth = goSmooth(data,smoothingDegree,fileName+"_tensors",dataDir)
  print "End Smoothing"
  WarpUtils.normalize(dataSmooth,-1.5,1.5)
  zm.apply(0.0,dataSmooth)
  plotUtils.plotPP3(data,s1=s1,s2=s2,s3=s3,title=fileName,label1=fileName+" time (s)",clips1=[-1.0,1.0],
          cbw=100,limits1=[0.0,3.0])
  plotUtils.plotPP3(dataSmooth,s1=s1,s2=s2,s3=s3,title=fileName+" Smooth",label1=fileName+" time (s)",
          clips1=[-1.0,1.0],cbw=100,limits1=[0.0,3.0])
  showTwo(s1,s2,s3,data,dataSmooth)
  return dataSmooth

def goTensors(f,fileName,dataDir,sigma1,sigma2,sigma3):
  print "Construct Local Orient Filter"
  lof = LocalOrientFilter(sigma1,sigma2,sigma3)
  print "Apply Local Orient Filter"
  et3 = lof.applyForTensors(f)
  et3.invertStructure(0.0,2.0,4.0)
  print "Write Tensors"
  writeTensors(et3,fileName,dataDir)

def goSmooth(f,smoothingDegree,fileName,dataDir):
  et3 = readTensors(fileName,dataDir)
  lsf = LocalSmoothingFilter()
  g = copy(f)
  lsf.apply(smoothingDegree,f,g)
  return g
  
def display(s1,s2,s3,f,fileName,dataDir):
  world = World()
  ipg = ImagePanelGroup(s1,s2,s3,f)
  ipg.setClips(-1.0,1.0)
  world.addChild(ipg)
  # ipg = addImageToWorld(world,f)
  et3 = readTensors(fileName,dataDir)
  plotUtils.addTensorsInImage(s1,s2,s3,ipg.getImagePanel(Axis.X),et3,30)
  plotUtils.addTensorsInImage(s1,s2,s3,ipg.getImagePanel(Axis.Y),et3,30)
  plotUtils.addTensorsInImage(s1,s2,s3,ipg.getImagePanel(Axis.Z),et3,30)
  frame = plotUtils.makeFrame(s1,s2,s3,world)

def showTwo(s1,s2,s3,g1,g2):
  sf = SimpleFrame()
  for g in [g1,g2]:
    ipg = sf.addImagePanels(s1,s2,s3,g)
    ipg.setClips(-1.0,1.0)
  sf.orbitView.setScale(1.0)
  sf.orbitView.setAxesScale(0.75,0.75,1.5)
  sf.setSize(1250,900)
  
def writeTensors(tensors,fileName,dataDir):
  fos = FileOutputStream(dataDir+fileName+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  oos.close()

def readTensors(fileName,dataDir):
  fis = FileInputStream(dataDir+fileName+".dat")
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  ois.close()
  return tensors

###############################################################################
#run(main)

