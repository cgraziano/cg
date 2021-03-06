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

def createAndWriteTensors(tensorFileDir,tensorFileName,datFileDir,datFileName,s1,s2,s3,sigma1,sigma2,sigma3):
  print "Find and write tensors:"
  goTensors(datFileDir,datFileName,s1,s2,s3,tensorFileDir,tensorFileName,sigma1,sigma2,sigma3)
  print "done"

def displayTensorsAndData(tensorFileDir,tensorFileName,datFileDir,datFileName,s1,s2,s3):
  print "Display tensors and data:"
  data = readImageFromDat(datFileDir,datFileName,s1.getCount(),s2.getCount(),s3.getCount())
  display(s1,s2,s3,data,tensorFileDir,tensorFileName)
  print "done"

def smoothDataWithTensorsAndWriteResult(tensorFileDir,tensorFileName,smoothFileDir,smoothFileName,datFileDir,datFileName,s1,s2,s3,smoothingDegree):
  print "Smooth data with tensors and write result:"
  data = readImageFromDat(datFileDir,datFileName,s1.getCount(),s2.getCount(),s3.getCount())
  zm = ZeroMask(data)
  dataSmooth = goSmooth(data,smoothingDegree,tensorFileDir,tensorFileName)
  #WarpUtils.normalize(dataSmooth,-1.5,1.5)#The data didn't look like seismic data after appling this.
  zm.apply(0.0,dataSmooth)
  writeImage(smoothFileDir,smoothFileName,dataSmooth)
  print "done"

def displaySmoothData(datDir,datFileName,smoothFileDir,smoothFileName,s1,s2,s3):
  print "Display smoothed data:"
  n1 = s1.getCount()
  dt1 = s1.getDelta()
  tmax = (n1-1)*dt1
  data = readImageFromDat(datDir,datFileName,s1.getCount(),s2.getCount(),s3.getCount())
  dataSmooth = readImageFromDat(smoothFileDir,smoothFileName,s1.getCount(),s2.getCount(),s3.getCount())
  clipSmooth = max(dataSmooth)*0.1
  clip = max(data)*0.05
  plotUtils.plotPP3(dataSmooth,s1=s1,s2=s2,s3=s3,title=datFileName+" Smooth",label1=datFileName+" time (s)",\
  clips1=[-clipSmooth,clipSmooth],cbw=100,limits1=[0.0,tmax])
  plotUtils.plotPP3(data,s1=s1,s2=s2,s3=s3,title=datFileName,label1=datFileName+" time (s)",\
  clips1=[-clip,clip],cbw=100,limits1=[0.0,tmax])
  #showTwo(s1,s2,s3,data,dataSmooth)
  print "done"


def goTensors(datFileDir,datFileName,s1,s2,s3,tensorFileDir,tensorFileName,sigma1,sigma2,sigma3):

  print "Construct Local Orient Filter"
  lof = LocalOrientFilter(sigma1,sigma2,sigma3)
  print "Apply Local Orient Filter"
  et3 = lof.applyForTensors(readImageFromDat(datFileDir,datFileName,s1.getCount(),s2.getCount(),s3.getCount()))
  print "Invert Structure"
  et3.invertStructure(0.0,2.0,4.0)
  print "Write Tensors"
  writeTensors(tensorFileDir,tensorFileName,et3)

def goSmooth(f,smoothingDegree,tensorFileDir,tensorFileName):
  et3 = readTensors(tensorFileDir,tensorFileName)
  lsf = LocalSmoothingFilter()
  g = copy(f)
  lsf.apply(et3,smoothingDegree,f,g)
  return g

def computeAndWriteSemblance(s1,s2,s3,datFileDir,datFileName,\
  tensorFileDir,tensorFileName,semblanceHalfWidth1,semblanceHalfWidth2,
  semblanceFileDir,semblanceFileName):
  f = readImageFromDat(datFileDir,datFileName,s1.getCount(),s2.getCount(),s3.getCount())
  tensors = readTensors(tensorFileDir,tensorFileName)
  g = copy(f)
  lsf = LocalSemblanceFilter(semblanceHalfWidth1, semblanceHalfWidth2)
  semf = lsf.semblance(LocalSemblanceFilter.Direction3.UVW,tensors,f,g)
  writeImage(semblanceFileDir,semblanceFileName,semf)
  
def display(s1,s2,s3,f,tensorFileDir,tensorFileName):
  world = World()
  s1 = Sampling(s1.getCount(),s1.getDelta()*1000,s1.getFirst())
  ipg = ImagePanelGroup(s1,s2,s3,f)
  clip = max(f)*0.1
  ipg.setClips(-clip,clip)
  world.addChild(ipg)
  # ipg = addImageToWorld(world,f)
  et3 = readTensors(tensorFileDir,tensorFileName)
  plotUtils.addTensorsInImage(s1,s2,s3,ipg.getImagePanel(Axis.X),et3,30)
  plotUtils.addTensorsInImage(s1,s2,s3,ipg.getImagePanel(Axis.Y),et3,30)
  plotUtils.addTensorsInImage(s1,s2,s3,ipg.getImagePanel(Axis.Z),et3,30)
  frame = plotUtils.makeFrame(s1,s2,s3,world)

def showTwo(s1,s2,s3,g1,g2):
  s1 = Sampling(s1.getCount(),s1.getDelta()*1000,s1.getFirst())
  sf = SimpleFrame()
  clip = max(g2)*0.1
  for g in [g1,g2]:
    ipg = sf.addImagePanels(s1,s2,s3,g)
    ipg.setClips(-clip,clip)
  sf.orbitView.setScale(1.0)
  sf.orbitView.setAxesScale(0.75,0.75,1.5)
  sf.setSize(1250,900)
  
def writeTensors(tensorFileDir,tensorFileName,tensors):
  fos = FileOutputStream(tensorFileDir+tensorFileName+".dat")
  oos = ObjectOutputStream(fos)
  oos.writeObject(tensors)
  oos.close()

def readTensors(tensorFileDir,tensorFileName):
  fis = FileInputStream(tensorFileDir+tensorFileName+".dat")
  ois = PythonObjectInputStream(fis)
  tensors = ois.readObject()
  ois.close()
  return tensors


def writeImage(datDir,fileName,x):
  print datDir+fileName+".dat"
  aos = ArrayOutputStream(datDir+fileName+".dat")
  aos.writeFloats(x)
  aos.close()

def readImageFromDat(datDir,fileName,n1,n2,n3):
  print fileName
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datDir+fileName+".dat")
  ais.readFloats(x)
  ais.close()
  return x


###############################################################################
#run(main)

