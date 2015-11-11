from imports import *
from warp import *
from dnp import *
from viewer import *
from wwarp import *
import tensors
import plotUtils

#In a location of your choosing, create folders called "dat", "sgy", "tensors", "smooth", and "writtenSgy".
#Put the .sgy file in the sgy folder.
#Below, specify the location of these folders below in foldersDir.
# Ex) foldersDir = "C:/data/BellCreek/"
#The folders "dat", "segy", "tensors", "smooth", and "writtenSgy" are located in the BellCreek folder.

foldersDir = "C:/data/BellCreek/"
datDir = foldersDir+"dat/"
sgyDir = foldersDir+"sgy/"
writtenSgyDir = foldersDir+"writtenSgy/"
tensorDir = foldersDir+"tensors/"
smoothDir = foldersDir+"smooth/"
semblanceDir = foldersDir+"semblance/"

#Specify the name of the PP and PS .sgy files.
sgyFileName = "3D_BASE_FF_AVOPSTMSTK_Unfiltered_0pt7Thru1pt5"

#Specify byte locations of the .sgy files that will be used.
#inlineByteLocation = 9
#crosslineByteLocation = 13
inlineByteLocation = 17#Subset
crosslineByteLocation = 25#Subset

#Halfwidth of windows in each dimension used to fine structure of image.
sigma1,sigma2,sigma3 = 8.0,8.0,8.0

#Degree of smoothing in structure tensor's direction.
smoothingDegree = 20.0

#Semblance
semblanceHalfWidth1 = 2
semblanceHalfWidth2 = 2

#The PP and PS data (everything except all sgy headers) will be written to 
#separate .dat files. These files will have the same name as the .sgy files
#specified above.
datFileName = sgyFileName

#For testing purposes, a new .dat file will be created from 
#the written .sgy file.
writtenSgyFileName = sgyFileName+"semblance_"+"semblanceHW1"+str(semblanceHalfWidth1)+"semblanceHW2"+str(semblanceHalfWidth2)+"UVM"
#writtenDatFileName = sgyFileName+"_Written"

#Tensor, smooth, and semblance files will be written because there is not enough RAM to hold all arrays needed for processing.
tensorFileName = sgyFileName+"_tensor"
smoothFileName = sgyFileName+"_smooth"
semblanceFileName = sgyFileName+"_semblance"

#Plotting parameters
bwr = ColorMap.HUE
gry = ColorMap.GRAY
jet = ColorMap.JET
iMap = gry

k1,k2,k3 = 345,20,20 # default indices for plotting
slices = [k1,k2,k3]
he0 = 332 # height of time slice panel in plots
iClips = [-1.5,1.5]
uClips = [0.0,0.75]
vClips = [1.9,3.5]

def main(args):
  srwc = SegyReaderWriterAndConverter(sgyDir,sgyFileName,inlineByteLocation,crosslineByteLocation,datDir,datFileName)
  s1 = srwc.getSampling1()
  s2 = srwc.getSampling2()
  s3 = srwc.getSampling3()

  ##Print header information
  srwc.printHeaderInformation()

  ##Convert .sgy files to .dat files.
  convertSgyFileToDat(srwc)

  ##Display PP and PS images
  #displayImageAccordingToXAndY(srwc)
  #displayImageAccordingToInlineAndCrossline(srwc)
  #displayImageIn3DView(srwc)
  
  ##Create Tensors
  #tensors.createAndWriteTensors(tensorDir,tensorFileName,datDir,datFileName,s1,s2,s3,sigma1,sigma2,sigma3)

  ##Smooth data 
  #tensors.smoothDataWithTensorsAndWriteResult(tensorDir,tensorFileName,smoothDir,smoothFileName,datDir,datFileName,s1,s2,s3,smoothingDegree)

  ##Display tensors and data
  #tensors.displayTensorsAndData(tensorDir,tensorFileName,datDir,datFileName,s1,s2,s3)
  
  ##Display smooth data
  #tensors.displaySmoothData(datDir,datFileName,smoothDir,smoothFileName,s1,s2,s3)

  ##computeSemblance
  tensors.computeAndWriteSemblance(s1,s2,s3,datDir,datFileName,tensorDir,tensorFileName,semblanceHalfWidth1,semblanceHalfWidth2,semblanceDir,semblanceFileName)

  ##Display semblance data
  tensors.displaySmoothData(datDir,datFileName,semblanceDir,semblanceFileName,s1,s2,s3)

  ##What .dat file do you want converted to sgy?
  whichDir = semblanceDir
  whichDatFileName = semblanceFileName

  #srwc.setDatDirAndFileName(whichDir,whichDatFileName)

  #srwc.convertDatToSgy(writtenSgyDir,writtenSgyFileName)



def readImageFromDat(fileDir,fileName,n1,n2,n3):
  print fileName
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileDir+fileName+".dat")
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileDir,fileName,x):
  aos = ArrayOutputStream(fileDir+fileName+".dat")
  aos.writeFloats(x)
  aos.close()


def gain(hw,f):
  """ normalize RMS amplitude within overlapping windows, half-width hw """
  #add epsilon value to eliminate any NaN.
  f = add(eps,f)
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  return div(f,sqrt(g),f)

def setPlotVars(dw,s2,s3):
  global ne1,nel,se,su,el,ul,xl,yl
  se = dw.getSampling1() # error sampling
  su = dw.getSamplingU() # shift sampling
  ne1 = se.getCount()
  el = [se.getFirst(),se.getLast()] # error limits for plotting
  ul = [su.getFirst(),su.getLast()] # shift limits for plotting 
  xl = [s2.getFirst(),s2.getLast()] # trace limits for plotting
  yl = [s3.getFirst(),s3.getLast()] # frame limits for plotting


def plot23(si):
  i2 = si.getI2sAsFloats()
  i3 = si.getI3sAsFloats()
  sp = SimplePlot()
  sp.setHLabel("inline sample index i2")
  sp.setVLabel("crossline sample index i3")
  pv = sp.addPoints(i2,i3)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(i2,i3)
  sp.setSize(w,h)

def plotXY(si):
  x = si.getXs()
  y = si.getYs()
  sp = SimplePlot()
  sp.setHLabel("x coordinate (km)")
  sp.setVLabel("y coordinate (km)")
  pv = sp.addPoints(x,y)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(x,y)
  sp.setSize(w,h)

def goodWidthHeight(x,y):
  xmin,xmax = min(x),max(x)
  ymin,ymax = min(y),max(y)
  w,h = 1000,1000
  if (xmax-xmin)>(ymax-ymin):
    h = int(h*(ymax-ymin)/(xmax-xmin))
  else:
    w = int(w*(xmax-xmin)/(ymax-ymin))
  return w,h

def show3d(f,clip=None):
  print "show3d: f min =",min(f)," max =",max(f)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(f)
  if clip:
    ipg.setClips(-clip,clip)
  frame.orbitView.setScale(2.0)
  frame.setSize(1000,1000)




#Uses the PP and PS SegyReaderWriterAndConverter to convert the PP and PS .sgy files to .dat files.
def convertSgyFileToDat(srwc):
  srwc.convertSgyToDat()

def displayImageAccordingToXAndY(srwc):
  print "Plot according to x and y."
  srwc.plotDataAccordingToXY()

def displayImageAccordingToInlineAndCrossline(srwc):
  print "Plot according to inline and crossline"
  srwc.plotDataAccordingToInlineCrossline()

def displayImageIn3DView(srwc):
  print "display PP and PS in 3D view"
  srwc.show3D()




  
  

















#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())


