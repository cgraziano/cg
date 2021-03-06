from imports import *

datadir = "C:/Users/Chris/Documents/CWP/Research/research/gbc/dat"

def main(args):
  print "Available Options: pp, ps1, ps2, shiftl, shiftc"
  view("pp")

def view(data):
  if data == "pp":
    f = readImage(datadir+"/pp.dat",2000,150,145)
  elif data == "ps1":
    f = readImage(datadir+"/ps1.dat",2000,150,145)
  elif data == "ps2":
    f = readImage(datadir+"/ps2.dat",2000,150,145)
  elif data == "shiftl":
    f = readImage(datadir+"/u_sag_linear.dat",2000,150,145)
  elif data == "shiftc":
    f = readImage(datadir+"/u_sag_monotonic.dat",2000,150,145)
  sf = SimpleFrame.asImagePanels(f)

def readImage(fileName,n1,n2,n3=1):
  if n3==1:
    x = zerofloat(n1,n2)
  else:
    x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
