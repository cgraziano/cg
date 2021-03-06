from imports import *
from linalgebra import ToeplitzRecursion
from edu.mines.jtk.util.ArrayMath import *
from edu.mines.jtk.dsp.Conv import *
from edu.mines.jtk.lapack import DMatrix
from random import *
 
def main(args):
  numtests = 50
  fcount = 0;
  pcount = 0;
  ppr = PlotPanel()
  ppr.setTitle("Plot of 1st row of R")
  ppg = PlotPanel()
  ppg.setTitle("Plot of g")
  ppf = PlotPanel()
  ppf.setTitle("Plot of f")
  for i in range(0,numtests):
    lowV = .01
    highV = 2000
    n = 1000
    r = zerofloat(n)
    f = zerofloat(n)
    #for i in range(0,n):
    #  r[i] = i
    #  f[i] = i
    #r = [1,2]
    #f = [100000000000000,90000000000000]

    #frg = toeplitzSystem(r,f)
    #frgj = ToeplitzRecursion.solveToeplitzSystem(r,f)
    frg = makeRandomToeplitzSystem(n,lowV,highV)
    print getDeterminantToeplitz(frg[1])
    fAct = frg[0]
    r = frg[1]
    g = frg[2]
    x = zerofloat(n)
    xcor(n,0,fAct,n,0,fAct,n,0,x)
    SimplePlot.asPoints(x)
    fCalc = ToeplitzRecursion.solve(r,g)
    diff = sub(fAct,fCalc)
    sqErr = pow(diff,2.0)
    meansqErr = sum(sqErr)/n
    rms = pow(meansqErr,0.5)
    print "rms error = "+str(rms)
    if rms<.1:
      print "test passed, rmsErr<.1"
      pvr = PointsView(r)
      pvr.setLineColor(Color.GREEN)
      pvg = PointsView(g)
      pvg.setLineColor(Color.GREEN)
      pvf = PointsView(fAct)
      pvf.setLineColor(Color.GREEN)
      pcount = pcount+1
    else:
      print "test failed, rmsErr>=.1"
      pvr = PointsView(r)
      pvr.setLineColor(Color.RED)
      pvg = PointsView(g)
      pvg.setLineColor(Color.RED)
      pvf = PointsView(fAct)
      pvf.setLineColor(Color.RED)
      fcount = fcount+1
    ppr.addTiledView(pvr)
    ppg.addTiledView(pvg)
    ppf.addTiledView(pvf)
  print "number failed "+str(fcount)
  print "number passed "+str(pcount)
  pfr = PlotFrame(ppr)
  pfr.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pfr.setVisible(True)
  pfg = PlotFrame(ppg)
  pfg.setVisible(True)
  pfg.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  pff = PlotFrame(ppf)
  pff.setVisible(True)
  pff.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)





#Creates Toeplitz System were
#    Rf = g and R and f are random between a range (a,b)
#    n determines the length of f and g as well as what size
#    r is(n x n)
#    |r0 r1 r2 rn||f0|  |g0|
#    |r1 r0 r1 r2||f1|  |g1|
#    |r2 r1 r0 r1||f2| =|g2|
#    |rn r2 r1 r0||fn|  |gn|
#Output: 2D array where frg[0] is f
#                       frg[1] is r
#                       frg[2] is g
def makeRandomToeplitzSystem(n,a,b):
  r = zerofloat(n)
  f = zerofloat(n)
  g = zerofloat(n)
  frg = zerofloat(n,3)
  for i in range(0,n):
    f[i] = uniform(a,b)
    r[i] = uniform(a,b)

  #quickSort(r)
  #r = reverse(r)

  for i1 in range(0,n):
    for i2 in range(0,n):
      g[i1] += r[abs(i1-i2)]*f[i2]
  frg[0] = f
  frg[1] = r
  frg[2] = g
  return frg

def toeplitzSystem(r,f):
  n = len(r)
  g = zerofloat(n)
  frg = zerofloat(n,3)

  for i1 in range(0,n):
    for i2 in range(0,n):
      g[i1] += r[abs(i1-i2)]*f[i2]
  frg[0] = f
  frg[1] = r
  frg[2] = g
  return frg

def getDeterminantToeplitz(r):
  n = len(r)
  d = DMatrix(n,n)
  for i1 in range(0,n):
    for i2 in range(0,n):
      d.set(i1,i2,r[abs(i1-i2)])
  return d.det()



#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
