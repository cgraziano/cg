#############################################################################
# Plotting Software

from imports import *
from edu.mines.jtk.sgl import *

############################################################################

#vmin and vmax are in units defined by the sampling.
#hsize and vsize are in pixels
#Default hsize and vsize 960,560
def plotTracesSideBySide(st, traces, 
  vlabel=None, vminmax=None, vint=None, 
  hlabel=None, hminmax=None, hint=None,
  hsize=None, vsize=None,
  tilespacing=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):

  ntraces = len(traces)
  pvs = list()
  #Build list of PointViews
  for itrace in range(0,ntraces):
    pvs.append(PointsView(st,traces[itrace]))
  for itrace in range(0,ntraces):
    pvs[itrace].setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  pp = PlotPanel(1,ntraces,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.LEFT_TOP)
  if tilespacing:
    mosaic = pp.getMosaic()
    mosaic.setWidthTileSpacing(tilespacing)
  if vlabel:
    pp.setVLabel(vlabel)
  if vminmax:
    pp.setVLimits(vminmax[0],vminmax[1])
  if vint:
    pp.setVInterval(vint)
  for itrace in range(0,ntraces):
    if hlabel:
      pp.setHLabel(itrace,hlabel[itrace])
    if hminmax:
      pp.setHLimits(itrace,hminmax[itrace][0],hminmax[itrace][1])
    if hint[itrace]:
      pp.setHInterval(itrace,hint[itrace])
    pp.addTiledView(0,itrace,pvs[itrace])
  pf = PlotFrame(pp)
  if hsize and vsize:
    pf.setSize(hsize,vsize)
  else:
    pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

#vmin and vmax are in units defined by the sampling.
#hsize and vsize are in pixels
#Default hsize and vsize 960,560
def plot2TracesInSamePlotSideBySideWithOtherPlots(st, tracePairs, 
  vlabel=None, vminmax=None, vint=None, 
  hlabel=None, hminmax=None, hint=None,
  hsize=None, vsize=None,
  color=None,
  tilespacing=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):
  nTracePairs = len(tracePairs)
  pvs1 = list()
  pvs2 = list()
  #Build list of PointViews
  print "hint"+str(hint)
  for itrace in range(0,nTracePairs):
    pvs1.append(PointsView(st,tracePairs[itrace][0]))
    pvs2.append(PointsView(st,tracePairs[itrace][1]))
  for itrace in range(0,nTracePairs):
    pvs1[itrace].setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
    pvs2[itrace].setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
    pvs1[itrace].setLineWidth(2.0)
    pvs1[itrace].setLineColor(color[0])
    pvs2[itrace].setLineColor(color[1])
    #pvs2[itrace].setLineWidth(1.5)
    #pvs2[itrace].setLineStyle(PointsView.Line.DASH)
  pp = PlotPanel(1,nTracePairs,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.LEFT_TOP)
  if tilespacing:
    mosaic = pp.getMosaic()
    mosaic.setWidthTileSpacing(tilespacing)
  if vlabel:
    pp.setVLabel(vlabel)
  if vminmax:
    pp.setVLimits(vminmax[0],vminmax[1])
  if vint:
    pp.setVInterval(vint)
  for itrace in range(0,nTracePairs):
    if hlabel:
      pp.setHLabel(itrace,hlabel[itrace])
    if hminmax:
      pp.setHLimits(itrace,hminmax[itrace][0],hminmax[itrace][1])
    if hint[itrace]:
      pp.setHInterval(itrace,hint[itrace])
    pp.addTiledView(0,itrace,pvs1[itrace])
    pp.addTiledView(0,itrace,pvs2[itrace])
  pf = PlotFrame(pp)
  if hsize and vsize:
    pf.setSize(hsize,vsize)
  else:
    pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)


def plotMeasOnTopOfEachOther(sm, measures, 
  color=None,
  vlabel=None, vminmax=None, vint=None, 
  hlabel=None, hminmax=None, hint=None,
  hsize=None, vsize=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):

  nmeasures = len(measures)
  pvs = list()
  #Build list of PointViews
  for imeasure in range(0,nmeasures):
    pvs.append(PointsView(sm,measures[imeasure]))
  for imeasure in range(0,nmeasures):
    pvs[imeasure].setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
    if color:
      pvs[imeasure].setLineColor(color[imeasure])
  pp = PlotPanel(nmeasures,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_BOTTOM)
  if hlabel:
    pp.setHLabel(hlabel)
  if hminmax:
    pp.setHLimits(hminmax[0],hminmax[1])
  if hint:
    pp.setHInterval(hint)
  for imeasure in range(0,nmeasures):
    pp.addTiledView(imeasure,0,pvs[imeasure])
    print "imeasure = "+str(imeasure)
    if (vminmax[imeasure]):
      pp.setVLimits(imeasure,vminmax[imeasure][0],vminmax[imeasure][1])
    if vint[imeasure]:
      pp.setVInterval(imeasure,vint[imeasure])
      print "vint = "+str(imeasure)
    if vlabel:
      pp.setVLabel(imeasure,vlabel[imeasure])
      print "vlabel = "+str(imeasure)
  pf = PlotFrame(pp)
  if hsize and vsize:
    pf.setSize(hsize,vsize)
  else:
    pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plotMeasInSamePlot(sm, measures, 
  color=None,
  vlabel=None, vminmax=None, vint=None, 
  hlabel=None, hminmax=None, hint=None,
  hsize=None, vsize=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):

  nmeasures = len(measures)
  pvs = list()
  #Build list of PointViews
  for imeasure in range(0,nmeasures):
    pvs.append(PointsView(sm,measures[imeasure]))
  for imeasure in range(0,nmeasures):
    pvs[imeasure].setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
    pvs[imeasure].setLineWidth(1.0);
    pvs[imeasure].setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pvs[imeasure].setMarkSize(4.0);
    if color:
      pvs[imeasure].setLineColor(color[imeasure])
      pvs[imeasure].setMarkColor(color[imeasure])
  pp = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_BOTTOM)
  if hlabel:
    pp.setHLabel(hlabel)
  if hminmax:
    pp.setHLimits(hminmax[0],hminmax[1])
  if hint:
    pp.setHInterval(hint)
  if (vminmax):
    pp.setVLimits(vminmax[0],vminmax[1])
  if vint:
    pp.setVInterval(vint)
  if vlabel:
    pp.setVLabel(vlabel)
  for imeasure in range(0,nmeasures):
    pp.addTiledView(0,0,pvs[imeasure])
    print "imeasure = "+str(imeasure)
  pf = PlotFrame(pp)
  if hsize and vsize:
    pf.setSize(hsize,vsize)
  else:
    pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plotMeasInSamePlotNoise(sm, measures, 
  color=None,
  vlabel=None, vminmax=None, vint=None, 
  hlabel=None, hminmax=None, hint=None,
  hsize=None, vsize=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):

  nmeasures = len(measures)
  pvs = list()
  #Build list of PointViews
  for imeasure in range(0,nmeasures):
    pvs.append(PointsView(sm,measures[imeasure]))
  for imeasure in range(0,nmeasures):
    pvs[imeasure].setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
    pvs[imeasure].setLineStyle(PointsView.Line.NONE)
    pvs[imeasure].setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pvs[imeasure].setMarkSize(5.0);
    if color:
      pvs[imeasure].setLineColor(color[imeasure])
      pvs[imeasure].setMarkColor(color[imeasure])
  pp = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_BOTTOM)
  if hlabel:
    pp.setHLabel(hlabel)
  if hminmax:
    pp.setHLimits(hminmax[0],hminmax[1])
  if hint:
    pp.setHInterval(hint)
  if (vminmax):
    pp.setVLimits(vminmax[0],vminmax[1])
  if vint:
    pp.setVInterval(vint)
  if vlabel:
    pp.setVLabel(vlabel)
  for imeasure in range(0,nmeasures):
    pp.addTiledView(0,0,pvs[imeasure])
    print "imeasure = "+str(imeasure)
  pf = PlotFrame(pp)
  if hsize and vsize:
    pf.setSize(hsize,vsize)
  else:
    pf.setSize(960,360)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)


def plotMeasInSamePlotLargeMarks(sm, measures, 
  color=None,
  vlabel=None, vminmax=None, vint=None, 
  hlabel=None, hminmax=None, hint=None,
  hsize=None, vsize=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):

  nmeasures = len(measures)
  pvs = list()
  #Build list of PointViews
  for imeasure in range(0,nmeasures):
    pvs.append(PointsView(sm,measures[imeasure]))
  for imeasure in range(0,nmeasures):
    pvs[imeasure].setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
    pvs[imeasure].setLineWidth(1.0);
    pvs[imeasure].setMarkStyle(PointsView.Mark.FILLED_CIRCLE);
    pvs[imeasure].setMarkSize(10.0);
    if color:
      pvs[imeasure].setLineColor(color[imeasure])
      pvs[imeasure].setMarkColor(color[imeasure])
  pp = PlotPanel(1,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_BOTTOM)
  if hlabel:
    pp.setHLabel(hlabel)
  if hminmax:
    pp.setHLimits(hminmax[0],hminmax[1])
  if hint:
    pp.setHInterval(hint)
  if (vminmax):
    pp.setVLimits(vminmax[0],vminmax[1])
  if vint:
    pp.setVInterval(vint)
  if vlabel:
    pp.setVLabel(vlabel)
  for imeasure in range(0,nmeasures):
    pp.addTiledView(0,0,pvs[imeasure])
    print "imeasure = "+str(imeasure)
  pf = PlotFrame(pp)
  if hsize and vsize:
    pf.setSize(hsize,vsize)
  else:
    pf.setSize(960,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)



#Default hsize and vsize 720,400
def plotWavelets(st, wavelets, 
  hint=None,
  linestyle=None,linecolor=None,markstyle=None,
  dotsonsticks=None,
  hsize=None, vsize=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):
  sp = SimplePlot()
  if linestyle:
    ls = linestyle
  else:
    ls = [PointsView.Line.SOLID,PointsView.Line.SOLID,PointsView.Line.SOLID]
  if linecolor:
    lc = linecolor
  else:
    lc = [Color.BLUE,Color.RED,Color.GREEN]
  if markstyle:
    ms = markstyle
  else:
    ms = [PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.NONE,PointsView.Mark.FILLED_CIRCLE]
  nc = len(wavelets)
  hsmax = 1.0
  for ih in range(nc):
    if ih==0:
      if dotsonsticks:
        pv = sp.addSequence(st,wavelets[ih])
      else:
        pv = sp.addPoints(st,wavelets[ih])
        pv.setLineStyle(ls[ih])
        pv.setLineColor(lc[ih])
        pv.setMarkStyle(ms[ih])
        pv.setMarkSize(3.4)
      hsmax = max(hsmax,abs(max(wavelets[ih])),abs(min(wavelets[ih])))
    if ih==1:
      if dotsonsticks:
        pv = sp.addSequence(st,wavelets[ih])
      else:
        pv = sp.addPoints(st,wavelets[ih])
        pv.setLineStyle(ls[ih])
        pv.setLineColor(lc[ih])
        pv.setLineWidth(1)
      hsmax = max(hsmax,abs(max(wavelets[ih])),abs(min(wavelets[ih])))
    if ih==2:
      if dotsonsticks:
        pv = sp.addSequence(st,wavelets[ih])
      else:
        pv = sp.addPoints(st,wavelets[ih])
        pv.setLineStyle(ls[ih])
        pv.setLineColor(lc[ih])
        pv.setMarkStyle(ms[ih])
        pv.setMarkSize(3.4)
      hsmax = max(hsmax,abs(max(wavelets[ih])),abs(min(wavelets[ih])))
  hmax = hsmax*1.05
  sp.setVLimits(-1.1,1.1)
  if hint:
    sp.setHInterval(hint)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Amplitude (normalized)")
  if hsize and vsize:
    sp.setSize(hsize,vsize)
  else:
    sp.setSize(960,560)
  if title:
    if pngDir==None:
      sp.setTitle(title)
  if pngDir:
    if slide:
      print "zz = "+str(aspectRatio)
      sp.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      sp.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        sp.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        sp.paintToPng(720.0,3.08,pngDir)
      if twocol:
        sp.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        sp.paintToPng(720.0,6.51,pngDir)

def plotWaveletsNoise(st, wavelets, 
  hint=None,
  linestyle=None,linecolor=None,markstyle=None,
  dotsonsticks=None,
  hsize=None, vsize=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):
  sp = SimplePlot()
  if linestyle:
    ls = linestyle
  else:
    ls = [PointsView.Line.SOLID,PointsView.Line.NONE,PointsView.Line.NONE,PointsView.Line.NONE,PointsView.Line.NONE,PointsView.Line.NONE,PointsView.Line.NONE]
  if linecolor:
    lc = linecolor
  else:
    lc = [Color.RED,Color.BLACK,Color.GREEN,Color.BLUE,Color.CYAN,Color.PINK,Color.ORANGE]
  if markstyle:
    ms = markstyle
  else:
    ms = [PointsView.Mark.NONE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE,PointsView.Mark.FILLED_CIRCLE]
  nc = len(wavelets)
  hsmax = 1.0
  for ih in range(nc):
    print "ih = "+str(ih)
    pv = sp.addPoints(st,wavelets[ih])
    pv.setLineStyle(ls[ih])
    if ih==0:
      pv.setLineColor(lc[ih])
      pv.setLineWidth(3.0)
    else:
      pv.setMarkColor(lc[ih])
    pv.setMarkStyle(ms[ih])
    pv.setMarkSize(3.4)
    hsmax = max(hsmax,abs(max(wavelets[ih])),abs(min(wavelets[ih])))
  hmax = hsmax*1.05
  sp.setVLimits(-1.1,1.1)
  if hint:
    sp.setHInterval(hint)
  sp.setHLabel("Time (s)")
  sp.setVLabel("Amplitude (normalized)")
  if hsize and vsize:
    sp.setSize(hsize,vsize)
  else:
    sp.setSize(960,560)
  if title:
    if pngDir==None:
      sp.setTitle(title)
  if pngDir:
    if slide:
      print "zz = "+str(aspectRatio)
      sp.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      sp.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        sp.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        sp.paintToPng(720.0,3.08,pngDir)
      if twocol:
        sp.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        sp.paintToPng(720.0,6.51,pngDir)


#Example of parameters to set
"""
pngDir = "./report15/synthetics/"
title= "testsino"
vmin,vmax = tmin*dt,tmax*dt
hmin,hmax = None,None
vlabel,vminmax,vint = "time (s)",[vmin,vmax],1.0
hlabel,hminmax,hint = ["f","csbg","hsg","sg"],None,0.015
clip = 2
"""
def plotImagesSideBySide(st, sx, images, 
  vlabel=None, vminmax=None, vint=None, 
  hlabel=None, hminmax=None, hint=None,
  tilespacing=None,
  hsize=None, vsize=None,
  clip=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):

  nimages = len(images)
  pvs = list()
  #Build list of PixelsViews
  for itrace in range(0,nimages):
    pvs.append(PixelsView(st,sx,images[itrace]))
  for itrace in range(0,nimages):
    pvs[itrace].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)
    if clip:
      pvs[itrace].setClips(-clip,clip)
  pp = PlotPanel(1,nimages,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.LEFT_TOP)
  if tilespacing:
    mosaic = pp.getMosaic()
    mosaic.setWidthTileSpacing(tilespacing)
  if vlabel:
    pp.setVLabel(vlabel)
  if vminmax:
    pp.setVLimits(vminmax[0],vminmax[1])
  if vint:
    pp.setVInterval(vint)
  for itrace in range(0,nimages):
    if hlabel:
      pp.setHLabel(itrace,hlabel[itrace])
    if hminmax:
      pp.setHLimits(itrace,hminmax[itrace][0],hminmax[itrace][1])
    if hint:
      pp.setHInterval(itrace,hint)
    pp.addTiledView(0,itrace,pvs[itrace])
  pf = PlotFrame(pp)
  if hsize and vsize:
    pf.setSize(hsize,vsize)
  else:
    pf.setSize(740,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(10.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)

def plotImagesSideBySideDU(st, sx, images, 
  vlabel=None, vminmax=None, vint=None, 
  hlabel=None, hminmax=None, hint=None,
  clipmin=None,clipmax=None,
  tilespacing=None,
  hsize=None, vsize=None,
  clip=None,
  title=None, pngDir=None,
  slide=None, fracWidth=None, fracHeight=None, aspectRatio=None, 
  paper=None, onecol=None, twocol=None):

  nimages = len(images)
  pvs = list()
  #Build list of PixelsViews
  for itrace in range(0,nimages):
    pvs.append(PixelsView(st,sx,images[itrace]))
  for itrace in range(0,nimages):
    pvs[itrace].setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT)

    pvs[itrace].setColorModel(ColorMap.JET)
    pvs[itrace].setClips(clipmin,clipmax)
    if clip:
      pvs[itrace].setClips(-clip,clip)
  pp = PlotPanel(1,nimages,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.LEFT_TOP)
  pp.addColorBar()
  if tilespacing:
    mosaic = pp.getMosaic()
    mosaic.setWidthTileSpacing(tilespacing)
  if vlabel:
    pp.setVLabel(vlabel)
  if vminmax:
    pp.setVLimits(vminmax[0],vminmax[1])
  if vint:
    pp.setVInterval(vint)
  for itrace in range(0,nimages):
    if hlabel:
      pp.setHLabel(itrace,hlabel[itrace])
    if hminmax:
      pp.setHLimits(itrace,hminmax[itrace][0],hminmax[itrace][1])
    if hint:
      pp.setHInterval(itrace,hint)
    pp.addTiledView(0,itrace,pvs[itrace])
  pf = PlotFrame(pp)
  if hsize and vsize:
    pf.setSize(hsize,vsize)
  else:
    pf.setSize(740,560)
  if title:
    if pngDir==None:
      pp.setTitle(title)
  if pngDir:
    if slide:
      pf.setFontSizeForSlide(fracWidth,fracHeight,aspectRatio)
      pngDir = pngDir+title+"w"+str(fracWidth)+"h"+str(fracHeight)+"slide.png"
      pf.paintToPng(720.0,3.0,pngDir)
    if paper:
      if onecol:
        pf.setFontSizeForPrint(8.0,222.0)
        pngDir = pngDir+title+"paper"+"onecol.png"
        pf.paintToPng(720.0,3.08,pngDir)
      if twocol:
        pf.setFontSizeForPrint(10.0,469.0)
        pngDir = pngDir+title+"paper"+"twocol.png"
        pf.paintToPng(720.0,6.51,pngDir)
  pf.setVisible(True)
  pf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)








