package warp;

import java.io.*;
import java.nio.*;
import java.util.Arrays;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.sgl.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class SegyReaderWriterAndConverter {

  

  //Do not include extensions in datFileName or sgyFileName.
  //Reads number of time samples and minimum and maximum inline and crossline from trace headers. 
  //Reads the number of samples from headers.
  //Samples scaled by 0.0001
  /**
   * A converter that converts a .sgy file to a .dat files by stripping the headers of the .sgy file
   * and only writing the data to the .dat file. This converter will also convert the .dat file 
   * to a .sgy file by using the data in the .dat file and the header from a specified .sgy file to
   * write the new .sgy file. Do not include extensions the any file names.
   * @param sgyDir the directory of a .sgy file that will be assumed to be the .sgy file that will need to 
                   be converted to a .dat file. This directory will not be the location of a newly written .sgy 
                   file. 
   * @param sgyFileName the file name of a .sgy file that will be assumed to be the .sgy file that will need to 
                   be converted to a .dat file. The header of this .sgy file will be included in any written
                   .sgy files. The name of the new .sgy file can be changed.
   * @param inlineByteLocation the inline byte location in the trace header of the .sgy file specified above.
   * @param crosslineByteLocation the inline byte location in the trace header of the .sgy file specified above.
   * @param datDir the directory of the new .dat file.
   * @param datFileName the file name of the new .dat file.
   */
  public SegyReaderWriterAndConverter(String sgyDir, String sgyFileName, 
                                      int inlineByteLocation, int crosslineByteLocation,
                                      String datDir, String datFileName)
  {
    _sgyDir = sgyDir;
    _sgyFileName = sgyFileName;
    _inlineByteLocation = inlineByteLocation;
    _crosslineByteLocation = crosslineByteLocation;
    _datDir = datDir;
    _datFileName = datFileName;
    _si = new SegyImage(_sgyDir+_sgyFileName+".sgy");
    _si.setInlineXlineBytes(_inlineByteLocation,_crosslineByteLocation);
    _n1 = _si.getN1();
    _d1 = _si.getD1();
    _i1min = 0;
    _i1max = _n1-1;
    _i2min = _si.getI2Min();
    _i2max = _si.getI2Max();
    _i3min= _si.getI3Min();
    _i3max = _si.getI3Max();
    _n2 = 1+_i2max-_i2min;
    _n3 = 1+_i3max-_i3min;
  }

  /**
   * Changes the directory and file name of the .sgy file to be written and the .sgy file, whose headers will
   * be used for any written .sgy files. Do not include extensions the any file names.
   * @param sgyDir the directory of a .sgy file that will be assumed to be the .sgy file that will need to 
                   be converted to a .dat file. This directory will not be the location of a newly written .sgy 
                   file. 
   * @param sgyFileName the file name of a .sgy file that will be assumed to be the .sgy file that will need to 
                   be converted to a .dat file. The header of this .sgy file will be included in any written
                   .sgy files. 
   * @param inlineByteLocation the inline byte location in the trace header of the .sgy file specified above.
   * @param crosslineByteLocation the inline byte location in the trace header of the .sgy file specified above.

   */
  public void setSgyDirAndFileName(String sgyDir, String sgyFileName, int inlineByteLocation, int crosslineByteLocation) 
  {
    _sgyDir = sgyDir;
    _sgyFileName = sgyFileName;
    _inlineByteLocation = inlineByteLocation;
    _crosslineByteLocation = crosslineByteLocation;
    _si = new SegyImage(_sgyDir+_sgyFileName+".sgy");
    _si.setInlineXlineBytes(_inlineByteLocation,_crosslineByteLocation);
    _n1 = _si.getN1();
    _d1 = _si.getD1();
    _i1min = 0;
    _i1max = _n1-1;
    _i2min = _si.getI2Min();
    _i2max = _si.getI2Max();
    _i3min= _si.getI3Min();
    _i3max = _si.getI3Max();
    _n2 = 1+_i2max-_i2min;
    _n3 = 1+_i3max-_i3min;
  }

  /**
   * Changes the directory and file name of the .dat file to be written.
   * Do not include extensions the any file names.
   * @param datDir the directory of the new .dat file.
   * @param datFileName the file name of the new .dat file.
   */
  public void setDatDirAndFileName(String datDir, String datFileName) {
    _datDir = datDir;
    _datFileName = datFileName;
  }
  
  /**
   * Prints the header information of the current .sgy file set.
   */
  public void printHeaderInformation() 
  {
    System.out.println("Printing Header Information");
    _si.printAllInfo();
  }

  /** 
   * Converts the .sgy file to a .dat file using the set 
   * directories and file names of the .dat and .sgy files.
   */
  public void convertSgyToDat() 
  {
    String datFile = _datDir+_datFileName+".dat";
    double sampleFactorIndex = 1.0;
    _si.writeFloats(datFile,sampleFactorIndex,_i1min,_i1max,_i2min,_i2max,_i3min,_i3max);
  }


  /** 
   * Converts the set .dat file to the .sgy file specified.
   * The output directory will be the 
   * @param writtenSgyDir the directory of the new .sgy file.
   * @param writtenSgyFileName the name of the new .sgy file.
   */
  public void convertDatToSgy(String writtenSgyDir, String writtenSgyFileName)
  {
    System.out.println("Convert .dat file to .sgy file: ");
    float[][][] f = convertDatToFloatArray(_datDir,_datFileName);
    try {
      ArrayFile afNewSgyWriter = new ArrayFile(writtenSgyDir+writtenSgyFileName+".sgy","rw",ByteOrder.BIG_ENDIAN,ByteOrder.BIG_ENDIAN);
      storeEBCDICHeader(_sgyDir,_sgyFileName);
      storeBinaryHeader(_sgyDir,_sgyFileName);

      trace("Get All Trace Headers");
      getAllTraceHeaders(_sgyDir,_sgyFileName);
      afNewSgyWriter.writeBytes(_ebcdicHeader);
      afNewSgyWriter.writeBytes(_binaryHeader);
      int[][] gridI2sgridI3s = _si.getGridIndices();
      int[] gridI2s = gridI2sgridI3s[0];
      int[] gridI3s = gridI2sgridI3s[1];
      int ntrace = _si.countTraces();
      int i2 = 0;
      int i3 = 0;
      for (int it=0; it<ntrace; ++it) {
        i2 = gridI2s[it];
        i3 = gridI3s[it];
        afNewSgyWriter.writeBytes(_traceHeaders[it]);
        afNewSgyWriter.writeInts(floatToIntBitsIBM(f[i3][i2]));
      }
    }
    catch (IOException e) {
      throw new RuntimeException("Cannot write .sgy file: "+writtenSgyDir+writtenSgyFileName);
    }
    System.out.print("Done");
  }

  
  /** 
   * Converts the set .dat file to an array of floats.
   */
  public float[][][] readDat() {
    return convertDatToFloatArray(_datDir,_datFileName);
  }

  public void plotDataAccordingToXY() 
  {
    double[] x = _si.getXs();
    double[] y = _si.getYs();
    SimplePlot sp = new SimplePlot();
    sp.setHLabel("x coordinate (km)");
    sp.setVLabel("y coordinate (km)");
    PointsView pv = sp.addPoints(x,y);
    pv.setMarkStyle(PointsView.Mark.POINT);
    pv.setLineStyle(PointsView.Line.NONE);
    int[] wh = goodWidthHeightD(x,y);
    sp.setSize(wh[0],wh[1]);
  }

  //The 2nd dimension is inline (i2) and the third dimension is crossline (i3).
  public void plotDataAccordingToInlineCrossline() 
  {
    float[] i2 = _si.getI2sAsFloats();
    float[] i3 = _si.getI3sAsFloats();
    SimplePlot sp = new SimplePlot();
    sp.setHLabel("inline sample index i2");
    sp.setVLabel("crossline sample index i3");
    PointsView pv = sp.addPoints(i2,i3);
    pv.setMarkStyle(PointsView.Mark.POINT);
    pv.setLineStyle(PointsView.Line.NONE);
    int[] wh = goodWidthHeightF(i2,i3);
    sp.setSize(wh[0],wh[1]);
  }

  /**
   * Returns the sampling in the 1st dimension in terms of time.
   */
  public Sampling getSampling1() {
    return new Sampling(_n1,_d1,0.0);//assuming starting time of 0s.
  }

  /**
   * Returns the sampling in the 2nd dimension in terms of inline.
   */
  public Sampling getSampling2() {
    return new Sampling(_n2,_d2,_i2min);
  }

  /**
   * Returns the sampling in the 3nd dimension in terms of crossline.
   */
  public Sampling getSampling3() {
    return new Sampling(_n3,_d3,_i3min);
  }

  public void show3D() {
    float[][][] f = convertDatToFloatArray(_datDir,_datFileName);
    double clip = (double) (max(f))*0.02;
    SimpleFrame frame = new SimpleFrame();
    ImagePanelGroup ipg = frame.addImagePanels(f);
    ipg.setClips(-clip,clip);
    frame.getOrbitView().setScale(2.0);
    frame.setSize(1000,1000);
  }

  public void close() 
  {
    _si.close();
  }

  //Private
  ///////////////////////////////////////////////////////////////////////////////////////////
  private SegyImage _si; 
  private int _i1min, _i2min, _i3min;
  private int _i1max, _i2max, _i3max;
  private int _inlineByteLocation, _crosslineByteLocation;
  private double _d1;
  private double _d2 = 1;//Sampling interval in terms of inline
  private double _d3 = 1;//Sampling interval in terms of crossline
  private int _n1, _n2, _n3;
  private byte[] _ebcdicHeader, _binaryHeader;
  private byte[][] _traceHeaders;
  private String _sgyDir, _sgyFileName;
  private String _datDir, _datFileName;

  private float[][][] convertDatToFloatArray(String datDir, String datFileName)
  {
    float[][][] f = new float[_n3][_n2][_n1];
    try {
      ArrayFile af = new ArrayFile(datDir+datFileName+".dat","rw",ByteOrder.BIG_ENDIAN,ByteOrder.BIG_ENDIAN);
      af.readFloats(f); 
    }
    catch (IOException e) {
      throw new RuntimeException("Cannot read .dat file"+"datFileName");
    }
    return f;
  }

  private void getAllTraceHeaders(String sgyDir,String sgyFileName) {
    int ntrace = _si.countTraces();
    int bytesInTraceHeader = 240;
    byte[] traceHeader = new byte[bytesInTraceHeader];
    int bytesPerSample = _si.getBytesPerSample();
    _traceHeaders = new byte[ntrace][bytesInTraceHeader];
    try {
      ArrayFile af = new ArrayFile(sgyDir+sgyFileName+".sgy","r",ByteOrder.BIG_ENDIAN,ByteOrder.BIG_ENDIAN);
      af.seek(3600);
      for (int it=0; it<ntrace; ++it) {
        af.readBytes(traceHeader);
        _traceHeaders[it] = copy(traceHeader);
        af.skipBytes(_n1*bytesPerSample);
      }
      /*for (int it=0; it<ntrace; ++it) {
        byte[] iLBytes = new byte[4];
        byte[] xLBytes = new byte[4];
        for (int i=0; i<4; ++i) {
          xLBytes[i] = _traceHeaders[it][8+i]; 
          iLBytes[i] = _traceHeaders[it][12+i];
        }
        int iL = ByteBuffer.wrap(iLBytes).getInt();
        int xL = ByteBuffer.wrap(xLBytes).getInt();
        System.out.println("iL = "+iL);
        System.out.println("xL = "+xL);

      }
      */
    }
    catch (IOException e) {
      throw new RuntimeException("Cannot read trace headers: "+sgyDir+sgyFileName+".sgy");
    }
  }

  private void storeEBCDICHeader(String sgyDir, String sgyFileName) {
    int startByte = 0;
    int nByte = 3200;
    _ebcdicHeader = new byte[nByte];
    try {
      ArrayFile af = new ArrayFile(sgyDir+sgyFileName+".sgy","r",ByteOrder.BIG_ENDIAN,ByteOrder.BIG_ENDIAN);
      af.readBytes(_ebcdicHeader,startByte,nByte);
    }
    catch (IOException e) {
      throw new RuntimeException("Cannot read EBCDIC Header:"+sgyDir+sgyFileName+".sgy");
    }
  }


  private void storeBinaryHeader(String sgyDir, String sgyFileName) {
    int startByte = 3200;
    int nByte = 400;
    _binaryHeader = new byte[nByte];
    try {
      ArrayFile af = new ArrayFile(sgyDir+sgyFileName+".sgy","r",ByteOrder.BIG_ENDIAN,ByteOrder.BIG_ENDIAN);
      af.seek(startByte);
      af.readBytes(_binaryHeader,0,nByte);
    }
    catch (IOException e) {
      throw new RuntimeException("Cannot read Binary Header");
    }
  }

  private int[] goodWidthHeightD(double[] x, double[] y) 
  {
    double xmin = min(x);
    double ymin = min(y);
    double xmax = max(x);
    double ymax = max(y);
    int w = 1000;
    int h = 1000;
    if ((xmax-xmin)>(ymax-ymin)) {
      h = (int) (h*(ymax-ymin)/(xmax-xmin));
    }
    else {
      w = (int) (w*(xmax-xmin)/(ymax-ymin));
    }
    return new int[]{w,h};
  }

  private int[] goodWidthHeightF(float[] x, float[] y) 
  {
    float  xmin = min(x);
    float  ymin = min(y);
    float  xmax = max(x);
    float  ymax = max(y);
    int w = 1000;
    int h = 1000;
    if ((xmax-xmin)>(ymax-ymin)) {
      h = (int) (h*(ymax-ymin)/(xmax-xmin));
    }
    else {
      w = (int) (w*(xmax-xmin)/(ymax-ymin));
    }
    return new int[]{w,h};
  }
  public int floatToIntBitsIBM(float from)
	{
		int fconv, fmant, t;
		fconv = Float.floatToIntBits(from);
		if (fconv != 0) {
			fmant = (0x007fffff & fconv) | 0x00800000;
			t = (int) ((0x7f800000 & fconv) >> 23) - 126;
			while ((t & 0x3) != 0) {
				++t;
				fmant >>= 1;
			}
      fconv = (0x80000000 & fconv) | (((t>>2) + 64) << 24) | fmant;
		}
		return fconv;
	}

  public int[] floatToIntBitsIBM(float[] from) {
    int n = from.length;
    int[] result = new int[n];
    for (int i=0; i<n; ++i) {
      result[i] = floatToIntBitsIBM(from[i]);
    }
    return result;
  }



  private static void trace(String s) 
  {
    System.out.println(s);
  }

}
