Warping with wavelets code.

#########Software needed to run our code################################################
To use this code, you will need to install jython (to use our .py scripts) and the
java jdk (to run the java code). This code will also require the use of the mines jtk.

jython can be downloaded from http://www.jython.org/.
The mines jtk can be downloaded from https://github.com/dhale/jtk. Please read the 
corresponding readme.txt file.
Java SE JDK 7, which is the version used when writing the .java files, can be downloaded 
from http://www.oracle.com/technetwork/java/javase/downloads
########################################################################################

#########Organization of the folders####################################################
In the cg folder is my workbench (bench folder) for research. Inside the bench folder,
I have the corresponding src folder (obviously for our source code), the build folder 
contains the .class files corresponding to the .java files, and
the other folders and files are for building the .java files with gradle.

Inside the src folder, the only folder you will need is the wwarp folder, which
contains our main research. The other folders are other topics I have built to
increase my understanding of a specific topic, but these folders do not contain our research.
########################################################################################

#########Building our code##############################################################
To build the .java files for the warping with wavelets code, refer to “Building the Mines
JTK” section in the readme.txt file in the mines jtk.
For a Mac OS, the required folders and files to build the .java files are located in our bench folder.
If you are running a Mac OS, you will need to cd into the bench directory and type “sh gradlew” to build the .java files. If you
are using another operating system, please refer to “Building the Mines
JTK” section in the readme.txt file in the mines jtk.
########################################################################################

#########Research#######################################################################
Inside the wwarp folder, you will find files pertaining to our research. 

The java files that correspond to the code used in my 2015 presentation are
WaveletWarpingCBCyclic.java and WaveletWarpingCBGN.java, which correspond
to the cyclic search and Gauss-Newton method, respectively. 
Any code that is commented out in these two files corresponds to experimental code
that has not been fully tested.

Many python files exist that use these two java files, but the python file that was used to create
my images in the presentation are in cwppres15synthetics.py.

Note that my repository does not contain data, if you have data that you want to run 
through the existing scripts, modify the methods that deal with extracting the data
or feel free to create your own method.

If you have any questions, please do not hesitate to contact me at cgrazian@mines.edu.
########################################################################################













