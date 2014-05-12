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

Note, the png and pres14 folders are the destinations for images created for our
CWP report (png) and our presentation (pres14).

In each file, there should be a small description of the purpose of why the file exists
in the first few lines of the file.

The files that begin with cwppres… are the files used to generate the figures for
our wavelets and warping PS seismic images presentation, which was presented at the
2014 CWP annual meeting.

The files that begin with nmo deal with different tests of the warping with wavelets
algorithm with one or multiple CMP gathers.

The  main files that are related to warping PS seismic images are waveletwarping.py
and waveletwarpingha.py. 

The other files represent research in progress or files that I’m
not sure whether I will need later on.

Note that my repository does not contain data, if you have data that you want to run 
through the existing scripts, modify the methods that deal with extracting the data
or feel free to create your own method.

If you have any questions, please do not hesitate to contact me at cgrazian@mines.edu.
########################################################################################













