# MEXICAN NUMERICAL SIMULATIONS SCHOOL (Sept-Oct 2016), Instituto de Física, UNAM, Ciudad de México

This instruction manual is to help install the required packages, libraries and programs needed to run N-body
   simulations ON A LINUX OPERATING SYSTEM COMPUTER, the following was all done on a computer with Ubuntu 14.04 LTS
   
# Table of Contents

  * Installation for Preschool
 	* Ubuntu
 	* Mac
  * Preschool Projects
  * School

# Installations for Preschool (sept 29 & 30)

 - gfortran and gcc compiler
 - python yt library:         http://yt-project.org/
 - tipsy:                      http://www-hpcc.astro.washington.edu/tools/tipsy/tipsy.html
 - python plplot through NagBody
 - gadget viewer
 - MPI

### Installation instructions - Ubuntu

   1. Installing gcc/gfortran
      - Go to terminal
      - Type `sudo apt-get install gfortran gcc`
      - check installation with `gfortran --version` and `gcc --version`, 

   2. Installing yt library:
      - type `sudo pip install yt`, the output might say something like "successfully installed cython-0.24.1 yt-3.3.1", i 

   3. Installing tipsy:
      - `mkdir /home/tipsy`
      - `cd /home/tipsy`
      - `wget ftp://ftp-hpcc.astro.washington.edu/pub/hpcc/tipsy.tar.gz`
      - Extract all the files 
      		`tar xvf tipsy.tar.gz`
      - Move to the tipsy directory you just extracted with 
      		`cd tipsy-2.2.3d`
      - Go to the code directory:
      	    `cd code`
      - To install: 
      		`bash configure`
      - Now run a make to finish installation: 
      		`make`
      
   4. Install plplot library and gadgetviewer through NagBody
   
      -  Follow the INSTALL file at
      	*https://github.com/FavioVazquez/MexicanNumericalSimulationSchool/tree/master/preschool/Codes
      - Follow the install file "Installation_steps_Linux-Ubuntu16"
      - Download the "zip" folder that contains all the programs that will be installe    
      - Follow the instructions in "Installation_steps...." 
      - Plplot and nbody kit should be installed.

   5. MPI installation - source: https://www.open-mpi.org/
   
      - Download the version for your OS from
      		https://www.open- mpi.org/software/ompi/v2.0/
      - Move the compressed file into the directory of your preference
      - On the terminal, move to that directory and type 
      		`tar xvf the_compressed_file_name`
      - Move to the recently created directory: `cd openmpi-2.0.1`
      - Open the "INSTALL" document (which has detailed instructions on the installation) with the command `vi INSTALL`: "general build of OpenMPI is a combination of 'configure' and 'make' commands"
      - type `./configure`
      - `make all install`
      
### Installation instructions - Mac


### Projects





      

      


# School (Oct 4,5,6 & 7)
 TBA
