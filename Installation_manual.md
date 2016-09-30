# MEXICAN NUMERICAL SIMULATIONS SCHOOL (Sept-Oct 2016), Instituto de Física, UNAM, Ciudad de México

This instruction manual is to help install the required packages, libraries and programs needed to run N-body
   simulations ON A LINUX OPERATING SYSTEM COMPUTER, the following was all done on a computer with Ubuntu 14.04 LTS

I. Preschool (sept 29 & 30)

- The requirements are mentioned at the School's webpage: http://iac.edu.mx/mexsimschool/pre-school/reading-material/
 - gfortran compiler
 - gcc compiler
 - python yt library:         http://yt-project.org/
 - tipsy:                      http://www-hpcc.astro.washington.edu/tools/tipsy/tipsy.html
 - python plplot library
 - gadget viewer

- Previous requirements:
 - python
 - Maybe some python scientific libraries like scipy #correct me if i'm wrong. 
  To install python just type 'sudo apt-get install python' on the terminal, check version with 'python --version'

   1. To install gfortran gcc compiler:
      - go to terminal
      - type 'sudo apt-get install gfortran gcc'
      - if the packages are already installed, check version with 'gfortran --version' and 'gcc --version', 
       should be equivalent
   2. To install the yt library (for python):
      - type 'sudo pip install yt', the output might say something like "successfully installed cython-0.24.1 yt-3.3.1", i 
      really don't know what cython is, but who cares! yt is installed
   3. To install tipsy:
      - download this tar ball ftp://ftp-hpcc.astro.washington.edu/pub/hpcc/tipsy.tar.gz
      - Since tipsy will be installed on the folder you are located, my suggestion is to move the tar you just downloaded             into a /home subfolder (e.g. i did this by error inside the Dropbox folder and eventually it was hard for me to    locate where i'd installed the program)
      - Once the tipsy.tar.gz is in the folder, extract all the files (just type 'tar xvf tipsy.tar.gz')
      - Move to the tipsy directory you just extracted with 'cd tipsy-2.2.3d'
      - Go to the code directory: 'cd code'
      - To install: 'bash configure', this will create a Makefile in the same directory
      - Now run a make to finish installation: 'make'
   4. To install plplot library (for python) and gadget viewer, whilst installing Nbodykit
      - reffer to the instruction manual provided in the school's webpage:
        http://iac.edu.mx/mexsimschool/pre-school/reading-material/
      - you will find a link to a dropbox folder:
        https://www.dropbox.com/sh/wvh6vthsv13jia6/AAC3ZyOQNrHDgmdZgVwIM1O9a/Codes?dl=0 (check:active on sept 30)
      - Download the file: "Installation_steps_Linux-Ubuntu16"
      - Download the "zip" folder that contains all the programs that will be installe    
      - NOW JUST FOLLOW THE INSTRUCTIONS PROVIDED IN THE "Installation_steps...."
    #this will install plplot and nbodykit, I DONT KNOW IF GADGET VIEWER WAS INSTALLED WITH THIS STEPS, PLEASE CORRECT IF
     WRONG.

   5. Now, we require MPI to run Julio's codes, MPI's technical specs can be found here: https://www.open-mpi.org/
      - To install MPI go to the given webpage and click on Downloads or go directly to https://www.open- mpi.org/software/ompi/v2.0/
      - Move the compressed file into the directory of your preference
      - On the terminal, move to that directory and type 'tar xvf the_compressed_file_name'
      - Move to the recently created directory: 'cd openmpi-2.0.1'
      - Open the "INSTALL" document (which has detailed instructions on the installation) with the command 'vi INSTALL': "general build of OpenMPI is a combination of 'configure' and 'make' commands"
      - type './configure', wait for it to finish
      - now type 'make all install'


II. School (Oct 4,5,6 & 7)
 TBA
