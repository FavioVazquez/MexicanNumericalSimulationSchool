# MEXICAN NUMERICAL SIMULATIONS SCHOOL (Sept-Oct 2016), Instituto de Física, UNAM, Ciudad de México

This instruction manual is to help install the required packages, libraries and programs needed to run N-body
   simulations ON A LINUX OPERATING SYSTEM COMPUTER, the following was all done on a computer with Ubuntu 14.04 LTS
   
# Table of Contents

  * Installation for Preschool
 	* Ubuntu
 	* Mac

# Installations

 - gfortran and gcc compiler
 - python [yt]( http://yt-project.org/) library:
 - [tipsy](http://www-hpcc.astro.washington.edu/tools/tipsy/tipsy.html)
 - 14 codes through NagBody 
 - MPI
 - [Gadget Viewer](http://astro.dur.ac.uk/~jch/gadgetviewer/index.html)

### Installation instructions - Ubuntu

Note that you must have Python installed in order to install and run yt

   1. Installing gcc/gfortran
      - Go to terminal
      - Type `sudo apt-get install gfortran gcc`
      - check installation with `gfortran --version` and `gcc --version`, 

   2. Installing yt library:
      - type `sudo pip install yt`, the output might say something like "successfully installed cython-0.24.1 yt-3.3.1"

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
      
   4. Install the 14 codes necessary through NagBody -- [available here](https://www.dropbox.com/sh/wvh6vthsv13jia6/AAC3ZyOQNrHDgmdZgVwIM1O9a/Codes?dl=0)
   
      -  Follow this [INSTALL](https://github.com/FavioVazquez/MexicanNumericalSimulationSchool/blob/master/preschool/MarioNagBodyCodes/INSTALL) file
      - Follow the install file [Installation_steps_Linux-Ubuntu16](https://github.com/FavioVazquez/MexicanNumericalSimulationSchool/blob/master/preschool/MarioNagBodyCodes/Installation_steps_Linux-Ubuntu16)
      - Download the [tar packages](https://github.com/FavioVazquez/MexicanNumericalSimulationSchool/tree/master/preschool/MarioNagBodyCodes/zip)  i.e. all the programs that will be installed.
      - Follow the instructions in "Installation_steps...." 
      - Plplot and nbody kit should be installed.

   5. MPI installation - [source](https://www.open-mpi.org/)
   
      - Download the version for your OS from
      		https://www.open-mpi.org/software/ompi/v2.0/
      - Move the compressed file into the directory of your preference
      - On the terminal, move to that directory and type 
      		`tar xvf the_compressed_file_name`
      - Move to the recently created directory: `cd openmpi-2.0.1`
      - Open the "INSTALL" document (which has detailed instructions on the installation) with the command `vi INSTALL`: "general build of OpenMPI is a combination of 'configure' and 'make' commands"
      - type `./configure`
      - type `make all install`, this should take a while.
      
   6. Gadget Viewer installation - [source](http://astro.dur.ac.uk/~jch/gadgetviewer/index.html)
      - Download the tar ball from [website](http://astro.dur.ac.uk/~jch/gadgetviewer/index.html) or from our github [directory](https://github.com/FavioVazquez/MexicanNumericalSimulationSchool/blob/master/preschool/gadgetfileviewer/gadgetviewer-1.0.6.tar.gz).
      - Move to /home and unpack with `tar xvf the_compressed_file_name`
      - The [README](https://github.com/FavioVazquez/MexicanNumericalSimulationSchool/blob/master/preschool/gadgetfileviewer/README) file specifies that we must have C compiler and fortran90 compiler, tools that, by that point we must have installed if followed the guide. It also specifies that we must have GTK+ 2.0 GUI Library:
      - Now, it seems from GTK+'s [webpage](http://www.gtk.org/download/linux.php) that we need MORE packages to even build GTK+, this packages are glib, gobject-introspection, pango, gdk-pixbuf and atk, you can find all these packages @ the above webpage (we recommend downloading the latest version of each)
      - The same packages (versions from c. October 2016) are in [this]() link, download and unpack.
      - There is a dependency tree that prohibits you from installing Gadgetviewer from the start:
         1. Gadgetviewer dependencies: GTK+2.0 GUI library, HDF5, libpng, PLplot
         2. GTK+2.0 dependencies: GLib, GObject-instrospection, Pango, Gdk-Pixbuf, ATK
         3. GLib dependencies: libpcre
         
      - We'll start by install HDF5, we downloaded the latest stable version from [here](https://launchpad.net/ubuntu/xenial/+source/hdf5),version 1.8.1, unpack and open the terminal, type  `./configure --prefix=/usr/local/hdf5` , after the configuration is done, run `make`, then `make check` (this might take a while), `make install` and finally `make check-install`
      -Now, to install libpng google "libpng" or download from [here](https://sourceforge.net/projects/libpng/files/libpng16/1.6.25/libpng-1.6.25.tar.gz/download) version 1.6.25, it's again a combination of `./configure`, `make check` and `make install`
      -PLplot was previously installed in step 4
      -To install GTK+2.0 you need to install the dependencies in the next order: firstly, libpcre that can be downloaded from [here](https://sourceforge.net/projects/pcre/?source=typ_redirect), remember this is a dependency of glib. When installed, proceed to install GLib
      
      
### Installation instructions - Mac

First step is install Xcode, you can download it from App Store.
Then you have to install de command line tools by typing in the terminal 

```
xcode-select --install
```

The command line tools contain a gcc compiler and at this point you should have already installed python by default in OS X.

In order to install gfortran, you can download it and install it  from 
https://gcc.gnu.org/wiki/GFortranBinaries#MacOS

You must get quartz from
https://www.xquartz.org
Download it and install 

**NOTE**: There are several options to install the libraries and programs, the most common and recommended are Homebrew and Mac Ports, the major difference between them is the path that use to install, and even if you use both, you should not have problems with compatibility because they install in different directories. 

The instructions provided by Mario, are thought to install everything using MacPorts, and I recommend to follow in that way, but if you decide to use homebrew, I will give you the instructions below , but the plplot version that Mario's codes requiere is easily installed using MacPorts, (using homebrew I was not able to solve the errors that).

Note: At some point the Mario's instrucctions ask you to change you profile the name of that file can be .profile or .bash_profile

- USING HOMEBREW :

Install home-brew typing in terminal

```
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
```

You will get Open-MPI by typing in the terminal

```
brew install open-mpi
```

In case you prefer use Mpich type

```
brew install mpich
```

**Note**: be careful never install open-mpi and mpich, you have to choose only one.

This version works perfectly with the Cosmic Emulator that Salman Habib provide us. 

- Install GSL

```
brew install gsl
```
- Install gnuplot

```
brew install gnuplot
```

- Install Plplot

```
brew install plplot
```    
