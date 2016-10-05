# Salman Habib project

1. Plot the matter power spectra for WMAP7, WMAP9, and Planck cosmologies and see what the differences are.

2. Try to figure out the degeneracy directions in the parameter space. As one way to understand what is going on, carry out a sensitivity analysis by fixing all other parameters at the mid-point of the cosmological parameter space and by varying the chosen one over its full range.

3. Use the power spectrum emulator and the galaxy power spectrum to measure galaxy bias as the HOD parameters are changed. Make sure you match the cosmologies in the two emulators before you do this.

**Note**: A galaxy bias is given by the ratio of the galaxy power spectrum to the matter power spectrum (this equals the bias squared). The bias is a function of the wavenumber k, and should tend to a constant value at low k (i.e., become independent of k).

# Instructions to install the emulator

Visit this webpage http://www.hep.anl.gov/cosmology/CosmicEmu/emu.html. Get familiarized with the page and read the overview. 
Then download the FrankEmu (v. 2) from [here](http://www.hep.anl.gov/cosmology/CosmicEmu/CosmicEmu_v2.tar.gz). The code works currently in C. A future release will include a Python
wrapper. The Gnu Scientific Library (GSL) is required for compilation.

## Install [GSL](https://www.gnu.org/software/gsl/)

The GNU Scientific Library (GSL) is a numerical library for C and C++ programmers. It is free software under the GNU General Public License.

The library provides a wide range of mathematical routines such as random number generators, special functions and least-squares fitting. There are over 1000 functions in total with an extensive test suite.

To install this library in ubuntu 14.04 run:

```
sudo apt-get install gsl-bin
sudo apt-get install libgsl0ldbl
sudo apt-get install libgsl0-dev
```

In 16.04 they have erased *libgsl0ldbl* so you will have to install the development version:

```
sudo apt-get install gsl-bin
sudo apt-get install libgsl0-dev
```

The resources/docs can also be installed using: 

```
sudo apt-get install gsl-doc-info gsl-doc-pdf gsl-ref-html gsl-ref-psdoc
```

## Compile the code 

 A makefile is provided that assumes the GSL is 
available and stored in a usual place.  If you have stored it 
somewhere else, then you're probably more qualified than any of us to fix
the compilation so that it works. To compile the code just unpack the tar file you just downloaded, enter the 
folder any type 

```
make
```

Now you just have to run the emu.exe which will prompt you for an output file, parameter
values for each input, and the desired type of output.  The header
information contains the parameters and other information calculated
from the parameters.  The main output will be 2 columns. The first
column is k, the next column is the spectrum (in one of three forms)
for the chosen redshift z.

To run the emu.exe file just type 

```
./emu.exe
```

Make sure to read the readme.txt for further instructions and to see how to use the code.
