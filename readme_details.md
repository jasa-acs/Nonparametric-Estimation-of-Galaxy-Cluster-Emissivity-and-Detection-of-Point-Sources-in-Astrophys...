# ALIAS - Astrophysics Lasso Inverse Abel Solver

Astrophysicists are interested in recovering the 3D gas emissivity of a galaxy cluster from a 2D image taken by a telescope. A blurring phenomenon and presence of point sources make this inverse problem even harder to solve. The current state-of-the-art technique is two step: First identify the location of potential point sources, then mask these locations and deproject the data.  

We instead model the data as a Poisson generalized linear model (involving blurring, Abel and wavelets operators) regularized by two lasso penalties to induce sparse wavelet representation and sparse point sources. The amount of sparsity is controlled by two quantile universal thresholds. As a result, our method outperforms the existing one.

## How to use

### Prerequisites
To build the program, you will need the following tools installed on your linux machine:
* make >= 3.0
* g++ >= 8.0
* cfitsio >= 3.47
* CCfits >= 2.5

On Windows machines, you'll need g++ for Windows, e.g. mingw-w64. Support for Visual Studio is not provided, but feel free to adapt to project for your own needs.

### Build
On linux systems, just use the make tool to build the program:
```
make
```
By default, the program is build with:
```
CC := g++
CFLAGS := -std=c++17 -fopenmp -pedantic -Wall -m64 -march=native -mno-sse5
COPTFLAGS := -O3
CLIB := -lstdc++fs -lcfitsio -lCCfits
```
These variables can of course be changes by passing them to the make utility.
For example, to build with debug flags:
```
make COPTFLAGS="-Og -g"
```
Or to use a different compiler:
```
make CC=g++-9
```
This will create the ALIAS binary in the current directory.


On Windows systems, the easiest is to download Code::Blocks and open the project file. Once the compiler has been linked with Code::Blocks, you can build the project with the desired flavour (Debug, DebugOg, Release, etc...)

### Usage
```
ALIAS -f|--source SOURCE -e|--sensitivity SENSITIVITY -o|--background BACKGROUND
      -b|--blurring BLURRING -r|--result RESULT -s|--size SIZE -x|--bootstrap BOOTSTRAP
```
where
* SOURCE is the path to the source image
* SENSITIVITY is the path to the sensitivity image
* BACKGROUND is the path to the background image
* BLURRING is the path to the blurring filter, defaults to data/blurring.data
* RESULT is the path to the solution file
* SIZE is the width of the picture, assumed square and a power of 2
* BOOTSTRAP is the amount of bootstraps to perform, 1 for no bootstrap

File format currently supports raw binary files and fits format. For fits files, the file name must end in '.fits' (capitalized or not), all other extensions are considered raw binary files.
