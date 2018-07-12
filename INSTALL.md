Install Instructions                                                 {#install}
====================

## Build dependencies

C-M.E.S.S. uses the [CMake](http://www.cmake.org) build system to configure the
source code. In contrast to the GNU Autotools this is a relatively new tool,
but it is used by many, even large, projects like KDE or MySQL.
Unfortunately you need to install CMake before you are able to build C-M.E.S.S.

C-M.E.S.S. requires a C and Fortran Compiler supporting at least OpenMP 2.5
 * g++ + gcc + gfortran >= 4.8
 * icc+ifort >= Composer XE 2013
 * clang >= 3.9

The key library dependencies are:
 * LAPACK > 3.2.0
 * SuiteSparse > 3.7.0

Optionally C-M.E.S.S. depends on:
 * Python >= 2.7 with NumPy (>= 1.6) and SciPy (>= 0.9.0)
 * libmatio >= 1.5.0
 * Doxygen
 * ARPACK
 * Serveral compression libraries like zlib, bz2, lzma



## Preparation on Systems with Package-Managers

The dependencies can be installed as well using the package management:

On Debian(>=7.0)/Ubuntu(>=14.04):

    sudo apt-get install build-essential cmake gfortran libblas-dev liblapack-dev \
         libsuitesparse-dev zlib1g-dev  libbz2-dev liblzma-dev libarpack2-dev libmatio-dev \
         doxygen python3-dev python-dev python-scipy python3-scipy

If the cmake package installs an old (<3.0.2) version of cmake check if the `cmake3` is available.

On Redhat/Fedora/CentOS or Scientific Linux:

    sudo yum install gcc-gfortran cmake blas-devel lapack-devel suitesparse-devel zlib-devel
         bzip2-devel rpm-build lzma-devel arpack-devel python-devel numpy scipy

If the cmake package installs an old (<3.0.2) version of cmake check if the `cmake3` is available.

On FreeBSD:

    portinstall gcc cmake blas lapack suitesparse

or via the pkg system

    pkg install gcc cmake blas lapack suitesparse

On MacOSX:

    brew install cmake gcc suite-sparse

You have to specify the compilers on MacOSX, e.g.

    CC=gcc-7 CXX=g++-7 FC=gfortan-7 cmake ../ <OPTIONS>

For an optimal performance we recommend to install an optimized and threaded
BLAS library like Intel MKL, OpenBLAS or ATLAS.

## Build and install MESS

If you have installed or build all dependencies, you can build the MESS library.
We recommend an out-of-source build to keep the sources clean.

    mkdir mess-build
    cd mess-build

From inside this build directory you can now run CMake

    cmake ../ <OPTIONS>

and build and install MESS using

    make
    make install

## Options for CMake
Options that influences the behavior of the configuration process are:

| Parameter                     | Description
| ------------------------------| --------------------------------------------
| -DCMAKE_INSTALL_PREFIX=/path/ | Install path for make install. Default /usr/local/
| -DMESS64=ON/OFF               | Turn on 64-bit integers in MESS This requires all external libraries to work with 64bit integers too. Default: OFF
| -DDEBUG=ON/OFF                | Enables the debug build. Optimizations are turned off, verbose output of the library, debug symbols are included. Default: OFF
| -DCOVERAGE=ON/OFF             | Enable the code coverage, only useful with the DEBUG option and the automatic tests. Default: OFF
| -DTESTS=ON/OFF                | Build software tests. Default: ON
| -DLARGETESTS=ON/OFF           | Build long running software tests. Default: OFF
| -DHUGETESTS=ON/OFF            | Build an extensiv test. Try to test all possible combination of LRADI/LRNM options. Default: OFF
| -DPYTHON=ON/OFF               | Try to build the Py-M.E.S.S. interface. It searches for all available Python interpreters and if the support NumPy and SciPY. Default: OFF
| -DPYTHON_VERSION=x.y          | Search only for the given Python version x.y instead of all.
| -DMATIO=ON/OFF                | Enables MATIO support to read form Matlab files. If MATIO is compiled with HDF5 support, MESS can load Matlab v7.3 files. Default: OFF
| -DDOC=ON/OFF                  | Enable the <tt>make doc</tt> target to build the documentation. Default: OFF.
| -DDOCPDF=ON/OFF               | Enable the <tt>make doc</tt> target to build the documentation as PDF. Default: OFF.
| -DMATHJAX=ON/OFF              | Enable MATHJAX for documentation. Default: OFF.
| -DSUPERLU=ON/OFF              | Enable SuperLU support. Default OFF.
| -DSUPERLU_ROOT=/path/         | Give CMake a path, where to find your SuperLU  build directory.
| -DARPACK=ON/OFF               | Enable Arpack support. Default OFF.
| -DLD_MULTIPLICITY=NUMBER      | Set Leading Dimension of Dense Matrices to a multiple of NUMBER. Default 1.
| -DSUITESPARSE=/path/          | Give CMake a path, to your SuiteSparse build directory.
| -DBLAS_LIBRARIES=lib          | Specify a BLAS library.
| -DLAPACK_LIBRARIES=lib        | Specify a LAPACK library.
| -DMATLAB=ON/OFF               | Enable Mex-M.E.S.S. Interface
| -DMATLABROOT=/path/           | Give CMake a path, to your Matlab installation path

If CMake detects the wrong compilers you can specify them by setting the *CC*,
*CXX* and *FC* environment variables. In this way the Intel Compiler suite can
be enabled using

    CC=icc CXX=icpc FC=ifort cmake ../ <OPTIONS>

## Build the Documentation
In order to build the documentation you have to install
[doxygen](http://www.doxygen.org) before you run the configure step. If
`-DDOC=ON` is given as an argument to `CMake` the *make doc* target is enabled.
The documentation is build in `doc/html` by calling

    make doc

Caused by a problem in Doxygen the bibliography is only integrated in the HTML
page if the installed Doxygen version is higher than 1.8.3.

## Run the Software-Tets
If the tests are enabled you can verify your configuration after the build process using

    make test

The coverage can be enabled using the `COVERAGE` option and needs `gcov` and
`lcov` to be installed on the system. The code coverage is done using

    make coverage

If valgrind was found during the configuration, you can run all tests with
valgrind using

    make memtest


##  Py-M.E.S.S. Interface and Sphinx Documentation
In order to build the Py-M.E.S.S. Documentation you need `sphinx`, `numpydoc`, `matplotlib` and
a working `Py-M.E.S.S.` installation. You can build the `Py-M.E.S.S.` by setting `PYTHON=ON`.
`Py-M.E.S.S.` needs `numpy` and `scipy`. After a successfull CMake configuration type

    make

and change into the directory `python` and using the `setup.py` files to install the Python interface

    cd python/python_X.Y
    python setup.py build && sudo python setup.py install

For a local installation create a directory `mkdir ~/pymess_install` compile `C-M.E.S.S.` as described above and execute the following commands

    cd python/python_X.Y
    python setup.py build && python setup.py install --prefix=~/pymess_install

Do not forget to adapt the `PYTHONPATH` variable

    export PYTHONPATH=/home/daniels/pymess_install/lib/python2.7/site-packages:$PYTHONPATH


If `CMake` found all necessary tools for the Py-M.E.S.S. documentation a special target for your `Python` executable
is activated to generate the documentation:

    make pymess-python2.7-doc-html
    make pymess-python2.7-doc-pdf


##  Mex-M.E.S.S. Interface
To generate Mex-M.E.S.S the Matlab interface to C-M.E.S.S. set `MATLAB=ON` and
give your Matlab installation directory via `MATLABROOT=/path/`.
Make sure that you use compilers which are compatible with your Matlab installation

    CC=gcc-4.9 CXX=g++-4.9 FC=gfortran-4.9 cmake ../cmess -DMATLAB=ON -DMATLABROOT=/path/


If the mex-compiler complains abot your gcc version configure your Matlab mex setup properly.

Start Matlab and change to the C-M.E.S.S. `<build>/matlab` directory.
Execute `mess_path.m` script to add Mex-M.E.S.S. to the matlab path.
Afterwards you can test the Mex-M.E.S.S. interface using

    run_tests()
    run_examples()

function calls.


## Debian Packages
C-M.E.S.S. comes with all scripts required to build a debian packages. Therefore, use

    dpkg-buildpackage -us -uc

to create the `libmess`, `libmess-dev`, and `libmess-doc` package for your current system.




