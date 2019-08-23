# Hartree-Fock
This is a Hartree-Fock program that aims to show how Hartree-Fock really works in every process. It is planned to add comments on almost each line in the source code describing what this is going to perform, what this line is intended for, etc. Currently it only supports Restricted Hartree-Fock (RHF) method.

It is still at an eary stage, and there might be unexpected errors.

There are some features in this program:

  1. The programming style is completely in C - the only C++ feature used in the program is `new` and `delete`, which is much more convenient than `malloc` in C. Thus it should be of no problem that you can read the source code and understand how Hartree-Fock can be performed.

  2. The only dependency is GSL (GNU Scientific Library, https://www.gnu.org/software/gsl), utilizing its linear algebra functions (solving the eigensystem). Even the gaussian integrals are performed itself - libint is not used.

  3. It is (in principle) able to use all the bases that have been published in Basis Set Exchange (https://www.basissetexchange.org).

  4. No special trick is used all through the program - the only complicated thing should be the linked list, which should be a must to learn in C.

  5. In principle there (will be) comments at each line, which enables the program to be demonstrated by 'debugging', in other words you can use simple debugging programs like gdb or lldb to directly show the whole trace of the program in the source code level.

And something bad:

  1. Currently it uses the simplest mixing method - NO MIXING method, namely the coefficient matrix obtained in each iteration is directly used in the next iteration.

  2. BAD PROGRAMMING STYLE

  3. The output functions still seem to be disordered and scattered in the whole program, which may cause problem in trying to configure the output format and the total amount of information in it.

Currently the optimized program can be compiled using the `build/build.sh` script, but the path for GSL should be stated first in the script (the default install path for GSL should be /usr/local). It is also able to run `build/debug.sh` to get the debug version of the program, which can be used in debugging programs like gdb or lldb. Considering that the program itself is quite simple, you can also try to compile the program all by yourself - writing a makefile or cmakelist should be fine and easy to get it working.

Currently supported flags are:

  `-h` : print the help message
  
  `-f <filename>` : specify the input file
  
  `-F` : print fock matrix (It can also be specified in the input file.)
  
  `-C` : print coefficient matrix (It can also be specified in the input file.)
  
  `--input-help`: print a template input file
  
  P.S.: the input file should have `.in` postfix. (To be honest, any postfix would work, in current version, but it is still suggested not to do so.)
  
Known issues:
  1. the program will receive wrong results if it uses a large basis (e.g. cc-pVDZ), but everything seems to be fine with small bases (e.g. sto-3G and 6-31G)

Future plan should be achieving adding comments on each line. Also it is planned to write a document written in Chinese (English ver. may come later?) describing the mechanisms of each step, which enables understanding of Hartree Fock from a physical/mathematical background.

All comments on how to improve this program is much appreciated.
