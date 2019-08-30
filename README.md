# Hartree-Fock
This is a Hartree-Fock program that aims to show how Hartree-Fock really works in every process. It is planned to add comments on almost each line in the source code describing what this is going to perform, what this line is intended for, etc. Currently it only supports Restricted Hartree-Fock (RHF) method. Currently there is a Chinese version of ducumentation  `doc/main.tex` along with a compiled pdf file.

It is still at an eary stage, and there might be unexpected errors.

There are some features in this program:

  1. The programming style is completely in C - the only C++ feature used in the program is `new` and `delete`, which is much more convenient than `malloc` in C. Thus it should be of no problem that you can read the source code and understand how Hartree-Fock can be performed.

  2. The only dependency is GSL (GNU Scientific Library, https://www.gnu.org/software/gsl), utilizing its linear algebra functions (solving the eigensystem). Even the gaussian integrals are performed itself - libint is not used.

  3. It is (in principle) able to use all the bases that have been published in Basis Set Exchange (https://www.basissetexchange.org).

  4. No special trick is used all through the program - the only complicated thing should be the linked list, which should be a must to learn in C.

  5. In principle there (will be) comments at each line, which enables the program to be demonstrated by 'debugging', in other words you can use simple debugging programs like gdb or lldb to directly show the whole trace of the program in the source code level.

And something bad:

  1. Currently it can use either the simplest mixing method - NO MIXING method, namely the coefficient matrix obtained in each iteration is directly used in the next iteration, and the Anderson's mixing method 

  2. BAD PROGRAMMING STYLE

  3. The output functions still seem to be disordered and scattered in the whole program, which may cause problem in trying to configure the output format and the total amount of information in it.

Currently the optimized program can be compiled using the `build/build.sh` script, but the path for GSL should be stated first in the script (the default install path for GSL should be /usr/local). It is also able to run `build/debug.sh` to get the debug version of the program, which can be used in debugging programs like gdb or lldb. Considering that the program itself is quite simple, you can also try to compile the program all by yourself - writing a makefile or cmakelist should be fine and easy to get it working.

The program can also be fetched from release, containing a static program and corresponding basis sets. Only Mac OS X and Linux are supported (for I do not know how to compile this stuff with GSL dependency in Windows - and moreover, I guess it is really hard to use such an command-line-interface program in Windows).

Currently supported flags are:

  `-h` : print the help message
  
  `-f <filename>` : specify the input file
  
  `-F` : print fock matrix (It can also be specified in the input file.)
  
  `-C` : print coefficient matrix (It can also be specified in the input file.)
  
  `--input-help`: print a template input file
  
  P.S.: the input file should have `.in` postfix. (To be honest, any postfix would work, in current version, but it is still suggested not to do so.)
  
Known issues:
  1. the program will receive wrong results if it uses a large basis (e.g. cc-pVDZ), but everything seems to be fine with small bases (e.g. sto-3G and 6-31G)

All comments on how to improve this program is much appreciated.

Bibliography:

1. Obara, Shigeru, and A. Saika. "Efficient recursive computation of molecular integrals over Cartesian Gaussian functions." The Journal of chemical physics 84.7 (1986): 3963-3974.

2. Helgaker, Trygve, Poul Jorgensen, and Jeppe Olsen. Molecular electronic-structure theory. John Wiley & Sons, 2014.

3. May, Andrew James. Density fitting in explicitly correlated electronic structure theory. Diss. University of Bristol, 2006.

4. Johnson, Duane D. "Modified Broyden’s method for accelerating convergence in self-consistent calculations." Physical Review B 38.18 (1988): 12807.

5. Schuchardt, Karen L., et al. "Basis set exchange: a community database for computational sciences." Journal of chemical information and modeling 47.3 (2007): 1045-1052.

A typical result:
```
$ ./HF -f test/H2O.in (in STO-3G)
```
```
...(iteration information)..。

SCF converged.

Fock matrix:

 -0.437615 -0.404153 -1.133001 -1.022207 -0.415453  0.393827  0.000000
 -0.404153 -0.437615 -1.133001 -1.022207  0.415453  0.393827  0.000000
 -1.133001 -1.133001-20.146148 -5.145982 -0.000000  0.027780  0.000000
 -1.022207 -1.022207 -5.145982 -2.331408 -0.000000  0.111017  0.000000
 -0.415453  0.415453 -0.000000 -0.000000 -0.227014  0.000000  0.000000
  0.393827  0.393827  0.027780  0.111017  0.000000 -0.275262  0.000000
  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000 -0.296996



MOs:

MO_NUM: 0,    MO_ENERGY = -20.154493 , occ = 2
MO_NUM: 1,    MO_ENERGY = -1.208528 , occ = 2
MO_NUM: 2,    MO_ENERGY = -0.557811 , occ = 2
MO_NUM: 3,    MO_ENERGY = -0.379794 , occ = 2
MO_NUM: 4,    MO_ENERGY = -0.296996 , occ = 2
MO_NUM: 5,    MO_ENERGY = 0.963347 , occ = 0
MO_NUM: 6,    MO_ENERGY = 1.048381 , occ = 0

Total Energy: -75.003904

MO_LABEL:
[ H1s , H1s , O1s , O2s , O2py , O2pz , O2px ]

MO_COEFF:


    0.003630   -0.208665    0.422342   -0.230836    0.000000   -0.778788   -0.829289
    0.003630   -0.208665   -0.422342   -0.230836    0.000000   -0.778788    0.829289
   -0.994402    0.222023    0.000000   -0.106497    0.000000   -0.142506   -0.000000
   -0.024274   -0.772323   -0.000000    0.522462    0.000000    0.898434    0.000000
   -0.000000   -0.000000    0.643727    0.000000    0.000000   -0.000000    0.939536
    0.003189    0.128637   -0.000000    0.800802    0.000000   -0.715459   -0.000000
    0.000000    0.000000    0.000000    0.000000    1.000000    0.000000    0.000000
```
