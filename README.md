# Hartree-Fock
This is a Hartree-Fock program that aims to show how Hartree-Fock really works in every process. It is planned to add comments on almost each line in the source code describing what this is going to perform, what this line is intended for, etc.

It is at a VERY EARLY stage, TOO EARLY that it even CANNOT get the correct molecular orbital energy correctly.

There are some features in this program:

1. The programming style is completely in C - the only C++ feature used in the program is `new` and `delete`, which is much more convenient than `malloc` in C. Thus it should be of no problem that you can read the source code and understand how Hartree-Fock can be performed.

2. The only dependency is GSL (GNU Scientific Library, https://www.gnu.org/software/gsl), utilizing its linear algebra functions (solving the eigensystem). Even the gaussian integrals are performed itself - libint is not used.

3. It is (in principle) able to use all the bases that have been published in Basis Set Exchange (https://www.basissetexchange.org).

4. No special trick is used all through the program - the only complicated thing should be the linked list, which should be a must to learn in C.

5. In principle there (will be) comments at each line, which enables the program to be demonstrated by 'debugging', in other words you can use simple debugging programs like gdb or lldb to directly show the whole trace of the program in the source code level.

And something bad:

1. As it has been mentioned, the energy result itself is still WRONG currently.

2. Currently it uses the simplest mixing method - NO MIXING method, namely the coefficient matrix obtained in each iteration is directly used in the next iteration.

3. BAD PROGRAMMING STYLE

4. The output functions still seem to be disordered and scattered in the whole program, which may cause problem in trying to configure the output format and the total amount of information int it.

Currently the program can be compiled using the `build/build.sh` script, but the path for GSL should be stated first in the script (the default install path for GSL should be /usr/local). It is also able to run `build/debug.sh` to get the debug version of the program, which can be used in debugging programs like gdb or lldb. Considering that the program itself is quite simple, you can also try to compile the program all by yourself - writing a makefile or cmakelist should be fine and easy to get it working.

All comments on how to improve this program is much appreciated.
