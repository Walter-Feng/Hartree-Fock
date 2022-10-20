# Hartree-Fock
=================Three years later==================

Hey there. Thank you for having interest in this repository. 

While C is usually the first programming language we will learn - well ok maybe it's python nowadays, but at least back in my time :( - I have to admit, compared with other languages, has more technical difficulties reading / writing the codes due to restricted number of tools available. The pointers are indeed powerful, but ... just too hard.

I also confess that the 2-electron integral code here has some problems that makes it unable to calculate integrals that involve l > 2 orbitals correctly, i.e. d orbitals and f orbitals will fail.

There is a new repository which is written in C++ (https://github.com/Walter-Feng/Hartree-Fock-in-CPP) that addresses these problems. No python - C++ is my last resort X(  I still wish people who are interested in such level of programming can acquire some knowledge of static programming. This repository has also got a taste of software engineering, so slightly more sophisticated and more abstract, but better to read and extend the functionalities. I am also aiming to polish the documents and comments that can help more enthusiasts learn Hartree Fock.

Still, thank you to my repository, who gave me abundant experience how to write functioning codes and relatively large projects.

====================================================

This is a Hartree-Fock program that aims to show how Hartree-Fock really works in every process. It is planned to add comments on almost each line in the source code describing what this is going to perform, what this line is intended for, etc. Currently it only supports Restricted Hartree-Fock (RHF) method. Currently there is a Chinese version of ducumentation  `doc/main.tex` along with a compiled pdf file.

It is still at an eary stage, and there might be unexpected errors.

There are some features in this program:

  1. The programming style is completely in C - the only C++ feature used in the program is `new` and `delete`, which is much more convenient than `malloc` in C. Thus it should be of no problem that you can read the source code and understand how Hartree-Fock can be performed.

  2. The only dependency is GSL (GNU Scientific Library, https://www.gnu.org/software/gsl), utilizing its linear algebra functions (solving the eigensystem). Even the gaussian integrals are performed itself - libint is not used.

  3. It is (in principle) able to use all the bases that have been published in Basis Set Exchange (https://www.basissetexchange.org).

  4. No special trick is used all through the program - the only complicated thing should be the linked list, which should be a must to learn in C.

  5. In principle there (will be) comments at each line, which enables the program to be demonstrated by 'debugging', in other words you can use simple debugging programs like gdb or lldb to directly show the whole trace of the program in the source code level.
  
  6. Currently it can use either the simplest mixing method - NO MIXING method, namely the coefficient matrix obtained in each iteration is directly used in the next iteration, and the Anderson's mixing method.

And something bad:

  1. BAD PROGRAMMING STYLE

  2. The output functions still seem to be disordered and scattered in the whole program, which may cause problem in trying to configure the output format and the total amount of information in it.

Currently the optimized program can be compiled using `cmake`:
```
$ CMAKE_PREFIX_PATH=/the/path/to/your/gsl/library cmake .
$ make
```

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

 -0.541790 -0.375920 -1.109660 -0.972600 -0.378526  0.373677  0.000000
 -0.375920 -0.541790 -1.109660 -0.972600  0.378526  0.373677  0.000000
 -1.109660 -1.109660-20.242856 -5.163859 -0.000000  0.028711  0.000000
 -0.972600 -0.972600 -5.163859 -2.439862 -0.000000  0.118289  0.000000
 -0.378526  0.378526 -0.000000 -0.000000 -0.294245  0.000000  0.000000
  0.373677  0.373677  0.028711  0.118289  0.000000 -0.336572  0.000000
  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000 -0.392653



MOs:

MO_NUM: 0,    MO_ENERGY = -20.251526 , occ = 2
MO_NUM: 1,    MO_ENERGY = -1.257797 , occ = 2
MO_NUM: 2,    MO_ENERGY = -0.594117 , occ = 2
MO_NUM: 3,    MO_ENERGY = -0.459773 , occ = 2
MO_NUM: 4,    MO_ENERGY = -0.392653 , occ = 2
MO_NUM: 5,    MO_ENERGY = 0.582270 , occ = 0
MO_NUM: 6,    MO_ENERGY = 0.693198 , occ = 0

Total Energy: -74.965901

MO_LABEL:
[ H1s , H1s , O1s , O2s , O2py , O2pz , O2px ]

MO_COEFF:


    0.005589    0.155627   -0.449167    0.294977    0.000000   -0.769472    0.815072
    0.005589    0.155627    0.449167    0.294977    0.000000   -0.769472   -0.815072
   -0.994215   -0.233746   -0.000000    0.104073    0.000000   -0.125881   -0.000000
   -0.025856    0.844278    0.000000   -0.538407    0.000000    0.820893    0.000000
   -0.000000    0.000000   -0.612738   -0.000000    0.000000    0.000000   -0.960034
    0.004169   -0.122998    0.000000   -0.755924    0.000000   -0.763647   -0.000000
    0.000000    0.000000    0.000000    0.000000    1.000000    0.000000    0.000000
```
