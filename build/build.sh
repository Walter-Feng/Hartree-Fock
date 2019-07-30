GSL_HOME=/usr/local
export PATH=${GSL_HOME}/bin:$PATH
export LD_LIBRARY_PATH=${GSL_HOME}/lib:$LD_LIBRARY_PATH
export INCLUDE=${GSL_HOME}/include:$INCLUDE
export MANPATH=${GSL_HOME}/share/man:$MANPATH

g++ -c ../src/*.cpp -I${GSL_HOME}/include -g
g++ *.o -static -L${GSL_HOME}/lib -lgsl -lgslcblas -lm -o ../HF
