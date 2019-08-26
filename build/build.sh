GSL_HOME=/usr/local

g++ -c ../src/*.cpp -I${GSL_HOME}/include 
g++ *.o -O3 -static -L${GSL_HOME}/lib -lgsl -lgslcblas -lm -o ../HF
