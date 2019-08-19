GSL_HOME=/Users/dj19970/usr/local

g++ -c ../src/*.cpp -I${GSL_HOME}/include 
g++ *.o -O3 -L${GSL_HOME}/lib -lgsl -lgslcblas -lm -o ../HF
