GSL_HOME=/Users/dj19970/usr/local

g++ -c ../src/*.cpp -I${GSL_HOME}/include -g
g++ *.o -g -L${GSL_HOME}/lib -lgsl -lgslcblas -lm -o ../HF
