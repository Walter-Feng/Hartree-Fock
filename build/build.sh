g++ -c ../src/*.cpp -I /usr/local/include
g++ *.o -g -static -lgsl -lgslcblas -lm -o ../HF
