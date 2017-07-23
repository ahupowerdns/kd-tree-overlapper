CXXFLAGS=-Wall -O -std=c++14 -fpermissive -fopenmp		-I $(CURDIR)/../flann/src/cpp \
						-I $(CURDIR)/../seqan/include 
STATICFLAGS= -static -static-libstdc++ -static-libgcc 

all: kd

clean:
	rm -f kd *.o


kd: kd.o
	g++ -g  -std=c++14 \
		-Wl,-rpath=$(CURDIR)/../flann/lib/ -L $(CURDIR)/../flann/lib/ \
						 -fopenmp  -o kd kd.o -lz -lflann_s ${STATICFLAGS}		