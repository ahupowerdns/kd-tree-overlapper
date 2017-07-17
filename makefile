kd: kd.cpp
	g++ -g  -std=c++14 \
		-Wl,-rpath=$(CURDIR)/../flann/lib/  \
						-I $(CURDIR)/../flann/src/cpp \
						-I $(CURDIR)/../seqan/include \
						-L $(CURDIR)/../flann/lib/ \
						 -fopenmp -fpermissive -o kd kd.cpp -lz -lflann_s -static -static-libstdc++ -static-libgcc 
		