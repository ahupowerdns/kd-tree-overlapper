kd: kd.cpp
	g++ -g  -std=c++11 \
		-Wl,-rpath=$(HOME)/src/flann/build/lib/:$(HOME)/opt/hdf5/lib/:$(HOME)/opt/lib/ -lz -lflann -lhdf5 \
						-I $(HOME)/src/flann/src/cpp/ \
						-I $(HOME)/opt/hdf5/include/ \
						-I $(HOME)/src/seqan-library-1.4.2/include/ \
						-L$(HOME)/opt/hdf5/lib/ \
						-L$(HOME)/src/flann/build/lib/ \
						-L$(HOME)/opt/lib/ \
						-I$(HOME)/opt/include/ \
						 -fopenmp -fpermissive -o kd kd.cpp
		