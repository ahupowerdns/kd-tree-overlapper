kd: kd.cpp
	g++ -g  -std=c++14 \
		-Wl,-rpath=../flann/lib/ -lz -lflann  \
						-I ../flann/src/cpp \
						-I ../seqan/include \
						-L../flann/lib/ \
						 -fopenmp -fpermissive -o kd kd.cpp
		