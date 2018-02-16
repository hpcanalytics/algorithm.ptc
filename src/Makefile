#CC = xlc_r
CPP = g++
CFLAGS = -O2 -Wno-deprecated  #-qsmp #-O2 #-qarch=auto -qtune=auto  
BOOST = -I /vbi/packages/boost/include/boost-1_38
	
ptriangle: TrianglePar.cpp TrianglePar.hpp SGraph.hpp Set.hpp utility.hpp
	module load mvapich2/gcc;	mpicxx -o ptriangle TrianglePar.cpp
	

clean: 
	rm -f ptriangle
	