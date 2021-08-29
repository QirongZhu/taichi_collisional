#OPT += -DHOLD
OPT += -DFROST
OPT += -DFMM
OPT += -DOUTPUTPOT
#OPT += -DOUTPUTACC
OPT += -DSIMD_P2P
#OPT += -DDOUBLE_P2P
OPT += -DEXPANSION=16
OPT += -DMINIBALL
OPT += -DDEBUG
OPT += -DUSETBB

CXX = g++-11 -g -funroll-loops -O3 -march=native -fopenmp -fcx-limited-range  -mavx -mfma #-ltcmalloc #-mavx512f

#CXX = g++ -funroll-loops -Wfatal-errors -O3 -Wno-format -mavx -march=native -fopenmp -mavx -mavx2 -mfma

INCL += -I/usr/local/Cellar/tbb/2021.3.0/include
LIBS += -L/usr/local/Cellar/tbb/2021.3.0/lib -ltbb

INCL += -I/usr/local/include -DH5_USE_16_API
LIBS += -L/usr/local/lib -lhdf5 -lz

# hdf5 location for Odyssey 
#INCL += -I/n/sw/fasrcsw/apps/Core/hdf5/1.8.12-fasrc04/ -DH5_USE_16_API
#INCL +=	-L/n/sw/fasrcsw/apps/Core/hdf5/1.8.12-fasrc04/ -lhdf5 -lz

all:
	@make Taichi

Taichi: main.cxx
	$(CXX) $? -o $@ $(OPT) $(INCL) $(LIBS)
	#./Taichi 0 dat.10 500 1
	#./Taichi 1 snapshot_2.hdf5 500 1

clean:
	$(RM) ./*.o ./Taichi ./fmm ./*~
