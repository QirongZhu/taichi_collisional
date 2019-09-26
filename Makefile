OPT += -DHOLD
#OPT += -DFMM
#OPT += -DOUTPUTPOT
OPT += -DOUTPUTACC
OPT += -DSIMD_P2P
OPT += -DEXPANSION=20

CXX = /usr/local/bin/g++-9 -funroll-loops -Wfatal-errors -O3 -Wno-format -mavx -march=native -fopenmp -mavx2 -mfma

INCL += -I/usr/local/include -DH5_USE_16_API
INCL += -L/usr/local/lib -lhdf5 -lz

all:
	@make Taichi

Taichi: main.cxx
	$(CXX) $? -o $@ $(OPT) $(INCL)
	./Taichi 0 dat.10 500 1
	#./Taichi 1 snapshot_2.hdf5 500 1

clean:
	$(RM) ./*.o ./Taichi ./fmm ./*~
