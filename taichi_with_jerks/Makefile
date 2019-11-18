OPT += -DFMM
OPT += -DOUTPUTPOT
OPT += -DOUTPUTACC
OPT += -DSIMD_P2P
OPT += -DDOUBLE_P2P
OPT += -DSIMD_M2L
OPT += -DEXPANSION=24
OPT += -DMINIBALL
OPT += -DDEBUG

CXX = g++ -funroll-loops -Wfatal-errors -O3 -Wno-format -march=native -fopenmp -fcx-limited-range #-mavx -mavx2 -mfma -mavx512f

#CXX = g++ -funroll-loops -Wfatal-errors -O3 -Wno-format -mavx -march=native -fopenmp -mavx -mavx2 -mfma

#INCL += -I/usr/local/include -DH5_USE_16_API
#INCL += -L/usr/local/lib -lhdf5 -lz

# hdf5 location for Odyssey 
INCL += -I/n/sw/fasrcsw/apps/Core/hdf5/1.8.12-fasrc04/ -DH5_USE_16_API
INCL +=	-L/n/sw/fasrcsw/apps/Core/hdf5/1.8.12-fasrc04/ -lhdf5 -lz

all:
	@make Taichi

Taichi: taichi.cxx
	$(CXX) $? -o $@ $(OPT) $(INCL)	

clean:
	$(RM) ./*.o ./Taichi ./fmm ./*~
