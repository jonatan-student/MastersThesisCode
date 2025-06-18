
CC = mkoctfile
FLAGS = -v -Wall -W -Xpreprocessor -fopenmp -O3 -march=native -ffast-math\
        -I$(LIBOMP_PATH)/include \
        -L$(LIBOMP_PATH)/lib -lomp

LIBOMP_PATH = $(shell brew --prefix libomp)
OBJECTS= fdcds2d.oct fdfrwd2d.oct fdbckwd2d.oct \
		 fdscnd2d.oct \
		 fdsor2d.oct \
		 fdv2m2d.oct fdm2v2d.oct \
		 fdcplb1.oct fdcplb2.oct fdflowlb.oct \
		 fdcorrectn2d.oct fdstreamlb.oct fdboundaryextract.oct\
		 fdcollisionlb.oct fdboundaryapply.oct semi_implicit_update.oct\
		 computeQMap_mex.oct steadyStateRDA_mex.oct\

all: $(OBJECTS)

%.oct:%.cpp
	$(CC) $(FLAGS) $<

clean:
	rm -f *.oct
