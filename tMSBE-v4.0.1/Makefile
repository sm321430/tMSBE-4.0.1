
# OLD
#all:
#	g++ -O3 -std=c++11 -fopenmp main.cpp classEcyclic.cpp statcenter.cpp classDevice.cpp classModule.cpp classVCSEL.cpp classCavity.cpp classBoundary.cpp
#	g++ -O3 -fopenmp main.cpp classEcyclic.cpp classDevice.cpp classModule.cpp classVCSEL.cpp classCavity.cpp classBoundary.cpp
#	icc -O3 -openmp main.cpp classEcyclic.cpp classDevice.cpp classModule.cpp classVCSEL.cpp classCavity.cpp classBoundary.cpp
#	/usr/bin/gcc-4.8.2/bin/g++ -O3 -std=c++11 -fopenmp main.cpp classEcyclic.cpp statcenter.cpp classDevice.cpp classModule.cpp classVCSEL.cpp classCavity.cpp classBoundary.cpp
#	g++ -c -O3 -fopenmp main.cpp classEcyclic.cpp classDevice.cpp classModule.cpp classVCSEL.cpp classCavity.cpp classBoundary.cpp	

#echo KMP_AFFINITY
# export KMP_AFFINITY=verbose,scatter
# 

LD        := g++ -fopenmp
CC        := g++ -std=c++0x  -fopenmp -O3 -D EIGEN_DONT_PARALLELIZE
#CC        := g++ -std=c++0x -fopenmp -O3 -Wno-write-strings -Wno-unused-result

MODULES = tMSBE-v4.0
SRC_DIR = $(addprefix src/,$(MODULES))
BUILD_DIR := $(addprefix build/,$(MODULES))

SRC       := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
OBJ       := $(patsubst src/%.cpp,build/%.o,$(SRC))
INCLUDES  := $(addprefix -I,$(SRC_DIR))

vpath %.cpp $(SRC_DIR)

define make-goal
$1/%.o: %.cpp
	$(CC) $(INCLUDES) -c $$< -o $$@
endef

.PHONY: all checkdirs clean

all: checkdirs removeOld build/a.out run/a.out configIsak

removeOld:
	@rm -f build/a.out
	@rm -f run/a.out

build/a.out: $(OBJ)
	$(LD) $^ -o $@

checkdirs: $(BUILD_DIR) run_dir

$(BUILD_DIR):
	@mkdir -p $@

run_dir:
	@mkdir -p run/
	
configIsak: $(BUILD_DIR)
	@cp material/material_* run/
	@cp material/cavity.config run/
	@cp material/*.py run/
	@cp myjob.sh run/

run/a.out:
	@cp build/a.out run/

clean:
	@rm -rf $(BUILD_DIR)
	@rm -rf run/

restart:
	@rm -r build/*
	@rm run/a.out
	make

boole:
	make -f Makefile.boole
walton:
	make -f Makefile.walton
hamilton:
	make -f Makefile.hamilton

reflection:
	make -f Makefile.reflection

$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))
