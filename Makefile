
CFLAGS = $(shell gsl-config --cflags) -std=c++0x
##CFLAGS += -ggdb -gdwarf-2 -Wall -I include -O
CFLAGS +=-O3 -I include
MPCFLAG = -I /cm/shared/mpi/openmpi/2.1/intel/17.0/include
LIBS   = $(shell gsl-config --libs)

##CC     = icpc
##MPCC 	= mpiicpc

BDIR   = bin
ODIR   = obj
SDIR   = src
EDIR   = EXEs
EPDIR   = EXE_PAR
OS    := $(shell uname)
INTEL := $(shell which icpc)

ifndef INTEL
	CC = g++
	MPCC = mpicxx
else
	CC = icpc
	MPCC = mpiicpc
endif

##ifeq ($(OS),Linux)
##	_OBJS = $(shell find $(SDIR) -name "*.cpp" | sed -r 's/(\.cc|.cpp)/.o/')
##else
	_OBJS = $(shell find $(SDIR) -name "*.cpp" | sed -E 's/(\.cc|.cpp)/.o/' | sed -E 's/src/./')
##endif

_EXECUTABLES = 	nerdss \
##		template \
##		test_angles \
##		test_loops \
##		rd_gen_reweightPBC \
		rd_boundfraconly \

_PAREXECUTABLES = 	nerdss_mpi \



OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))

EXECUTABLES = $(patsubst %,$(BDIR)/%,$(_EXECUTABLES))

PAREXECUTABLES = $(patsubst %,$(BDIR)/%,$(_PAREXECUTABLES))



_SOURCES = ${_EXECUTABLES:=.cpp}
SOURCES = $(patsubst %,$(EDIR)/%,$(_SOURCES))

_PARSOURCES = ${_PAREXECUTABLES:=.cpp}
PARSOURCES = $(patsubst %,$(EPDIR)/%,$(_PARSOURCES))



all: dirs $(EXECUTABLES) 
#$(PAREXECUTABLES)

$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $<"
	$(CC) $(CFLAGS) $(CFLAGS2) -c $< -o $@  
	@echo "------------"

$(EXECUTABLES): $(OBJS)
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CC) $(CFLAGS) $(CFLAGS2) -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS)
	@echo "------------"

$(PAREXECUTABLES): $(OBJS)
	@echo "Compiling $(EPDIR)/$(@F).cpp"
	$(MPCC) $(CFLAGS) $(CFLAGS2) -o $@ $(EPDIR)/$(@F).cpp $(OBJS) $(LIBS)
	@echo "------------"


dirs: 
	mkdir -p bin
	mkdir -p obj
	mkdir -p obj/classes
	mkdir -p obj/math
	mkdir -p obj/parser
	mkdir -p obj/shared
	mkdir -p obj/reactions
	mkdir -p obj/system_setup
	mkdir -p obj/boundary_conditions
	mkdir -p obj/trajectory_functions
	mkdir -p obj/io


clean:
	rm -rf $(ODIR) bin


