#
#  Update 2020-01-06: 
#  o Use VPATH for finding cpp file in different directories -- this simplifies rules
#  o Abort if gsl-config isn't available
#  o Fixed (INTEL) compiler search 0=found | 1=notfound ; make conditional simple (ifeq 0|1)
#  o Also use conditional for GCC
#  o For objects, use basename to get base file name
#  o Clean up directory prefix (shorten variable names and group)
#  o Simplified obj and bin rule logic and readability.
#  o put rules in canonical order
#  -- now the makefile looks a bit cleaner                    Kent milfeld@tacc.utexas.edu
#
# TODO: use function to create VPATH
# TODO: Check this out on OSX
# TODO: Fix MPI after learning purpose
#
VPATH=src/boundary_conditions:src/classes:src/io:src/math:src/parser:src/reactions:src/system_setup:src/trajectory_functions

BDIR   = bin
ODIR   = obj
SDIR   = src
EDIR   = EXEs
EPDIR  = EXE_PAR

.PHONY: any

#               REQUIREMENTS: gls and directories

hasGSL = $(shell type gsl-config >/dev/null 2>&1; echo $$?)
ifeq ($(hasGSL),1)
	$(error " GSL must be installed, and gsl-config must be in path.")
else
	$(shell mkdir -p bin)
	$(shell mkdir -p obj)
endif

#               EXECUTABLE SETUP for serial, MPI, OpenMP (omp)
#
ifeq (serial,$(MAKECMDGOALS))
	_EXEC = nerdss
endif

ifeq (mpi,$(MAKECMDGOALS))
	_EXEC = nerdss_mpi
         DEFS = -DMPI
endif

ifeq (omp,$(MAKECMDGOALS))
	_EXEC  = nerdss_omp
         DEFS  = -DOMP
         PLANG = -fopenmp
endif

ifeq (clean,$(MAKECMDGOALS))
	MAKECMDGOALS = dummy
endif

         EXEC  = $(patsubst %,$(BDIR)/%,$(_EXEC))


OS    := $(shell uname)
INTEL  = $(shell type icpc       >/dev/null 2>&1; echo $$?)

CFLAGS = $(shell gsl-config --cflags) -std=c++0x
LIBS   = $(shell gsl-config --libs)


#---------------COMPILER SETUP

ifeq ($(GCC),0)          # Found/Use Intel icpc compiler
	CC      = g++
	MPCC    = mpicxx
	CFLAGS += -O3 -I include
#	MPCFLAG = -I /cm/shared/mpi/openmpi/2.1/intel/17.0/include
endif


ifeq ($(INTEL),0)          # Found/Use Intel icpc compiler
	CC      = icpc
	MPCC    = mpicxx
	CFLAGS += -O3 -I include
endif

#---------------OBJECT FILES

ifeq ($(OS),Linux)
       _OBJS = $(shell find $(SDIR) -name "*.cpp" | xargs -n 1 basename | sed -r 's/(\.cc|.cpp)/.o/')
else
       _OBJS = $(shell find $(SDIR) -name "*.cpp" | xargs -n 1 basename | sed -E 's/(\.cc|.cpp)/.o/')
endif

        OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


#---------------RULES

syntax:
	@echo "------------------------------------"
	@printf '\033[31m%s\033[0m\n' "   USAGE: make serial|mpi|omp"
	@echo "------------------------------------"
	exit 0

#             Rules: for $(MAKECMDGOALS)  serial,     mpi, or            omp            build 
#                        $(EXEC)          bin/nerdss, bin/nerdss_mpi or /binnerdss_omp
$(MAKECMDGOALS):$(EXEC)
	@echo "Finished making (re-)building $(MAKECMDGOALS) version, $(EXEC)."

$(EXEC): $(OBJS)
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CC) $(CFLAGS) $(CFLAGS2) -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS) $(PLANG)
	@echo "------------"

obj/%.o: %.cpp
	@echo "Compiling $<"
	$(CC) $(CFLAGS) $(CFLAGS2) -c $< -o $@  $(PLANG) $(DEFS)
	@echo "------------"

clean:
	rm -rf $(ODIR) bin

# Reference: https://www.gnu.org/software/make/manual/html_node/Quick-Reference.html
