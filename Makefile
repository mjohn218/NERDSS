#
#  Update 2020-01-29:
#  o Uses required argument of serial, omp, or mpi.
#  o Use VPATH for finding cpp file in different directories -- this simplifies rules
#  o Abort if gsl-config isn't available
#  o Fixed (INTEL) compiler search 0=found | 1=notfound ; make conditional simple (ifeq 0|1)
#  o Also use conditional for GCC
#  o For objects, use basename to get base file name
#  o Clean up directory prefix (shorten variable names and group)
#  o Simplified obj and bin rule logic and readability.
#  o put rules in canonical order
#  o Now has PROF for profiling. (This is by default overrided with empty PROF.)
#  o Now uses INCS. CXXFLAGS is used for C++ specific options.
#  o Make executables with suffixes ( nerdss_serial | nerdss_mpi | nerdss_omp).
#  --  a bit cleaner                                            Kent milfeld@tacc.utexas.edu
#
# TODO: use function to create VPATH
# TODO: Fix MPI after learning purpose
# TODO: Make rules for *.hpp's
#
# Set terminal width to 220 to avoid viewing wrapped lines in output. A width of 200 avoids most wrapping.
BDIR   = bin
ODIR   = obj
SDIR   = src
EDIR   = EXEs

PROF   =

.PHONY: any

#               REQUIREMENTS: gls and directories

hasGSL = $(shell type gsl-config >/dev/null 2>&1; echo $$?)
ifeq ($(hasGSL),1)
$(error " GSL must be installed, and gsl-config must be in path.")
else
$(shell mkdir -p bin)
$(shell mkdir -p obj)
endif

#               EXECUTABLE SETUP for serial, MPI
INCLUDE_FOLDERS = boundary_conditions classes error math parser reactions system_setup trajectory_functions io

ifeq (serial,$(MAKECMDGOALS))
	_EXEC = nerdss
endif

ifeq (mpi,$(MAKECMDGOALS))
	_EXEC = nerdss_mpi
         DEFS = -Dmpi_
		 INCLUDE_FOLDERS += debug io_mpi mpi
endif

ifeq (clean,$(MAKECMDGOALS))
	MAKECMDGOALS = dummy
endif

SRCS = $(foreach dir,$(INCLUDE_FOLDERS),$(wildcard $(SDIR)/$(dir)/*.cpp))

         EXEC  = $(patsubst %,$(BDIR)/%,$(_EXEC))


OS    := $(shell uname)
INTEL  = $(shell type icpc  >/dev/null 2>&1; echo $$?)
GCC    = $(shell type g++   >/dev/null 2>&1; echo $$?)

INCS    = $(shell gsl-config --cflags) -Iinclude
CXXFLAGS = -std=c++0x
LIBS     = $(shell gsl-config --libs)
#LIBS     += $(shell pkg-config --libs libprofiler)


#---------------COMPILER SETUP

#               comment out next statement to allow profiling with gprof
override PROF   =

ifeq ($(GCC),0)          # Will use GCC. (See Intel below.)
	CC      = g++
	ifeq (mpi,$(MAKECMDGOALS))
	CC    = mpicxx
	endif
	# CFLAGS  = -O3
	PROF    = -pg -g
#	MPCFLAG = -I /cm/shared/mpi/openmpi/2.1/intel/17.0/include
endif


ifeq ($(INTEL),0)        # Will use Intel, overrides GCC if both present.
	CC      = icpc
	ifeq (mpi,$(MAKECMDGOALS))
	CC    = mpicxx
	endif
	MPCC    = mpicxx
	# CFLAGS  = -O3
	PROF    = -pg -g
endif

ifdef ENABLE_PROFILING
	CFLAGS = -DENABLE_PROFILING -g
	LIBS   += -lprofiler
else
	CFLAGS = -O3 #-g
endif

#---------------OBJECT FILES
# Create the OBJS variable maintaining the directory structure
OBJS = $(patsubst $(SDIR)/%.cpp,$(ODIR)/%.o,$(SRCS))

# Debug prints (uncomment if needed)
# $(info SRCS = $(SRCS))
# $(info OBJS = $(OBJS))

#---------------RULES
syntax:
	@echo "------------------------------------"
	@printf '\033[31m%s\033[0m\n' " USAGE: make serial|mpi"
	@echo "------------------------------------"
	exit 0

# Rules: for $(MAKECMDGOALS) serial, mpi, or omp build
# $(EXEC) bin/nerdss, bin/nerdss_mpi or /binnerdss_omp
$(MAKECMDGOALS):$(EXEC)
	@echo "Finished making (re-)building $(MAKECMDGOALS) version, $(EXEC)."

$(EXEC): $(OBJS)
	@echo "Compiling $(EDIR)/$(@F).cpp"
	$(CC) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) -o $@ $(EDIR)/$(@F).cpp $(OBJS) $(LIBS) $(PLANG)
	@echo "------------"

# Updated rule for object files
$(ODIR)/%.o: $(SDIR)/%.cpp
	@echo "Compiling $< to $@"
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) $(CXXFLAGS) $(INCS) $(PROF) -c $< -o $@ $(PLANG) $(DEFS)
	@echo "------------"

clean:
	rm -rf $(ODIR) bin
# Reference: https://www.gnu.org/software/make/manual/html_node/Quick-Reference.html
#            https://www.gnu.org/software/make/
#            https://www.cmcrossroads.com/article/basics-vpath-and-vpath
#            https://www.gnu.org/software/make/manual/html_node/Implicit-Variables.html
