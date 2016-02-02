#############################################################################
#   
#	Makefile for DISSECT
#   
#############################################################################

#############################################################################
# Leave blank after "=" to disable using BOOST libraries; put "= 1" to enable
BOOSTLIB = 1

#############################################################################
# Select MPI implementation
MPI = MPICH2
# MPI = OPENMPI

#############################################################################
# Put C++ compiler here. It is platform specific. 
CXX = mpic++
# CXX = CC
# CXX = icc

#############################################################################
# scaLAPACK implementation
SCALAPACK = MKL
# SCALAPACK = LIBSCI
#SCALAPACK = NETLIB

#############################################################################
# PATHs
MKL_PATH = /opt/intel/mkl
BOOST_PATH = 

#############################################################################
# Dynamic linking
DYNAMIC = 1

#############################################################################
# Additional flags for the compiler
CXXFLAGS = -O3

CRM = rm

##############################################################################
#
#       Build script
#
##############################################################################

#----------------------------------------------------------------------

ifeq ($(MYHOSTNAME),eddie)
  MKL_PATH=/exports/applications/apps/SL6/intel/mkl
  CXX = mpic++
endif

ifeq ($(MYHOSTNAME),archer)
  MKL_PATH=/opt/intel/composerxe/mkl
  CXX = CC
# CXXFLAGS = -I $(MKL_PATH)/include -openmp 
  DYNAMIC = 
endif

#----------------------------------------------------------------------

CXXFLAGS += -fopenmp -m64
LIB = 

# Dynamic linking?
ifneq ($(DYNAMIC), 1)
  CXXFLAGS += -static
endif

ifneq ($(BOOST_PATH),)
  ifeq ($(BOOSTLIB), 1)
    CXXFLAGS += -I $(BOOST_PATH)
  endif
endif

# scaLAPACK implementation?
ifeq ($(SCALAPACK),MKL)
  CXXFLAGS += -I $(MKL_PATH)/include
  
  ifeq ($(MPI), OPENMPI)
    MKLMPILIB = mkl_blacs_openmpi_lp64
  else
    MKLMPILIB = mkl_blacs_intelmpi_lp64
  endif
  
  ifeq ($(DYNAMIC), 1)
    LIB += -L$(MKL_PATH)/lib/intel64/ -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -l$(MKLMPILIB) -Wl,-rpath=$(MKL_PATH)/lib/intel64/
  else
    LIB += $(MKL_PATH)/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group $(MKL_PATH)/lib/intel64/libmkl_intel_lp64.a $(MKL_PATH)/lib/intel64/libmkl_core.a $(MKL_PATH)/lib/intel64/libmkl_intel_thread.a -Wl,--end-group $(MKL_PATH)/lib/intel64/lib$(MKLMPILIB).a
  endif
endif

ifeq ($(SCALAPACK),NETLIB)
    LIB += -L/usr/lib64/openmpi/lib -lscalapack -lmpiblacs -Wl,-rpath=/usr/lib64/openmpi/lib
endif

LIB += -lpthread -lm -ldl

# CXXFLAGS += -I $(MKL_PATH)/include -fopenmp -m64 #-Wall 

OUTPUT = dissect

# Define some specific flags

ifeq ($(BOOSTLIB),1)
  CXXFLAGS += -DBOOSTLIB
endif

SRC = options.cpp communicator.cpp matrix.cpp misc.cpp genotype.cpp grm.cpp reml.cpp gwas.cpp covariancematrix.cpp phenotype.cpp \
	covariate.cpp auxiliar.cpp simulatephenotype.cpp range.cpp message.cpp pca.cpp results.cpp analysis.cpp predictphenotype.cpp \
	test.cpp main.cpp singlereml.cpp multireml.cpp kernel.cpp blockmatrix.cpp
HDR = global.h communicator.h main.h options.h matrix.h misc.h genotype.h grm.h reml.h gwas.h covariancematrix.h phenotype.h \
	covariate.h auxiliar.h simulatephenotype.h range.h message.h pca.h results.h analysis.h predictphenotype.h \
	test.h singlereml.h multireml.h kernel.h blockmatrix.h
OBJ = $(SRC:.cpp=.o)

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean:
	$(CRM) *.o *~ dissect -f

