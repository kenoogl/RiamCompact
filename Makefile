SHELL=/bin/sh

#------------------------------------------------------------------
# library install path
PMDIR=/usr/local/PMlib
TPDIR=/usr/local/TextParser
CPMDIR=/usr/local/CPMlib
CDMDIR=/usr/local/CDMlib

# for Intel Compiler
CXX        = mpic++
CXXFLAGS   = -O3 -openmp -shared-intel -D_CDM_OUTPUT
F90        = mpif90
F90FLAGS   = -O3 -openmp

# for K/FX10
#CXX        = mpiFCCpx
#CXXFLAGS   = -Kfast,openmp -D_CDM_OUTPUT
#F90        = mpifrtpx
#F90FLAGS   = -Kfast,openmp
#------------------------------------------------------------------

.SUFFIXES: .cpp .f90

CXXSRCS = \
          3dtopo-c-PMlib.cpp \
          cpp-module/riamc.cpp \
          cpp-module/riamc_Initialize.cpp \
          cpp-module/riamc_allocateArray.cpp

F90SRCS = \
          fortran-module/moduled-init.f90 \
          fortran-module/moduled-eddy-omp.f90 \
          fortran-module/moduled-interv-omp.f90 \
          fortran-module/moduled-contrav-omp.f90 \
          fortran-module/moduled-SOR-FFVC-omp.f90 \
          fortran-module/moduled-newp-omp.f90 \
          fortran-module/moduled-newv-omp.f90 \
          fortran-module/moduled-display.f90 \
          fortran-module/moduled-display2.f90 \
          fortran-module/copy_array.f90

CXXOBJS = $(CXXSRCS:.cpp=.o)
F90OBJS = $(F90SRCS:.f90=.o)
OBJS    = $(CXXOBJS) $(F90OBJS)

CXXFLAGS  += `$(PMDIR)/bin/pm-config --cflags` \
             `$(TPDIR)/bin/tp-config --cflags` \
             `$(CDMDIR)/bin/cdm-config --cflags` \
             `$(CPMDIR)/bin/cpm-config --cflags`
F90FLAGS  += `$(CPMDIR)/bin/cpm-config --fcflags`
LIBS       = `$(PMDIR)/bin/pm-config --libs` \
             `$(TPDIR)/bin/tp-config --libs` \
             `$(CPMDIR)/bin/cpm-config --libs` \
             `$(CDMDIR)/bin/cdm-config --libs`
LDFLAGS    =  

TARGET = riam-pmlib


$(TARGET):		$(OBJS)
	$(CXX) $(CXXFLAGS) -o $(@) $(OBJS) $(LDFLAGS) $(LIBS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $(@) -c $<

.f90.o:
	$(F90) $(F90FLAGS) -o $(@) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

