#
# Which compiler to use
#
#COMPILER = INTEL
COMPILER = GNU
#
#
ifeq "$(COMPILER)" "INTEL"
   MPIFORTRAN = mpiifort
   FORTRAN = ifort
   FLAGS = -O3
endif
ifeq "$(COMPILER)" "GNU"
   MPIFORTRAN = mpif90
   FORTRAN = gfortran
   FLAGS = -O3 -ffast-math
endif

# Directory for executables
BIN=../bin

ifeq ($(MAKECMDGOALS),devel)
  FLAGS = -Wall -fcheck=bounds
endif
 
###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################

BINARIES = $(BIN)/topolink $(BIN)/evalmodels $(BIN)/linkcorrelation $(BIN)/filtermodels

MODULES = functionpars.f90 \
          linkedcells.f90 \
          topolink_data.f90 \
          topolink_operations.f90 \
          ioformat.f90 \
          flashmod.f90

all : bindir $(BINARIES)
devel : bindir $(BINARIES)

$(BIN)/topolink: Makefile \
                 modules \
                 cgnewton.f90 \
                 computedpath.f90 \
                 computef.f90 \
                 computeg.f90 \
                 goverlapcell.f90 \
                 goverlap.f90 \
                 gstretch.f90 \
                 initcells.f90 \
                 initguess.f90 \
                 linkstatus.f90 \
                 linkconsistency.f90 \
                 keywords.f90 \
                 lnf.f90 \
                 overlapcell.f90 \
                 overlap.f90 \
                 printdata.f90 \
                 random.f90 \
                 stretch.f90 \
                 test_grad.f90 \
                 title.f90 \
                 topolink.f90
	$(FORTRAN) -o $(BIN)/topolink \
                 $(MODULES) \
                 cgnewton.f90 \
                 computedpath.f90 \
                 computef.f90 \
                 computeg.f90 \
                 goverlapcell.f90 \
                 goverlap.f90 \
                 gstretch.f90 \
                 initcells.f90 \
                 initguess.f90 \
                 linkstatus.f90 \
                 linkconsistency.f90 \
                 keywords.f90 \
                 lnf.f90 \
                 overlapcell.f90 \
                 overlap.f90 \
                 printdata.f90 \
                 random.f90 \
                 stretch.f90 \
                 test_grad.f90 \
                 title.f90 \
                 topolink.f90 \
                 $(FLAGS)

$(BIN)/linkcorrelation : Makefile modules linkcorrelation.f90
	$(FORTRAN) -o $(BIN)/linkcorrelation $(MODULES) linkcorrelation.f90 $(FLAGS)

$(BIN)/filtermodels : Makefile modules filtermodels.f90
	$(FORTRAN) -o $(BIN)/filtermodels $(MODULES) filtermodels.f90 $(FLAGS)

$(BIN)/evalmodels : Makefile modules flashsort.f90 evalmodels.f90
	$(FORTRAN) -o $(BIN)/evalmodels $(MODULES) flashsort.f90 evalmodels.f90 $(FLAGS)

# modules

modules : $(MODULES) 
	$(FORTRAN) -c $(MODULES) $(FLAGS)

bindir : 
	@mkdir -p $(BIN)

clean : 
	rm -f ./*.o ./*.mod $(BIN)/*




