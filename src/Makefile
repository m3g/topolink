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

BINARIES = $(BIN)/topolink \
           $(BIN)/evalmodels \
           $(BIN)/linkcorrelation \
           $(BIN)/filtermodels \
           $(BIN)/maxmingyr \
           $(BIN)/modelmaps \
           $(BIN)/pbias

MODULES = functionpars.o \
          linkedcells.o \
          topolink_data.o \
          topolink_operations.o \
          ioformat.o \
          flashmod.o

all : bindir $(BINARIES)
devel : bindir $(BINARIES)

$(BIN)/topolink: Makefile $(MODULES) \
                 topolink.o \
                 cgnewton.o \
                 computedpath.o \
                 computef.o \
                 computeg.o \
                 goverlapcell.o \
                 goverlap.o \
                 gstretch.o \
                 initcells.o \
                 initguess.o \
                 linkstatus.o \
                 linkconsistency.o \
                 keywords.o \
                 lnf.o \
                 overlapcell.o \
                 overlap.o \
                 printdata.o \
                 random.o \
                 stretch.o \
                 test_grad.o \
                 title.o
	$(FORTRAN) -o $(BIN)/topolink $(MODULES) \
                 topolink.o \
                 cgnewton.o \
                 computedpath.o \
                 computef.o \
                 computeg.o \
                 goverlapcell.o \
                 goverlap.o \
                 gstretch.o \
                 initcells.o \
                 initguess.o \
                 linkstatus.o \
                 linkconsistency.o \
                 keywords.o \
                 lnf.o \
                 overlapcell.o \
                 overlap.o \
                 printdata.o \
                 random.o \
                 stretch.o \
                 test_grad.o \
                 title.o \
                 $(FLAGS)

$(BIN)/linkcorrelation : Makefile $(MODULES) linkcorrelation.o
	$(FORTRAN) -o $(BIN)/linkcorrelation $(MODULES) linkcorrelation.o $(FLAGS)

$(BIN)/maxmingyr : Makefile $(MODULES) maxmingyr.o
	$(FORTRAN) -o $(BIN)/maxmingyr maxmingyr.o $(FLAGS)

$(BIN)/modelmaps: Makefile $(MODULES) modelmaps.o
	$(FORTRAN) -o $(BIN)/modelmaps modelmaps.o $(FLAGS)

$(BIN)/pbias : Makefile $(MODULES) pbias.o
	$(FORTRAN) -o $(BIN)/pbias pbias.o $(FLAGS)

$(BIN)/filtermodels : Makefile $(MODULES) filtermodels.o
	$(FORTRAN) -o $(BIN)/filtermodels $(MODULES) filtermodels.o $(FLAGS)

$(BIN)/evalmodels : Makefile $(MODULES) flashsort.o evalmodels.o
	$(FORTRAN) -o $(BIN)/evalmodels $(MODULES) flashsort.o evalmodels.o $(FLAGS)

#
# Main programs
#

topolink.o : Makefile topolink.f90
	$(FORTRAN) -c topolink.f90 $(FLAGS)

linkcorrelation.o : Makefile linkcorrelation.f90 
	$(FORTRAN) -c linkcorrelation.f90 $(FLAGS)

filtermodels.o : Makefile filtermodels.f90 
	$(FORTRAN) -c filtermodels.f90 $(FLAGS)

maxmingyr.o : Makefile maxmingyr.f90 
	$(FORTRAN) -c maxmingyr.f90 $(FLAGS)

modelmaps.o : Makefile modelmaps.f90 
	$(FORTRAN) -c modelmaps.f90 $(FLAGS)

pbias.o : Makefile pbias.f90 
	$(FORTRAN) -c pbias.f90 $(FLAGS)

evalmodels.o : Makefile evalmodels.f90 
	$(FORTRAN) -c evalmodels.f90 $(FLAGS)

#
# modules
#

modules : $(MODULES) 

functionpars.o : functionpars.f90
	$(FORTRAN) -c functionpars.f90 $(FLAGS)

linkedcells.o : linkedcells.f90
	$(FORTRAN) -c linkedcells.f90 $(FLAGS)

topolink_data.o : topolink_data.f90
	$(FORTRAN) -c topolink_data.f90 $(FLAGS)

topolink_operations.o : topolink_operations.f90
	$(FORTRAN) -c topolink_operations.f90 $(FLAGS)

ioformat.o : ioformat.f90
	$(FORTRAN) -c ioformat.f90 $(FLAGS)

flashmod.o : flashmod.f90
	$(FORTRAN) -c flashmod.f90 $(FLAGS)

#
# All objects
#

cgnewton.o : cgnewton.f90 
	$(FORTRAN) -c cgnewton.f90 $(FLAGS)

computedpath.o : computedpath.f90
	$(FORTRAN) -c computedpath.f90 $(FLAGS)

computef.o : computef.f90
	$(FORTRAN) -c computef.f90 $(FLAGS)

computeg.o : computeg.f90
	$(FORTRAN) -c computeg.f90 $(FLAGS)

goverlapcell.o : goverlapcell.f90
	$(FORTRAN) -c goverlapcell.f90 $(FLAGS)

goverlap.o : goverlap.f90
	$(FORTRAN) -c goverlap.f90 $(FLAGS)

gstretch.o : gstretch.f90 
	$(FORTRAN) -c gstretch.f90 $(FLAGS)

initcells.o : initcells.f90
	$(FORTRAN) -c initcells.f90 $(FLAGS)

initguess.o : initguess.f90
	$(FORTRAN) -c initguess.f90 $(FLAGS)

linkstatus.o : linkstatus.f90
	$(FORTRAN) -c linkstatus.f90 $(FLAGS)

linkconsistency.o : linkconsistency.f90
	$(FORTRAN) -c linkconsistency.f90 $(FLAGS)

keywords.o : keywords.f90
	$(FORTRAN) -c keywords.f90 $(FLAGS)

lnf.o : lnf.f90
	$(FORTRAN) -c lnf.f90 $(FLAGS)

overlapcell.o : overlapcell.f90
	$(FORTRAN) -c overlapcell.f90 $(FLAGS)

overlap.o : overlap.f90
	$(FORTRAN) -c overlap.f90 $(FLAGS)

printdata.o : printdata.f90
	$(FORTRAN) -c printdata.f90 $(FLAGS)

random.o : random.f90
	$(FORTRAN) -c random.f90 $(FLAGS)

stretch.o : stretch.f90
	$(FORTRAN) -c stretch.f90 $(FLAGS)

flashsort.o : flashsort.f90 flashmod.o
	$(FORTRAN) -c flashsort.f90 $(FLAGS)

test_grad.o : test_grad.f90
	$(FORTRAN) -c test_grad.f90 $(FLAGS)

title.o : title.f90
	$(FORTRAN) -c title.f90 $(FLAGS)

bindir : 
	@mkdir -p $(BIN)

clean : 
	rm -f ./*.o ./*.mod $(BIN)/*




