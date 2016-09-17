#
# Which compiler to use
#
COMPILER = GNU
#
#
ifeq "$(COMPILER)" "INTEL"
   FORTRAN = ifort
   FLAGS = -O3
endif
ifeq "$(COMPILER)" "GNU"
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
           $(BIN)/linkensemble \
           $(BIN)/linkcorrelation \
           $(BIN)/comparecols \
           $(BIN)/filtermodels

MODULES = functionpars.o \
          linkedcells.o \
          topolink_data.o \
          topolink_operations.o \
          ioformat.o \
          string_operations.o

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
                 strclean.o \
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
                 strclean.o \
                 lnf.o \
                 overlapcell.o \
                 overlap.o \
                 printdata.o \
                 random.o \
                 stretch.o \
                 test_grad.o \
                 title.o \
                 $(FLAGS)

$(BIN)/linkcorrelation : Makefile $(MODULES) linkcorrelation.o title.o
	$(FORTRAN) -o $(BIN)/linkcorrelation $(MODULES) linkcorrelation.o title.o $(FLAGS)

$(BIN)/filtermodels : Makefile $(MODULES) filtermodels.o title.o
	$(FORTRAN) -o $(BIN)/filtermodels $(MODULES) filtermodels.o title.o $(FLAGS)

$(BIN)/evalmodels : Makefile $(MODULES) evalmodels.o title.o sort.o model_index.o
	$(FORTRAN) -o $(BIN)/evalmodels $(MODULES) evalmodels.o title.o \
                                      sort.o model_index.o $(FLAGS)

$(BIN)/linkensemble : Makefile $(MODULES) linkensemble.o title.o sort.o model_index.o
	$(FORTRAN) -o $(BIN)/linkensemble $(MODULES) linkensemble.o title.o \
                                      sort.o model_index.o $(FLAGS)

$(BIN)/comparecols : Makefile $(MODULES) comparecols.o title.o
	$(FORTRAN) -o $(BIN)/comparecols $(MODULES) comparecols.o title.o $(FLAGS)


#
# Main programs
#

topolink.o : Makefile topolink.f90
	$(FORTRAN) -c topolink.f90 $(FLAGS)

linkcorrelation.o : Makefile linkcorrelation.f90
	$(FORTRAN) -c linkcorrelation.f90 $(FLAGS)

filtermodels.o : Makefile filtermodels.f90 
	$(FORTRAN) -c filtermodels.f90 $(FLAGS)

evalmodels.o : Makefile evalmodels.f90 
	$(FORTRAN) -c evalmodels.f90 $(FLAGS)

linkensemble.o : Makefile linkensemble.f90 
	$(FORTRAN) -c linkensemble.f90 $(FLAGS)

comparecols.o : Makefile comparecols.f90
	$(FORTRAN) -c comparecols.f90 $(FLAGS)

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

string_operations.o : string_operations.f90
	$(FORTRAN) -c string_operations.f90 $(FLAGS)

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

strclean.o : strclean.f90
	$(FORTRAN) -c strclean.f90 $(FLAGS)

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

test_grad.o : test_grad.f90
	$(FORTRAN) -c test_grad.f90 $(FLAGS)

sort.o : sort.f90
	$(FORTRAN) -c sort.f90 $(FLAGS)

model_index.o : model_index.f90
	$(FORTRAN) -c model_index.f90 $(FLAGS)

title.o : title.f90
	$(FORTRAN) -c title.f90 $(FLAGS)

bindir : 
	@mkdir -p $(BIN)

clean : 
	rm -f ./*.o ./*.mod $(BIN)/*




