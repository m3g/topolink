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
#   FLAGS = -fbounds-check -Wunused
#   FLAGS = -fbounds-check
#   FLAGS = -g -pg 
endif

# Directory for executables
BIN=../bin
 
###################################################################
#                                                                 #
# Generally no modifications are required after this.             #
#                                                                 #
###################################################################

all : bindir $(BIN)/topolink

$(BIN)/topolink: Makefile \
                 functionpars.mod \
                 linkedcells.mod \
                 topolink_data.mod \
                 topolink_operations.mod \
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
                          cgnewton.f90 \
                          computedpath.f90 \
                          computef.f90 \
                          computeg.f90 \
                          functionpars.f90 \
                          topolink_data.f90 \
                          topolink_operations.f90 \
                          goverlapcell.f90 \
                          goverlap.f90 \
                          gstretch.f90 \
                          initcells.f90 \
                          initguess.f90 \
                          linkstatus.f90 \
                          keywords.f90 \
                          linkedcells.f90 \
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

functionpars.mod : functionpars.f90
	$(FORTRAN) -c functionpars.f90 $(FLAGS)

linkedcells.mod : linkedcells.f90
	$(FORTRAN) -c linkedcells.f90 $(FLAGS)

topolink_operations.mod : topolink_operations.f90 topolink_data.mod
	$(FORTRAN) -c topolink_operations.f90 $(FLAGS)

topolink_data.mod : topolink_data.f90
	$(FORTRAN) -c topolink_data.f90 $(FLAGS)

bindir : 
	@mkdir -p $(BIN)

clean : 
	rm -f ./*.o ./*.mod $(BIN)/*




