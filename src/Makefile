##
## Plane wave makefile 
## Jan Dettmer jand@uvic.ca 2007
##
## With full optimization KISPIOX
#FC      = mpif90 -Kieee -C -mp -Mchkfpstk -traceback
FC       = mpif90 -fast -O4 -Kieee -tp=sandybridge #-C -Mchkfpstk -traceback

## With full optimization HOMATHKO
#FC      = mpif90 -Kieee
## With full optimization STIKINE
#FC      = mpif90 -O3 -vec-report0 -axS -xS -opt-mem-bandwidth2 -ip -std03

## Without optimization
#FC       = mpif90 -Kieee -O0 -g -gopt -Mbounds -traceback
#FC      = mpif90 -O1 -vec-report0 -std03 -debug
## Gfortran 
#DB       = -Wall -Wextra -ffpe-trap=zero,overflow,underflow -fbounds-check \
	   -fbacktrace -O0 -g #Warning flags
#FC       = mpif90 $(DB)


##############################################################################
##
##  Define binary, object and include paths:
##
BIN = ./bin
OBJ = ./obj
INC = include
##############################################################################

all: bin/psrjmh_ais_auto

psrjmhreais:  bin/psrjmh_reais_auto

psrjmhais:  bin/psrjmh_ais_auto

clean: clean.done
	rm -f *.done
clean.done:
	rm -f bin/*
	rm -f *.o *.mod $(OBJ)/*
	rm -f $(INC)/*.mod
	rm -f make.log
##############################################################################
#
# Targets for linking of programs:
#
# PSRJMHREAIS:
#
$(BIN)/psrjmh_reais_auto: $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o $(OBJ)/svdcmp.o \
                      $(OBJ)/pythag.o $(OBJ)/median.o $(OBJ)/psrjmh_reais_auto.o
	$(FC)  -o $(BIN)/psrjmh_reais_auto $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o \
                  $(OBJ)/psrjmh_reais_auto.o $(OBJ)/svdcmp.o $(OBJ)/pythag.o $(OBJ)/median.o

$(OBJ)/psrjmh_reais_auto.o: psrjmh_reais_auto.f90 
	$(FC) -c -o $(OBJ)/psrjmh_reais_auto.o psrjmh_reais_auto.f90 -I$(INC)

##############################################################################
#
# PSRJMHAIS:
#
$(BIN)/psrjmh_ais_auto: $(OBJ)/data_type.o $(OBJ)/rjmcmc_com.o $(OBJ)/utils.o \
	                $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o $(OBJ)/svdcmp.o \
                        $(OBJ)/pythag.o $(OBJ)/median.o $(OBJ)/psrjmh_ais_auto.o
	$(FC)  -o $(BIN)/psrjmh_ais_auto $(OBJ)/data_type.o $(OBJ)/rjmcmc_com.o $(OBJ)/utils.o \
		  $(OBJ)/nrtype.o $(OBJ)/nr.o $(OBJ)/nrutil.o \
                  $(OBJ)/svdcmp.o $(OBJ)/pythag.o $(OBJ)/median.o \
	          -I$(INC) $(OBJ)/psrjmh_ais_auto.o

$(OBJ)/psrjmh_ais_auto.o: psrjmh_ais_auto.f90 
	$(FC) -c -o $(OBJ)/psrjmh_ais_auto.o psrjmh_ais_auto.f90 -I$(INC)

##############################################################################
#
# Targets for separate modules:
$(OBJ)/data_type.o: data_type.f90
	$(FC) -c -o $(OBJ)/data_type.o data_type.f90
	mv data_type.mod $(INC) 
$(OBJ)/rjmcmc_com.o: rjmcmc_com.f90 $(OBJ)/data_type.o $(INC)/data_type.mod
	$(FC) -c -o $(OBJ)/rjmcmc_com.o rjmcmc_com.f90 -I$(INC)
	mv rjmcmc_com.mod $(INC)
$(OBJ)/utils.o: utils.f90 $(OBJ)/data_type.o $(INC)/data_type.mod $(OBJ)/rjmcmc_com.o $(INC)/rjmcmc_com.mod
	$(FC) -c -o $(OBJ)/utils.o utils.f90 -I$(INC)
	mv utils.mod $(INC)
##############################################################################
#
# Targets for headers:
$(INC)/data_type.mod: data_type.f90
	$(FC) -c -o $(OBJ)/data_type.o data_type.f90
	mv data_type.mod $(INC) 
$(INC)/rjmcmc_com.mod: rjmcmc_com.f90 $(OBJ)/data_type.o $(INC)/data_type.mod
	$(FC) -c -o $(OBJ)/rjmcmc_com.o rjmcmc_com.f90 -I$(INC)
	mv rjmcmc_com.mod $(INC)
$(INC)/utils.mod: utils.f90 $(OBJ)/data_type.o $(INC)/data_type.mod $(OBJ)/rjmcmc_com.o $(INC)/rjmcmc_com.mod
	$(FC) -c -o $(OBJ)/utils.o utils.f90 -I$(INC)
	mv utils.mod $(INC)
##############################################################################


##############################################################################
#
# Targets for separate subroutines:
#
#$(OBJ)/slice.o: slice.f90
#	$(FC) -c -o $(OBJ)/slice.o slice.f90

$(OBJ)/median.o: median.f90
	$(FC) -c -o $(OBJ)/median.o median.f90
	mv m_valmed.mod $(INC)
##############################################################################
#
# Targets to compile NR subroutines that are in separate files:
#
$(OBJ)/ludcmp.o: ludcmp.f90
	$(FC) -c -o $(OBJ)/ludcmp.o ludcmp.f90 -I$(INC)

$(OBJ)/svdcmp.o: svdcmp.f90
	$(FC) -c -o $(OBJ)/svdcmp.o svdcmp.f90 -I$(INC)

$(OBJ)/pythag.o: pythag.f90
	$(FC) -c -o $(OBJ)/pythag.o pythag.f90 -I$(INC)

##############################################################################
#
#  Targets for numerical recipes modules:
#
$(OBJ)/nrtype.o: nrtype.f90 
	$(FC) -c -o $(OBJ)/nrtype.o nrtype.f90
	mv nrtype.mod $(INC)

$(OBJ)/nr.o: nr.f90 $(OBJ)/nrtype.o
	$(FC) -c -o $(OBJ)/nr.o nr.f90 -I$(INC)
	mv nr.mod $(INC)

$(OBJ)/nrutil.o: nrutil.f90 $(OBJ)/nrtype.o
	$(FC) -c -o $(OBJ)/nrutil.o nrutil.f90 -I$(INC)
	mv nrutil.mod $(INC)

##############################################################################
# ...this is the end my fiend
#
# EOF
