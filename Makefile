#
# --- Jayant
#
F90=gfortran
F90FLAGS=-g -O0 -fbacktrace
# -fcheck=all -Wall

#F90=pgf90
#F90FLAGS=-g -O2 -traceback


EXEC=./bin/cldMicDropAct.$(F90)

LDFLAGS=

OBJ =	mo_fio_utils.o \
	Micro_HugMorr.o \
	aerosolnact.o \
	main.o

main:	$(OBJ)
	$(F90) -o $(EXEC) $(LDFLAGS) $(OBJ)

mo_fio_utils.o	:	src/mo_fio_utils.f90
	$(F90) -c $(F90FLAGS) $(LDFLAGS) src/mo_fio_utils.f90

Micro_HugMorr.o	:	src/Micro_HugMorr.f90
	$(F90) -c $(F90FLAGS) $(LDFLAGS) src/Micro_HugMorr.f90

aerosolnact.o	:	src/aerosolnact.f90
	$(F90) -c $(F90FLAGS) $(LDFLAGS) src/aerosolnact.f90

main.o	:	src/main.f90 mo_fio_utils.o Micro_HugMorr.o
	$(F90) -c $(F90FLAGS) $(LDFLAGS) src/main.f90

clean:
	-rm $(OBJ)
	-rm $(EXEC)
	-rm *.mod
