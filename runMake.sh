#!/bin/sh

export pdir=.
export inpdir=${pdir}/input
export bindir=${pdir}/bin

#export F90=pgf90
export F90=gfortran

export EXEC=./bin/cldMicDropAct.${F90}

. ${pdir}/env.setup

MAKE=make

#${MAKE} clean
${MAKE}

${EXEC}
