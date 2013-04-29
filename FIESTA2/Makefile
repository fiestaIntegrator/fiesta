############# MathLink related paths and settings begin: ####################
#### Specify here the fulll path to the Mathematica MathLink Developer Kit:
#MLINKDIR = /home/MDK
#MLINKDIR = /export/local/Wolfram/Mathematica/7.0/SystemFiles/Links/MathLink/DeveloperKit
MLINKDIR = /export/pc/mathematica/7.0/SystemFiles/Links/MathLink/DeveloperKit

#### Specify the system SYS and the library name MLIB. 
### For 64-bit Linux on x86 platform:
#SYS = Linux-x86-64
#MLIB = -lML64i3

### For 32-bit Linux on x86 platform:
SYS = Linux
MLIB = -lML32i3

### For 32-bit Windows on x86 platform:
#SYS = Windows
#MLIB = -lML32i3
############# MathLink related paths and settings end. ####################

############# CUBA related paths and settings begin: ####################
#Set this two varables to full paths to "cuba.h" and to the CUBAL library
#file. If you have CUBA installed in your system, leave this variables empty:
CUBAINCLUDEDIR = -I/x1/asmirnov/Cuba32
CUBALIBDIR = -L/x1/asmirnov/Cuba32

#CUBAINCLUDEDIR = -I/cygdrive/d/Cuba/
#CUBALIBDIR = -L/cygdrive/d/Cuba/
############# CUBA related paths and settings end. ####################

############# MPFR related paths and settings begin: ####################
#Set this two varables to full paths to "cuba.h" and to the mpfr library
#file. If you have mpfr installed in your system, leave this variables empty:
MPFRINCLUDEDIR = -I/x1/asmirnov/mpfr32
MPFRLIBDIR = -L/x1/asmirnov/mpfr32/.libs
# ATTENTION! Here we assume that the GMP library is installed on the system.
# If you have non-standard GMP installation, edit the file  ./mpfr/Makefile and 
# add -I<path_to_gmp_includes> to  INCDIR and -L<path_to_gmp_libraries> to LIBDIRS
############# MPFR related paths and settings end. ####################

############# gcc related paths and settings begin: ####################
CC=gcc

### For Windows:
#LFLAGS = 
#LIBCSTATICDIR =
#MPFRLINKAGE = 
#LIBS = -lpthread  -lcuba ${MLIB} -lm -lmpfr -mwindows

### MPFR and Mathlink are statically linked (preferable for LINUX):
LFLAGS =
LIBCSTATICDIR=
MPFRLINKAGE = -Wl,-Bstatic -lmpfr -lgmp -Wl,-Bdynamic
LIBS = -pthread  -lcuba  -Wl,-Bstatic ${MLIB} -lstdc++ -Wl,-Bdynamic -lrt -lm

### Everything is dynamic.
### ATTENTION! For this case you need dynamical libraries MLIB, MPFR and GMP
### to be installed on the system or LD_LIBRARY_PATH to be set properly!:
#LFLAGS =
#LIBCSTATICDIR=
#MPFRLINKAGE = -lmpfr -lgmp
#LIBS = -pthread  -lcuba  ${MLIB} -lstdc++ -lrt -lm

LIBDIRS = ${LIBCSTATICDIR} ${CUBALIBDIR}
############# gcc related paths and settings end. ####################

CFLAGS= -O3 -Wall -Wno-char-subscripts
#CFLAGS= -g -Wall -Wno-char-subscripts

######################################################################
###### Do not change anything below this line. #######################
######################################################################

MAKEDEPENDNAME = makedepend -w20 -Y -s '\# ******** DEPENDENCIES: *********'

CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions
INCDIR = -I${CADDSDIR} ${CUBAINCLUDEDIR} 
MLIBDIR = -L${CADDSDIR}
MPREP = ${CADDSDIR}/mprep

ifeq ($(SYS),Windows)
CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions/cygwin
INCDIR = -I${CADDSDIR}/include ${CUBAINCLUDEDIR}
MLIBDIR = -L${CADDSDIR}/lib
MPREP = ${CADDSDIR}/bin/mprep.exe
endif

OBJ = collectstr.o scanner.o  runline.o Link.o hash.o integrators.o CIntegrate.o

export CFLAGS
export CC
export LFLAGS 
export LIBS
export LIBDIRS
export MAKEDEPENDNAME
export MLIBDIR
export OBJ
export INCDIR
export MPFRINCLUDEDIR
export MPFRLIBDIR
export MPFRLINKAGE

vpath %.c ./sources

all: Link.c
	$(MAKE) -C ./mpfr
	$(MAKE) -C ./native
	mv -f ./mpfr/CIntegrateMP .
	mv -f ./native/CIntegrate .

clean:
	rm -f ./mpfr/*.o
	rm -f ./native/*.o
	rm -f ./sources/Link.c

Link.c: Link.tm
	${MPREP} -o Link.c Link.tm
	mv Link.c ./sources

dep: depend

depend:
	$(MAKE) -C ./mpfr depend
	$(MAKE) -C ./native depend
