
############# MathLink related paths and settings start: ####################
#### Specify here the fulll path to the Mathematica MathLink Developer Kit:
MLINKDIR = /export/local/Wolfram/Mathematica/6.0/SystemFiles/Links/MathLink/DeveloperKit/

#### Specify the system SYS and the library name MLIB. 
### For 64-bit Linux on x86 platform:
SYS = Linux-x86-64
##Note, in Mathematica 6.0 the name for Linux-x86-64 was changed from ML to ML64i3:
MLIB = -lML64i3

### For 32-bit Linux on x86 platform:
#SYS = Linux
#MLIB = -lML
############# MathLink related paths and settings end. ####################

# You must specify either icc, or gcc.
############# icc related paths and settings start: ####################
#CC = icc
#FORT = ifort
#LIBS = -pthread  -Wl,-Bstatic ${MLIB} -lifport -lifcore -ldl -lifport -lrt -Wl,-Bdynamic
#LIBDIRS = -L/usr/local/lib  -L/export/local/intel/fce/9.0/lib
#LFLAGS = -i-static
############# icc related paths and settings end ####################
############# gcc related paths and settings start: ####################
# You must specify either g77, or gfortran, with the corresponding LIBS.
CC=gcc
FORT=g77
### LIBS for g77:
LIBS = -pthread -Wl,-Bstatic ${MLIB} -lg2c -lm -lrt -Wl,-Bdynamic
#FORT=gfortran
### LIBS for gfortran:
#LIBS = -pthread -Wl,-Bstatic ${MLIB} -lgfortranbegin -lgfortran -lm -lrt -Wl,-Bdynamic
LIBDIRS = -L/usr/local/lib

############# gcc related paths and settings end ####################

CFLAGS= -O3

###### Do not change anything below this line.


CADDSDIR = ${MLINKDIR}/${SYS}/CompilerAdditions
INCDIR = ${CADDSDIR}
MLIBDIR = ${CADDSDIR}
MPREP = ${CADDSDIR}/mprep

OBJ = vegas.o collectstr.o scanner.o  runline.o Link.o hash.o CIntegrate.o

CIntegrate$ : $(OBJ)
	$(CC) ${LFLAGS} $(OBJ) -L${MLIBDIR} ${LIBDIRS} ${LIBS} -o $@

Link.c : Link.tm
	${MPREP} -o Link.c Link.tm

vegas.o: vegas.f
	$(FORT) -c vegas.f

.c.o:
	$(CC) $(CFLAGS) -I${INCDIR} -c $<


# ******** DEPENDENCIES: *********

CIntegrate.o: scanner.h
CIntegrate.o: comdef.h
CIntegrate.o: collectstr.h
CIntegrate.o: hash.h
CIntegrate.o: runline.h
CIntegrate.o: vegas.h
collectstr.o: collectstr.h
collectstr.o: comdef.h
hash.o: comdef.h
hash.o: scanner.h
hash.o: collectstr.h
hash.o: hash.h
runline.o: scanner.h
runline.o: comdef.h
runline.o: collectstr.h
runline.o: hash.h
runline.o: runline.h
scanner.o: scanner.h
scanner.o: comdef.h
scanner.o: collectstr.h
scanner.o: hash.h
scanner.o: queue.h
