############# MathLink related paths and settings start: ####################
#### Specify here the fulll path to the Mathematica MathLink Developer Kit:
MLINKDIR = C:\Program Files\Wolfram Research\Mathematica\6.0\SystemFiles\Links\MathLink\DeveloperKit
MATHLINK = $(MLINKDIR)\$(SYS)\CompilerAdditions\mldev$(BIT)
MLIBDIR  = "$(MATHLINK)\lib"
MINCDIR = "$(MATHLINK)\Include"

#### Specify the system SYS
## for the 64bit platform:
## SYS = Windows-x86-64
## BIT = 64
## for the 32bit platform:
SYS = Windows
BIT = 32

### Specify the mathlink library name
## default setting
MLLIB = ml$(BIT)i3m.lib
## in case you compile with the mathlink that came with Mathematica 5.2 (only 32-bit version)
## MLLIB = ml32i2c.lib
############# MathLink related paths and settings end. ####################

###### Please specify the path to the software developer kit ########
SDK = "C:\Program Files\Microsoft Platform SDK for Windows Server 2003 R2"
LIBDIR = $(SDK)\Lib
INCDIR = $(SDK)\Include

##### The names of the compilers and linkers
## The mathlink compiler tm -> c
MPREP    = "$(MATHLINK)\bin\mprep"
## the main compiler and related settings
CL = cl
CFLAGS =  /TC /nologo /O2  -I$(INCDIR) -I$(MINCDIR) 
## the fortran compiler
IFORT_COMPILER_DIR = C:\Program Files\Intel\Compiler\Fortran\10.1.011
FORTRAN = "$(IFORT_COMPILER_DIR)\ia$(BIT)\bin\ifort"
FORTRAN_LIB_DIR = "$(IFORT_COMPILER_DIR)\ia$(BIT)\lib"
## the linker and related settings
LINK = link
LFLAGS = /NOLOGO /SUBSYSTEM:Console /INCREMENTAL:no /LIBPATH:$(LIBDIR) /LIBPATH:$(MLIBDIR) /LIBPATH:$(FORTRAN_LIB_DIR)




##############################################################################
###### Do not change anything below this line.

LIBS = $(MLLIB) user$(BIT).lib kernel$(BIT).lib gdi$(BIT).lib
 
OBJECTS = Link.obj CIntegrate.obj collectstr.obj runline.obj scanner.obj vegas.obj hash.obj

CIntegrate.exe : $(OBJECTS)
    $(LINK) $(LFLAGS) $(LIBS) $(OBJECTS) /OUT:CIntegrate.exe

Link.c : Link.tm
    $(MPREP) -o Link.c Link.tm 

vegas.obj : vegas.f
    $(FORTRAN) -c vegas.f 

.c.obj:
    $(CL) /c $(CFLAGS) $<

# ******** DEPENDENCIES: *********
CIntegrate.obj: scanner.h
CIntegrate.obj: comdef.h
CIntegrate.obj: collectstr.h
CIntegrate.obj: hash.h
CIntegrate.obj: runline.h
CIntegrate.obj: vegas.h
collectstr.obj: collectstr.h
collectstr.obj: comdef.h
hash.obj: comdef.h
hash.obj: scanner.h
hash.obj: collectstr.h
hash.obj: hash.h
runline.obj: scanner.h
runline.obj: comdef.h
runline.obj: collectstr.h
runline.obj: hash.h
runline.obj: runline.h
scanner.obj: scanner.h
scanner.obj: comdef.h
scanner.obj: collectstr.h
scanner.obj: hash.h
scanner.obj: queue.h
