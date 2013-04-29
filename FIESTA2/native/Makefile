vpath %.c ../sources
vpath %.h ../sources

CIntegrate : $(OBJ)
	$(CC) ${LFLAGS} $(OBJ) ${MLIBDIR} ${LIBDIRS} ${LIBS} -o $@

.c.o:
	$(CC) $(CFLAGS) -DARITHMETIC=NATIVE ${INCDIR} -c $<

dep: depend

depend:
	cd ../sources && $(MAKEDEPENDNAME) -f ../native/Makefile $(OBJ:.o=.c)

# ******** DEPENDENCIES: *********

collectstr.o: collectstr.h
collectstr.o: comdef.h
scanner.o: scanner.h
scanner.o: comdef.h
scanner.o: triesi.h
scanner.o: collectstr.h
scanner.o: hash.h
scanner.o: tries.h
scanner.o: constants.h
scanner.o: queue.h
runline.o: scanner.h
runline.o: comdef.h
runline.o: triesi.h
runline.o: collectstr.h
runline.o: hash.h
runline.o: runline.h
hash.o: comdef.h
hash.o: scanner.h
hash.o: triesi.h
hash.o: collectstr.h
hash.o: hash.h
integrators.o: integrators.h
integrators.o: scanner.h
integrators.o: comdef.h
integrators.o: triesi.h
integrators.o: collectstr.h
integrators.o: hash.h
integrators.o: runline.h
CIntegrate.o: scanner.h
CIntegrate.o: comdef.h
CIntegrate.o: triesi.h
CIntegrate.o: collectstr.h
CIntegrate.o: hash.h
CIntegrate.o: runline.h
CIntegrate.o: integrators.h
