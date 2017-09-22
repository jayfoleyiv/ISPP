LIB        = -L. 
INCLUDE    = -I.
CFLAGS     = -O3
EXEC       = ISPP.exe
CXX        = g++

${EXEC}: ISPP.cpp
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} ISPP.cpp -o ${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<
