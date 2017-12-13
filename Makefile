LIB        = -L. 
INCLUDE    = -I.
CFLAGS     = -O3
EXEC       = QF.exe
CXX        = g++

${EXEC}: QF.cpp
	${CXX} ${CFLAGS} ${INCLUDE} ${LIB} QF.cpp -o ${EXEC}

clean:
	rm -f *.o

%.o: $.cpp
	${CXX} -c ${CFLAGS} ${INCL} -cpp -o $*.o $<
