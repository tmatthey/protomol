OFILES= registerForceExemplarsBonded.o registerForceExemplars.o registerForceExemplarsCutoff.o registerForceExemplarsFastElectrostatic.o registerForceExemplarsFull.o registerForceExemplarsOther.o registerForceExemplarsSimpleFull.o registerIntegratorExemplars.o registerOutputExemplars.o registerTopologyExemplars.o protomol.o

INCLUDES= -I../../framework/frontend -I../../framework/base -I../../framework/forces -I../../framework/imd -I../../framework/integrators -I../../framework/io -I../../framework/parallel -I../../framework/topology  -I../../applications/liteprotomol-app -I../../

FLAGS= -DLITE_SWITCH -DLITE_PERIODIC

CFLAGS = -DHAVE_CONFIG_H -Wall -O9 -ffast-math -finline-functions -funroll-loops -DNDEBUG ${INCLUDES} ${FLAGS} -c -o 

protomol: ${OFILES}
	g++  -Wall -O9 -ffast-math -finline-functions -funroll-loops -DNDEBUG   -o protomol ${OFILES} -L../../framework/base -L../../framework/io -L../../framework/parallel -L../../framework/topology -L../../framework/forces -L../../framework/imd -L../../framework/integrators -L../../framework/frontend -L../../framework/factories -lfrontend -lintegrators -limd -lforces -ltopology -lparallel -lio -lbase  -llapack -lblas -lg2c -lpng  -lm -lnsl

protomol.o: ../../applications/protomol-app/protomol.cpp
	g++ ${CFLAGS} protomol.o ../../applications/protomol-app/protomol.cpp

registerForceExemplarsBonded.o: registerForceExemplarsBonded.cpp
	g++ ${CFLAGS} registerForceExemplarsBonded.o registerForceExemplarsBonded.cpp 

registerForceExemplars.o: registerForceExemplars.cpp
	g++ ${CFLAGS} registerForceExemplars.o registerForceExemplars.cpp

registerForceExemplarsCutoff.o: registerForceExemplarsCutoff.cpp
	g++ ${CFLAGS} registerForceExemplarsCutoff.o registerForceExemplarsCutoff.cpp

registerForceExemplarsFastElectrostatic.o: registerForceExemplarsFastElectrostatic.cpp
	g++ ${CFLAGS} registerForceExemplarsFastElectrostatic.o registerForceExemplarsFastElectrostatic.cpp

registerForceExemplarsFull.o: registerForceExemplarsFull.cpp
	g++ ${CFLAGS} registerForceExemplarsFull.o registerForceExemplarsFull.cpp

registerForceExemplarsOther.o: registerForceExemplarsOther.cpp
	g++ ${CFLAGS} registerForceExemplarsOther.o registerForceExemplarsOther.cpp 

registerForceExemplarsSimpleFull.o: registerForceExemplarsSimpleFull.cpp
	g++ ${CFLAGS} registerForceExemplarsSimpleFull.o registerForceExemplarsSimpleFull.cpp

registerIntegratorExemplars.o: registerIntegratorExemplars.cpp
	g++ ${CFLAGS} registerIntegratorExemplars.o registerIntegratorExemplars.cpp

registerOutputExemplars.o: registerOutputExemplars.cpp
	g++ ${CFLAGS} registerOutputExemplars.o registerOutputExemplars.cpp

registerTopologyExemplars.o: registerTopologyExemplars.cpp
	g++ ${CFLAGS} registerTopologyExemplars.o registerTopologyExemplars.cpp

clean:
	/bin/rm *.o protomol