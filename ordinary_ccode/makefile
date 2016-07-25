NAULINK=nauty/gtools.o nauty/nauty1.o nauty/nautil1.o nauty/naugraph1.o nauty/naugroup.o
CC=gcc
CPP=g++
CFLAGS= -O4 
CPPFLAGS=

all: gradiff gradiff_hairy liegradiff geng2 makemat makemat_hairy brackettest genwheelgra findgraindex mattohg6 some compbracket filter;

geng.o: nauty/gtools.h nauty/geng.c
	${CC} ${CFLAGS} -c -o geng.o -DMAXN=32 -DOUTPROC=myoutproc -DGENG_MAIN=geng_main2 -DPRUNE=mypruner -DPREPRUNE=myprepruner  \
		nauty/geng.c

tgra.o: tgra.h tgra.C nauty/gtools.h
	${CPP} ${CPPFLAGS} -c -DMAXN=32 -o tgra.o tgra.C

gradiff: nauty/gtools.h gradiff.C geng.o tgra.o nauty/gtools.o nauty/nauty1.o nauty/nautil1.o nauty/naugraph1.o
	${CPP} ${CPPFLAGS} -o gradiff -DMAXN=32 -DOUTPROC=myoutproc -DGENG_MAIN=geng_main2 -DPRUNE=mypruner -DPREPRUNE=myprepruner -Wno-write-strings \
		gradiff.C geng.o  tgra.o ${NAULINK}

gradiff_hairy: nauty/gtools.h gradiff_hairy.C geng.o tgra.o nauty/gtools.o nauty/nauty1.o nauty/nautil1.o nauty/naugraph1.o
	${CPP} ${CPPFLAGS} -o gradiff_hairy -DMAXN=32 -DOUTPROC=myoutproc -DGENG_MAIN=geng_main2 -DPRUNE=mypruner -DPREPRUNE=myprepruner -Wno-write-strings \
		gradiff_hairy.C geng.o  tgra.o ${NAULINK}

liegradiff: nauty/gtools.h liegradiff.C geng.o tgra.o nauty/gtools.o nauty/nauty1.o nauty/nautil1.o nauty/naugraph1.o
	${CPP} ${CPPFLAGS} -o liegradiff -DMAXN=32 -DOUTPROC=myoutproc -DGENG_MAIN=geng_main2 -DPRUNE=mypruner -DPREPRUNE=myprepruner -Wno-write-strings \
		liegradiff.C geng.o  tgra.o ${NAULINK}

geng2: nauty/gtools.h geng2.C geng.o tgra.o nauty/gtools.o nauty/nauty1.o nauty/nautil1.o nauty/naugraph1.o
	${CPP} ${CPPFLAGS} -o geng2 -DMAXN=32 -DOUTPROC=myoutproc -DGENG_MAIN=geng_main2 -DPRUNE=mypruner -DPREPRUNE=myprepruner -Wno-write-strings \
		geng2.C geng.o  tgra.o ${NAULINK}
		
makemat: makemat.C
	${CPP} ${CPPFLAGS} -o makemat makemat.C
	
makemat_hairy: makemat_hairy.C
	${CPP} ${CPPFLAGS} -o makemat_hairy makemat_hairy.C

brackettest: brackettest.C tgra.o
	${CPP} ${CPPFLAGS} -o brackettest -DMAXN=32 brackettest.C ${NAULINK} tgra.o

genwheelgra: genwheelgra.C tgra.o 
	${CPP} ${CPPFLAGS} -o genwheelgra -DMAXN=32 genwheelgra.C ${NAULINK} tgra.o

findgraindex: findgraindex.C tgra.o 
	${CPP} ${CPPFLAGS} -o findgraindex -DMAXN=32 findgraindex.C ${NAULINK} tgra.o

mattohg6: mattohg6.C tgra.o 
	${CPP} ${CPPFLAGS} -o mattohg6 -DMAXN=32 mattohg6.C ${NAULINK} tgra.o

some: some.C tgra.o
	${CPP} ${CPPFLAGS} -o some -DMAXN=32 -Wno-write-strings some.C ${NAULINK} tgra.o
	
compbracket: compbracket.C tgra.o
	${CPP} ${CPPFLAGS} -o compbracket -DMAXN=32 -Wno-write-strings compbracket.C ${NAULINK} tgra.o
	
filter: filter.C tgra.o
	${CPP} ${CPPFLAGS} -o filter -DMAXN=32 -Wno-write-strings filter.C ${NAULINK} tgra.o

clean:
	rm -f *.o

veryclean:
	 rm -f *.g6 *.hg6 *.txt
