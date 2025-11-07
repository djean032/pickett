CXX=g++
CFLAGS=-g -Wall
#CFLAGS=-O2 -Wall
EXEQ=spfit spcat calmrg dpfit dpcat
EXEA=${EXEQ} iamcalc moiam stark termval sortn calbak
#next line for atlas blas
#BLASLIB=-lcblas -latlas
#next line for fortran blas and cblas wrappers
#BLASLIB=-lcblas -lblas
#next line for supplied routines
BLASLIB=
ifndef ($(BLASLIB))
	LBLAS=dblas.o
endif
default: ${EXEQ} 
all: ${EXEA}
test: spfit
	python3 ./tests/test_spfit.py
install:  
	-mv ${EXEA} /usr/local/bin 
	-chmod o+rx /usr/local/bin/*
dpfit: calfit.o subfit.o dpi.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
dpcat: calcat.o sortsub.o dpi.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
spfit: calfit.o subfit.o spinv.o spinit.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
spcat: calcat.o sortsub.o spinv.o spinit.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
calmrg: calmrg.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
calbak: calbak.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
termval: termval.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
stark: stark.o splib.a ; g++ -o $@ $^ $(BLASLIB) -lm
moiam: moiam.o ftran.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
iamcalc: iamcalc.o ftran.o splib.a; g++ -o $@ $^ $(BLASLIB) -lm
calbak: calbak.o splib.a ; g++ -o $@ $^ $(BLASLIB) -lm
cnvwn: cnvwn.o splib.a ; g++ -o $@ $^ $(BLASLIB) -lm

splib.a: ulib.o cnjj.o slibg++.o catutil.o lsqfit.o $(LBLAS)
	ar r splib.a $^
	ranlib splib.a

calfit.o:calfit.cpp calpgm.h
subfit.o:subfit.cpp calpgm.h
lsqfit.o:lsqfit.cpp lsqfit.h
calcat.o:calcat.cpp calpgm.h
sortsub.o: sortsub.cpp calpgm.h
calmrg.o:calmrg.cpp calpgm.h
termval.o:termval.cpp calpgm.h
stark.o:stark.cpp calpgm.h
iamcalc.o:iamcalc.cpp calpgm.h
moiam.o:moiam.cpp calpgm.h
ulib.o:ulib.cpp calpgm.h
cnjj.o:cnjj.cpp cnjj.h
slibg++.o:slibg++.cpp calpgm.h
spinv.o:spinv.cpp calpgm.h spinit.h
spinit.o:spinit.cpp calpgm.h spinit.h
dpi.o:dpi.cpp calpgm.h
ftran.o:ftran.cpp
sortn.o: sortn.cpp
sortsub.o: sortsub.cpp
dblas.o: dblas.cpp
util.o:util.cpp
calbak.o:calbak.cpp
cnvwn.o:cnvwn.cpp
sortn: sortn.o sortsub.o; g++ -o $@ $^
