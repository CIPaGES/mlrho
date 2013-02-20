CC=gcc
CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 \
	-I/opt/local/include/ -L/opt/local/lib/   -g  #-p  #-m64

# The source files, object files, libraries and executable name.
SRCFILES= mlRho.c eprintf.c stringUtil.c interface.c profileTree.c tab.c mlComp.c piComp.c \
	deltaComp.c rhoComp.c queue.c 
OBJFILES= mlRho.o eprintf.o stringUtil.o interface.o profileTree.o tab.o mlComp.o piComp.o \
	deltaComp.o rhoComp.o queue.o 
LIBS= -lm -lgsl -lgslcblas -lz -lm

EXECFILE= mlRho
DIRECTORY= MlRho

VERSION= 1.21

# The make rule for the executable
.PHONY : all
all : $(EXECFILE)

$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)	

test:	testPi testDelta
testPi:
	@cd TestMlPi; make
	@echo Output from testMlPi:
	@./TestMlPi/testMlPi -i 1 -d test.dat
	@echo Output from mlRho:
	@./mlRho -M 0 test.dat
testDelta:
	@cd TestMlDelta; make
	@echo Output from testMlDelta:
	@./TestMlDelta/testMlDelta -i 1 -d test.dat -s 10000 -l
	@echo Output from mlRho:
	@./mlRho -M 1 -l -f -T test.dat 

doc:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more

clean:
	rm -f *.o *~
realClean:
	rm -f *.o *~
	cd TestMlPi; make clean
	cd TestMlDelta; make clean

tarfile:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)
	cp -rf $(SRCFILES) valgrind.sh *.h  Scripts Makefile README COPYING \
	../Doc/mlRhoDoc.pdf $(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)/TestMlPi 
	cp TestMlPi/*.c TestMlPi/*.h TestMlPi/README TestMlPi/Makefile $(DIRECTORY)_$(VERSION)/TestMlPi
	mkdir $(DIRECTORY)_$(VERSION)/TestMlDelta
	cp TestMlDelta/*.c TestMlDelta/*.h TestMlDelta/Makefile $(DIRECTORY)_$(VERSION)/TestMlDelta
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -rf $(DIRECTORY)_$(VERSION)
