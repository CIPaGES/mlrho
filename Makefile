# Comment with icc
CC=gcc
CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 \
	-I/opt/local/include/ -L/opt/local/lib/   #-g  #-p  #-m64

# Comment with gcc
#CC=icc
#CFLAGS= -O3 -fast -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 \
	-I/opt/local/include/ -L/opt/local/lib/   #-g  #-p  #-m64

# The source files, object files, libraries and executable name.
SRCFILES= mlRho.c eprintf.c stringUtil.c interface.c mlComp.c piComp.c profile.c ld.c profileTree.c deltaComp.c
OBJFILES= mlRho.o eprintf.o stringUtil.o interface.o mlComp.o piComp.o profile.o ld.o profileTree.o deltaComp.o
LIBS= -lm -lgsl -lgslcblas -lz -lm

EXECFILE= mlRho
DIRECTORY= MlRho

VERSION= 2.0

# The make rule for the executable
.PHONY : all
all : $(EXECFILE)

$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)	

doc:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more

clean:
	rm -f *.o *~
realClean:
	rm -f *.o *~
	cd SimPro; make clean
tarfile:
	cd ../Doc; make clean; make pdf; cd ../$(DIRECTORY)_$(VERSION)
	mkdir $(DIRECTORY)_$(VERSION)
	cp -rf $(SRCFILES) valgrind.sh *.h  Scripts Makefile README COPYING \
	../Doc/mlRhoDoc.pdf $(DIRECTORY)_$(VERSION)
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -rf $(DIRECTORY)_$(VERSION)
