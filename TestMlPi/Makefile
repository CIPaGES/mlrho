CC=gcc
CFLAGS= -O3 -Wall -Wshadow -pedantic -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64 -DVER32 -I/opt/local/include/ -L/opt/local/lib/#  -m64

# The source files, object files, libraries and executable name.
SRCFILES= testMlPi.c interface.c eprintf.c ranNum.c
OBJFILES= testMlPi.o interface.o eprintf.o ranNum.o

LIBS= -lm

EXECFILE= testMlPi
DIRECTORY= TestMlPi
VERSION= 0.4

# The make rule for the executable
.PHONY : all
all : $(EXECFILE)

$(EXECFILE) : $(OBJFILES)
	$(CC) $(CFLAGS) -o $(EXECFILE) $(OBJFILES) $(LIBS)	

# Other Standard make rules
lint : 
	lint $(SRCFILES) | more
clean:
	rm -f *.o *~
tarfile:
	mkdir $(DIRECTORY)_$(VERSION)
	cp $(SRCFILES) *.h Makefile $(DIRECTORY)_$(VERSION)
	tar cvzfh $(EXECFILE)_$(VERSION).tgz $(DIRECTORY)_$(VERSION)
	mv $(EXECFILE)_$(VERSION).tgz ../
	/bin/rm -r $(DIRECTORY)_$(VERSION)
