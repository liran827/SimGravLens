CC=gcc -std=c99 

CFLAGS=-fPIC -O3
FFTINC = -I/usr/local/include 
FFTLIB = -L/usr/local/lib -L/usr/local/Cellar/gsl/2.4/lib/
#GSLINC = -I/usr/local/include/gsl

LIBS = -lm -lfftw3 -lgsl -lgslcblas
OPT = $(FFTINC) $(FFTLIB) #$(GSLINC)

lens.so:lens.o convolution.o
	$(CC) -shared  -o liblens.so  $(CFLAGS) lens.o convolution.o $(LIBS) $(OPT)

lens.o: lens.c lens.h
	$(CC) -o lens.o $(CFLAGS) -c lens.c

convolution.o: convolution.h
	$(CC) -o convolution.o $(CFLAGS) $(OPT) -c convolution.c 

iLens: main.o  icosmo.o  lens.o convolution.o
	$(CC) -o iLens $(CFLAGS) main.o lens.o convolution.o  icosmo.o $(LIBS) $(OPT)

main:
	$(CC) -o main.o $(CFLAGS) $(OPT) $(LIBS) -c main.c

icosmo: icosmo.h
	$(CC) -o icosmo.o $(CFLAGS) -c icosmo.c  $(OPT) $(LIBS) 


.PHONY: clean
clean:
	rm -f *.o *.so *.out *.fits

