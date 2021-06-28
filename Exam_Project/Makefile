CC = gcc
CFLAGS = $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: all clean

all: out.txt TimeTest.txt

out.txt: main
	./$< > $@
	cat out.txt
	@echo
	cat TimeTest.txt
	@echo

main: main.o utilities.o jacobiSVD.o

clean:
	$(RM) main out.txt TimeTest.txt *.o
