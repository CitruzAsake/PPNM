CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: tal
	./$< > $@

clean:
	$(RM) tal out.txt

