CC= gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

all: out.txt
	cat out.txt

out.txt: math
	./math > out.txt

math.o: math.c
	$(CC) $(CFLAGS) -c math.c -o math.o
math: math.c
	$(CC) $(LDFLAGS) $^ -o $@ $(LBLIBS)

clean:
	$(RM) *.o out* math
