CC= gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

out.text: math
	./math > out.txt

math: math.c
	$(CC) $(LDFLAGS) $^ -o $@ $(LBLIBS)

clean:
	$(RM) *.o out* math
