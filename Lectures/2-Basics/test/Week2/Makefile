CC = gcc
CFLAGS = -std=gnu99 -O -Wall
LDLIBS = -lm

out.txt: hello math
	./hello > out.txt
	./math >> out.txt
hello.o: hello.c
	$(CC) $(CFLAGS) -c hello.c -o hello.o

hello: hello.o world.o hello.h
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS) 

world.o: world.c hello hello.h
	$(CC) $(CFLAGS) -c $< -o $@

clean: 
	$(RM) *.o out* hello


