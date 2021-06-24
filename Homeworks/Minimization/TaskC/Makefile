CC = clang
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: output.txt
	@echo
	@echo '              TASK C                   '
	@echo '-----------------------------------------'
	@cat  output.txt 

output.txt: main
	./$< > $@

main: main.o simplexmethod.c utilities.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png utilities.o

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)

