CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: output.txt
	@echo
	@cat output.txt
	@echo 
	@echo 'Number 2 and number 5 results comparison to see the effect of the Clenshaw-Curtis transformation.'
	@echo 'The most accurate is the result with the Clenshaw-Curtis transformation'
	@echo 'Though compared to the GSL integration routines the accuracy is smaller'
output.txt: main
	./$< > $@

main: main.o RAI.o CC.o

clean:
	$(RM) main *.o *.txt  
