CFLAGS= -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

out.myfunc.txt:myfunc
	./$< >$@ 


.PHONEY:clean
clean:
	$(RM) out* myfunc

