CFLAGS= -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)
CC = gcc

.PHONEY:default
default: out.myfunc.txt

out.myfunc.txt:myfunc
	./$< > $@

.PHONEY:clean
clean:
	$(RM) out* myfunc

.PHONEY:test
test:
	echo FLAGS
	echo $(CFLAGS)
	echo $(LDLIBS)
