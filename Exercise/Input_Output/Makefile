CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm


output.txt error.txt outputfiles.txt: input.txt cmdline stdio file 
	./cmdline $$(cat input.txt) 1>output.txt 2>error.txt
	./stdio < input.txt > out.file.txt
	./file input.txt outputfiles.txt

input.txt: Makefile 
	cat /dev/null > $@
	for x in $$(seq 0 1 5); do echo "$$x" >> $@; done

clean: 
	$(RM) input.txt out* error.txt out.file.txt outputfiles.txt cmdline stdio file
