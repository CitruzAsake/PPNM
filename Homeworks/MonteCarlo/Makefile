CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: output.txt errorscaling.pdf
	@echo
	@cat output.txt

errorscaling.pdf: errorscaling.ppl
	pyxplot $< 1> /dev/null

errorscaling.ppl: errorscaling.txt Makefile
	@echo 'set terminal pdf dpi 800'						>  $@
	@echo 'set output "errorscaling.pdf"'						>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$\ln(N)$$"'								>> $@
	@echo 'set ylabel "$$\ln$$(Error)"'								>> $@
	@echo 'set title "Error Values for our problem"'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"errorscaling.txt" using (ln($$1)):(ln($$2)) w l lw 2 color Blue  title "Value for Psuedo" 				\\'>> $@
	@echo ',"errorscaling.txt" using (ln($$1)):(ln($$3)) w l lw 2 color Black  title "Values for Quasi" 	\\'>> $@

output.txt errorscaling.txt: main
	./$< > $@

main: main.o pMonteCarlo.o pseudoMC.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png 

