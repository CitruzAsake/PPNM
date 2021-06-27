CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
LDLIBS += $(shell gsl-config --libs)
CFLAGS += $(shell gsl-config --cflags)


all: output.txt errorscale.pdf
	@echo
	@cat output.txt

errorscale.pdf: errorscale.ppl
	pyxplot $< 1> /dev/null

errorscale.ppl: errorscale.txt Makefile
	@echo 'set terminal pdf dpi 800'						>  $@
	@echo 'set output "errorscale.pdf"'						>> $@
	@echo 'set key outside'								>> $@
	@echo 'set grid'								>> $@
	@echo 'set xlabel "$$\ln(N_values)$$"'							>> $@
	@echo 'set ylabel "$$\ln$$(error_values)"'						>> $@
	@echo 'set title "Error for $$\int_{0}^{1}$$ d$$x"'
	>> $@
	@echo 'plot \\' 								>> $@
	@echo '"errorscale.txt" using (ln($$1)):(ln($$2)) w l lw 2 color Blue title "Values for Psuedo" 				\\'>> $@
	@echo ',"errorscale.txt" using (ln($$1)):(ln($$3)) w l lw 2 color Black  title "Values for Quasi" 	\\'>> $@

output.txt errorscale.txt: main
	./$< > $@

main: main.o plainMC.o pseudoMC.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf
