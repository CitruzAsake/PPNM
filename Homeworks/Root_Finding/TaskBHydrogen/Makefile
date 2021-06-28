CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)


Hydrogen.pdf: hydrogen.ppl
	pyxplot $< 1> /dev/null

hydrogen.ppl: hydrogen.txt Makefile
	@echo 'set terminal pdf dpi 800'						>  $@
	@echo 'set output "Hydrogen.pdf"'						>> $@
	@echo 'set key top' 					>> $@
	@echo 'set xlabel "$$Some X value$$"'								>> $@
	@echo 'set ylabel "$$Some Y value$$"'								>> $@
	@echo 'set title r"""Numerical integration Schr\"odinger equation"""'  >>$@
	@echo 'plot \\'				>> $@
	@echo '"hydrogen.txt" using 1:2 w l lt 1 lw 2 color Black   title "Numerical Data" 				\\'>> $@
	@echo ',"hydrogen.txt" using 1:3 w l lt 5 lw 2 color Yellow  title "Analytical" 	\\'>> $@

hydrogen.txt: main
	./$< > $@

main: main.o RK.o ../NewtonJacobi.o ../GramSchmidt.o ../utilities.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png
