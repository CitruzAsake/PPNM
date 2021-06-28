CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
LDLIBS += $(shell gsl-config --libs)
CFLAGS += $(shell gsl-config --cflags)

all: outfile.txt fitplot.pdf

	@echo '--------------------------------'
	@echo '-------------Task---------------'
	@echo '--------------------------------'
	@cat outfile.txt 
	@echo

fitplot.pdf: fitplot.ppl
	pyxplot $< 1> /dev/null

fitplot.ppl: outfile.txt Makefile

	@echo 'set terminal pdf dpi 800'						>  $@

	@echo 'set output "fitplot.pdf"'						>> $@

	@echo 'set key outside'								>> $@

	@echo 'set xlabel "$$x$$"'							>> $@

	@echo 'set ylabel "$$y$$"'							>> $@

	@echo 'set title "Artificial Neural Network different fits"'			>> $@

	@echo 'plot \\'									>> $@

	@echo '"outfile.txt"  index 1 using 1:1 w l lw 1 color black  title "Integral fit func" \\'>> $@

	@echo ',"outfile.txt" index 1 using 1:2 w l lw 1 color cyan title "Curve fit func" \\'>> $@

	@echo ',"outfile.txt" index 1 using 1:3 w l lw 1 color SpringGreen title "Derivative fit func" \\'>> $@

	@echo ',exp(-x*x)*(-5*sin(5*x-1)-2*x*cos(5*x-1)) w l lt 2 lw 2 color Red title "Analytical derivative function" \\' >> $@

	@echo ',"outfile.txt" index 0 using 1:2 w p pt 31 color Black title "Data"\\'>> $@

outfile.txt: main
	./$< 1> $@ 2>log

main: main.o quasiNewton.o ann.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.log
