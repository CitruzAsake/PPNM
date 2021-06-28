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

	@echo 'set key top'								>> $@

	@echo 'set xlabel "$$x$$"'							>> $@

	@echo 'set ylabel "$$y$$"'							>> $@

	@echo 'set title "Artificial Neural Network different fits"'			>> $@

	@echo 'plot \\'									>> $@

	@echo '"outfile.txt"  index 1 using 1:1 w l lw 2 color black  title "Integral fit func" \\'>> $@

	@echo ',"outfile.txt" index 1 using 1:2 w l lw 2 color cyan title "Curve fit func" \\'>> $@

	@echo ',"outfile.txt" index 1 using 1:3 w l lw 2 color SpringGreen title "Derivative fit func" \\'>> $@

	@echo ',exp(-x*x)*(-5*sin(5*x-1)-2*x*cos(5*x-1)) w l lt 2 lw 2 color blue title "Analytical derivative function" \\' >> $@

	@echo ',-0.1877184794-0.045046858263+cos(1.0)*x+5.0/2.0*sin(1.0)*x**2-9.0/2.0*cos(1.0)*x**3-155.0/24.0*sin(1.0)*x**4+937.0/120.0*cos(1.0)*x**5 w l lt 5 lw 2 color blue title "Analytical Taylor expansion function" \\' >> $@

	@echo ',"outfile.txt" index 0 using 1:2 w p pt 31 color Red title "Data"\\'>> $@

outfile.txt: main
	./$< 1> $@ 2>log

main: main.o quasiNewton.o ann.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.log
