CC = gcc
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: quadraturespline.pdf
	@echo '30 data points made with sin(i*i)'

output.txt: main
	./$< > $@

quadraturespline.pdf: quadraturespline.ppl
	pyxplot $< 1> /dev/null

quadraturespline.ppl: output.txt Makefile
	@echo 'set terminal pdf dpi 800'						>> $@
	@echo 'set output "quadraturespline.pdf"'		>> $@
	@echo 'set key outside' 								>> $@
	@echo 'set grid'											>> $@
	@echo 'set xlabel "$$x$$"'								>> $@
	@echo 'set ylabel "$$y$$"'								>> $@
	@echo 'set title "Interpolation by Quadraturespline "'	>> $@
	@echo 'plot \\'																								>> $@
	@echo '"output.txt" index 0 using 1:2 w p pt 31 ps 1 color yellow title "Data" 		\\'>> $@
	@echo ',"output.txt" index 1 using 1:2 w l lt 1 lw 1 ps 1 color black  title "Interpolation data" 		\\'>> $@
	@echo ',"output.txt" index 2 using 1:2 w l lt 5 lw 1 ps 0.75 color Cyan  title "Integral data" 		\\'>> $@
	@echo ',"output.txt" index 3 using 1:2 w l lt 6 lw 1 ps 0.75 color Mahogany  title "Derivative" 		\\'>> $@
main: main.o utilities.o quadraturespline.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png
