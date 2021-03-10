CFLAGS = -O -Wall -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: exp1.pdf func.pdf

exp1.pdf: data.txt Makefile
	echo '\
	set terminal pdf;\
	set output "$@";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Exponential function";\
	plot \
	"$<" using 1:2 with line title "expo func"\
	' | gnuplot

data.txt: main
	./$< > $@

.PHONEY: clean
clean: 
	$(RM) main *.png data.txt

func.pdf: func.tex
	pdflatex $<
	pdflatex $<
