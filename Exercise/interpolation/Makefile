CFLAGS = -Wall
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all:Out.plot.pdf

Out.plot.pdf: Out.data.txt Makefile
	echo '\
	set term pdfcairo noenhanced ;\
	set title "Illustration of gsl_interp_eval and gsl_interp_eval_integ" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
 "$<" index 0 with points pointtype 4 linecolor rgb "black" title "data" \
 ,"$<" index 1 using 1:2 with lines linecolor rgb "red" title "linear-spline" \
	' | gnuplot

Out.data.txt: main
	./$< > $@

clean:
	$(RM) Out* main
