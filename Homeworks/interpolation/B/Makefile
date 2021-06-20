CFLAGS = -Wall -O
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

.PHONEY: all
all: out.qdata.txt out.qxy.txt quadratic_plot.pdf

quadratic_plot.pdf: out.qdata.txt out.qxy.txt Makefile
		echo '\
	set term pdfcairo noenhanced ;\
	set title "Illustration of gsl_interp_eval and gsl_interp_eval_integ" ;\
	set key left ;\
	set xlabel "x" ;\
	set ylabel "y" ;\
	set out "$@" ;\
	plot \
 "$<" index 0 with points pointtype 4 linecolor rgb "black" title "interp value" \
,"$<" index 1 with points linecolor rgb "black" title "integral value +1" \
,"$<" index 2 using 1:2 with lines linecolor rgb "red" title "diff value" \
,"$<" index 3 using 1:2 with lines linecolor rgb "green" title "sin(x) values" \
,"$<" index 4 using 1:2 with lines linecolor rgb "blue" title "cos(x) values" \
	' | gnuplot

out.qdata.txt out.qxy.txt: Quadraturemain
	./$<

PHONEY: clean
clean:
	$(RM) main *.o out* Quadraturemain *.pdf
