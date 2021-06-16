CC = clang
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: plot.png logplot.png output.txt
	@echo
	@cat output.txt
	@echo 'Compare to modern half-life for Ra-224 (from Wikipedia): (3.6319 +/- 0.0023) days = [3.6296, 3.6342] days'
	@echo 'My result does not agree within one standard deviation but it is not too far off!'
	@echo 'Check .png plots!'
	@echo

plot.png: plot.ppl
	pyxplot $< 1> /dev/null 

plot.ppl: data.txt Makefile
	@echo 'set terminal png dpi 600'						>  $@ 
	@echo 'set output "plot.png"'							>> $@           		
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'											>> $@           		
	@echo 'set xlabel "$$y$$"'								>> $@      		
	@echo 'set ylabel "$$t/\mathrm{days}$$"'			>> $@     		
	@echo 'set title "Raw data with error bars"'		>> $@	
	@echo 'plot \\'											>> $@
	@echo ' "data.txt" using 1:2:3 with yerrorbars color black title "Measurements" 	\\'>> $@

logplot.png: logplot.ppl
	pyxplot $< 1> /dev/null

logplot.ppl: logdata.txt fit.txt params.txt Makefile
	@echo 'set terminal png dpi 600'						>  $@ 
	@echo 'set output "logplot.png"'						>> $@           		
	@echo 'params = open("params.txt", "r")'        >> $@
	@echo 'A = params.readline()' 						>> $@
	@echo 'B = params.readline()'                   >> $@
	@echo 'params.close()'                          >> $@
	@echo 'lna = A.split()[0]'                      >> $@
	@echo 'dlna = A.split()[1]'                     >> $@
	@echo 'lambda  = B.split()[0]'						>> $@
	@echo 'dlambda = B.split()[1]'                  >> $@
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'											>> $@           		
	@echo 'set xlabel "$$\ln(y)$$"'						>> $@      		
	@echo 'set ylabel "$$\ln(t/\mathrm{days})$$"'	>> $@     		
	@echo 'set title "Ordinary least squares"'		>> $@	
	@echo 'plot \\'											>> $@
	@echo ' "fit.txt" using 1:2 w l lw 1 lt 5 color black title "$$\ln(y) = (%.2f \pm %.2f) + (%.3f \pm %.3f)t$$"%(lna, dlna, lambda, dlambda) \\'>> $@
	@echo ',"fit.txt" using 1:3 w l lw 1 lt 3 color black title "" 					\\'>> $@
	@echo ',"fit.txt" using 1:4 w l lw 1 lt 3 color black title "" 					\\'>> $@
	@echo ',"logdata.txt" using 1:2:3 with yerrorbars color blue title "Measurements" 	\\'>> $@

logdata.txt fit.txt params.txt output.txt: main time.txt activity.txt unceractivity.txt
	./$< > output.txt

data.txt:
	@echo  1  117     5.85    > $@
	@echo  2  100     5.0     >> $@
	@echo  3   88     4.4     >> $@
	@echo  4   72     3.6     >> $@
	@echo  6   53     2.65    >> $@
	@echo  9   29.5   1.475   >> $@
	@echo 10   25.2   1.26    >> $@
	@echo 13   15.2   0.76    >> $@
	@echo 15   11.1   0.555   >> $@ 

time.txt act.txt dact.txt: data.txt
	cut -d " " -f 1 data.txt > time.txt
	cut -d " " -f 2 data.txt > activity.txt
	cut -d " " -f 3 data.txt > unceractivity.txt

main: main.o gs.o leastsquares.o ../util.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png ../util.o

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)

