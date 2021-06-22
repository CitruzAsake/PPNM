CC = clang
CFLAGS = -Wall -O1 -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

all: cosine.png sir1.png sir2.png sir3.png sir4.png
	@echo '--------------TASKS A and B----------------- '
	@echo 'Check .png plots!'
	@echo 'WHAT we can see is that the effect of Tc on the number of infectious people is noticeable (when the other parameters are fixed)!'
	@echo

cosine.png: cosine.ppl
	pyxplot $< 1> /dev/null 

cosine.ppl: cosine.txt Makefile
	@echo 'set terminal png dpi 600'						>  $@ 
	@echo 'set output "cosine.png"'						>> $@           		
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'											>> $@           		
	@echo 'set yrange [-1.1:1.1]'							>> $@           		
	@echo 'set xlabel "$$x$$"'								>> $@      		
	@echo 'set ylabel "$$y$$"'								>> $@     		
	@echo 'set title "Integration of $$u^{\prime\prime} = u$$ with $$u(0) = 1, \; u^{\prime}(0) = 0$$ "'	>> $@	
	@echo 'plot \\'																								>> $@
	@echo ' cos(x)  					 w l lw 3 color black title "$$\mathrm{cos}(x)$$" 		\\'>> $@
	@echo ',"cosine.txt" using 1:2 w l lw 1 color cyan   title "$$y_0 = u$$" 				\\'>> $@
	@echo ',-sin(x) 					 w l lw 3 color black title "$$-\mathrm{sin}(x)$$" 	\\'>> $@
	@echo ',"cosine.txt" using 1:3 w l lw 1 color red  title "$$y_1 = u^{\prime}$$" 	\\'>> $@

sir1.png: sir1.ppl
	pyxplot $< 1> /dev/null 

sir1.ppl: sir1.txt Makefile
	@echo 'set terminal png dpi 600'						>  $@ 
	@echo 'set output "sir1.png"'							>> $@           		
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'											>> $@           		
	@echo 'set xlabel "$$t/\mathrm{days}$$"'			>> $@      		
	@echo 'set ylabel "$$N$$"'								>> $@     		
	@echo 'set title "SIR model with $$N = 5.8 \\times 10^6, \; T_r = 7, \; T_c = 1.0, \; I_0 = 50$$ "'	>> $@	
	@echo 'plot \\'																								>> $@
	@echo ' "sir1.txt" using 1:2 w l lw 1 lt 1 color black title "Susceptible $$(S)$$" 	\\'>> $@
	@echo ',"sir1.txt" using 1:3 w l lw 1 lt 5 color black title "Infectious $$(I)$$" 	\\'>> $@
	@echo ',"sir1.txt" using 1:4 w l lw 1 lt 3 color black title "Removed $$(R)$$" 	 	\\'>> $@


sir2.png: sir2.ppl
	pyxplot $< 1> /dev/null 

sir2.ppl: sir2.txt Makefile
	@echo 'set terminal png dpi 600'						>  $@ 
	@echo 'set output "sir2.png"'							>> $@           		
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'											>> $@           		
	@echo 'set xlabel "$$t/\mathrm{days}$$"'			>> $@      		
	@echo 'set ylabel "$$N$$"'								>> $@     		
	@echo 'set title "SIR model with $$N = 5.8 \\times 10^6, \; T_r = 7, \; T_c = 1.5, \; I_0 = 50$$ "'	>> $@	
	@echo 'plot \\'																								>> $@
	@echo ' "sir2.txt" using 1:2 w l lw 1 lt 1 color black title "Susceptible $$(S)$$" 	\\'>> $@
	@echo ',"sir2.txt" using 1:3 w l lw 1 lt 5 color black title "Infectious $$(I)$$" 	\\'>> $@
	@echo ',"sir2.txt" using 1:4 w l lw 1 lt 3 color black title "Removed $$(R)$$" 	 	\\'>> $@


sir3.png: sir3.ppl
	pyxplot $< 1> /dev/null 

sir3.ppl: sir3.txt Makefile
	@echo 'set terminal png dpi 600'						>  $@ 
	@echo 'set output "sir3.png"'							>> $@           		
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'											>> $@           		
	@echo 'set xlabel "$$t/\mathrm{days}$$"'			>> $@      		
	@echo 'set ylabel "$$N$$"'								>> $@     		
	@echo 'set title "SIR model with $$N = 5.8 \\times 10^6, \; T_r = 7, \; T_c = 2.5, \; I_0 = 50$$ "'	>> $@	
	@echo 'plot \\'																								>> $@
	@echo ' "sir3.txt" using 1:2 w l lw 1 lt 1 color black title "Susceptible $$(S)$$" 	\\'>> $@
	@echo ',"sir3.txt" using 1:3 w l lw 1 lt 5 color black title "Infectious $$(I)$$" 	\\'>> $@
	@echo ',"sir3.txt" using 1:4 w l lw 1 lt 3 color black title "Removed $$(R)$$" 	 	\\'>> $@


sir4.png: sir4.ppl
	pyxplot $< 1> /dev/null 

sir4.ppl: sir4.txt Makefile
	@echo 'set terminal png dpi 600'						>  $@ 
	@echo 'set output "sir4.png"'							>> $@           		
	@echo 'set key outside' 								>> $@           		
	@echo 'set grid'											>> $@           		
	@echo 'set xlabel "$$t/\mathrm{days}$$"'			>> $@      		
	@echo 'set ylabel "$$N$$"'								>> $@     		
	@echo 'set title "SIR model with $$N = 5.8 \\times 10^6, \; T_r = 7, \; T_c = 4.0, \; I_0 = 50$$ "'	>> $@	
	@echo 'plot \\'																								>> $@
	@echo ' "sir4.txt" using 1:2 w l lw 1 lt 1 color black title "Susceptible $$(S)$$" 	\\'>> $@
	@echo ',"sir4.txt" using 1:3 w l lw 1 lt 5 color black title "Infectious $$(I)$$" 	\\'>> $@
	@echo ',"sir4.txt" using 1:4 w l lw 1 lt 3 color black title "Removed $$(R)$$" 	 	\\'>> $@


%.txt: main
	./$<

main: main.o RK.o ../utilities.o

clean:
	$(RM) main *.o *.txt *.ppl *.pdf *.png ../utilities.o

test:
	@echo $(CFLAGS)
	@echo $(LDLIBS)

