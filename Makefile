keys = -o3 -fmax-errors=10   -Waliasing -Warray-temporaries -Wc-binding-type -Wconversion -Wconversion-extra -Wsurprising -Wunderflow  -Wrealloc-lhs-all  -Wunused-parameter
compiler = sfort -O3 -xautopar#gfortran -O3 #-fbounds-check -Wsurprising -Wunderflow -Wconversion-extra#sfort -O3 -xautopar #-xloopinfo##-xopenmp #  #-parallel #-Werror# gfortran -Wunused -Wsurprising -Wrealloc-lhs#$(keys)#
#-fopenmp # -o3#
files = fundamental_data.f90 global_module.f90 integrating_module.f90 funct_module.f90 main.f90 # 
#objects_buf = $(patsubst %.f90, %.o, $(files))
objects = $(patsubst %.f90, %.o, $(files)) #$(objects_buf)

main: $(objects)
	$(compiler) -o run $(objects) 

#%.o:  %.f 
#	$(compiler) -c $< -o $@ 

%.o:  %.f90 
	$(compiler) -c $< -o $@ 
	
clear:
	rm *.o
	#rm fundamental_data.o fundamental_data.mod main.o
	rm ./run
	#rm ./graph
	#rm ./slices/*
	rm *.mod
	#rm ./pict/buf/*.gnu
	#rm ./pict/*.eps
	rm *.png

plot: 
#	pdflatex out/graph.tex 		# plot all graphs
	gnuplot plot.gnu 		# plot data for H-function
	gnuplot plot_cont.gnu		# plot reflected full continuum
	gnuplot plot_refl_line.gnu	# plot reflected full line
	
	


