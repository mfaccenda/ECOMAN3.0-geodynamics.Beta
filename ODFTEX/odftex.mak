odftex.x 	: odftex.o          
	  gfortran -O odftex.o -o odftex.x 
odftex.o	: odftex.f90       
	  gfortran -O  -c odftex.f90
