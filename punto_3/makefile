#all: transformada.txt puntos_python

all: transformada.txt

transformada.txt: fourier datos.txt
	./fourier datos.txt


fourier: f.cpp
	#gcc 3_fourier.c -lm -o fourier
	c++ f.cpp -o fourier

puntos_python:
	python 2_filtro.py monalisa.png alto
	python 2_filtro.py monalisa.png bajo
	python 1_suave.py monalisa.png 10