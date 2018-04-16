all: transformada.txt

transformada.txt: fourier datos.txt
	./fourier datos.txt


fourier: f.cpp
	c++ f.cpp -o fourier
