import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys

parametrso = sys.argv
I = Image.open(parametrso[1])
I = np.array(I)
I = I[60:80,80:130]

IF = I*0 + 1j*0

# Transformada de Fourier
i0 = 0;
while i0<np.size(I,0):
    j0 = 0;
    while j0<np.size(I,1):
        i1 = 0;
        while i1<np.size(I,0):
            j1 = 0;
            while j1<np.size(I,1):
                IF[i0,j0] += I[i1,j1]*np.cos(-2*3.1415*(i0*i1/np.size(I,0) + j0*j1/np.size(I,1)))
                IF[i0,j0] += 1j*I[i1,j1]*np.sin(-2*3.1415*(i0*i1/np.size(I,0) + j0*j1/np.size(I,1)))
                j1 += 1
            i1 += 1
        j0 += 1
    i0 += 1


Filtro_bajas = IF*0 + 0.0
Filtro_altas = IF*0 + 0.0

w1 = 5
Dw = 5
w2 = w1+Dw

# LLenar las matrices de los filtros
for i in range(np.size(I,0)):
    for j in range(np.size(I,1)):
        fij = 0
        if i < np.size(I,0)/2:
            if j < np.size(I,1)/2:
                fij = np.sqrt(i**2+j**2)
            else:
                fij = np.sqrt(i**2+(np.size(I,1)-j)**2)
        else:
            if j < np.size(I,1)/2:
                fij = np.sqrt((np.size(I,0)-i)**2+j**2)
            else:
                fij = np.sqrt((np.size(I,0)-i)**2+(np.size(I,1)-j)**2)
        
        if fij > w1 and fij<w2: # decaimiento con ecuacion de segundo grado
            # En la banda de transicion
            Filtro_bajas[i,j] = ((fij-w2)/Dw)**2
            Filtro_altas[i,j] = 1.0 -((fij-w2)/Dw)**2
        if fij < w1:
            Filtro_bajas[i,j] = 1.0
        if fij > w2:
            Filtro_altas[i,j] = 1.0

espectro_bajas = I*0 + 0.0j
espectro_altas = I*0 + 0.0j            
for i in range(np.size(I,0)):
    for j in range(np.size(I,1)):
        espectro_bajas[i,j] = IF[i,j]*Filtro_bajas[i,j]
        espectro_altas[i,j] = IF[i,j]*Filtro_altas[i,j]


imagen_bajas = IF*0.0
imagen_altas = IF*0.0

if parametrso[2] == "alto":
	# Transformada de Fourier inversa
	i0 = 0;
	while i0<np.size(I,0):
	    j0 = 0;
	    while j0<np.size(I,1):
	        i1 = 0;
	        while i1<np.size(I,0):
	            j1 = 0;
	            while j1<np.size(I,1):            	            	
	                imagen_altas[i0,j0] += espectro_altas[i1,j1]*np.cos(2*3.1415*(i0*i1/np.size(I,0) + j0*j1/np.size(I,1)))
	                imagen_altas[i0,j0] += 1j*espectro_altas[i1,j1]*np.sin(2*3.1415*(i0*i1/np.size(I,0) + j0*j1/np.size(I,1)))
	                j1 += 1
	            i1 += 1
	        j0 += 1
	    i0 += 1
if parametrso[2] == "bajo":
	# Transformada de Fourier inversa
	i0 = 0;
	while i0<np.size(I,0):
	    j0 = 0;
	    while j0<np.size(I,1):
	        i1 = 0;
	        while i1<np.size(I,0):
	            j1 = 0;
	            while j1<np.size(I,1):            	            	
	                imagen_bajas[i0,j0] += espectro_bajas[i1,j1]*np.cos(2*3.1415*(i0*i1/np.size(I,0) + j0*j1/np.size(I,1)))
	                imagen_bajas[i0,j0] += 1j*espectro_bajas[i1,j1]*np.sin(2*3.1415*(i0*i1/np.size(I,0) + j0*j1/np.size(I,1)))
	                j1 += 1
	            i1 += 1
	        j0 += 1
	    i0 += 1

if parametrso[2] == "alto":
	plt.imshow(imagen_altas.real, cmap='gray'); plt.savefig('altas.png')
if parametrso[2] == "bajo":
	plt.imshow(imagen_bajas.real, cmap='gray'); plt.savefig('bajas.png')