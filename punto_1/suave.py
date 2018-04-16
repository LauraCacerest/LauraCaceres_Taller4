import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import sys

parametrso = sys.argv
I = Image.open(parametrso[1])
I = np.array(I)
I = I[55:85,80:130]

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


# Crear kernel gaussiano
nn = float(parametrso[2])
n0 = -nn/2
n1 = nn/2


m = np.size(I,0)
n = np.size(I,1)
K = np.zeros((m,n))

for i in range(int(n1-n0)):
    for j in range(int(n1-n0)):
        i1 = i + n0
        j1 = j + n0
        K[i,j] = np.exp(-1.0*(1.0/(nn/2)*(i1**2+j1**2)))

# Transformada del kernel
KF = K*0 + 1j*0
i0 = 0;
while i0<np.size(K,0):
    j0 = 0;
    while j0<np.size(K,1):
        i1 = 0;
        while i1<np.size(K,0):
            j1 = 0;
            while j1<np.size(K,1):
                KF[i0,j0] += K[i1,j1]*np.cos(-2*3.1415*(i0*i1/np.size(K,0) + j0*j1/np.size(K,1)))
                KF[i0,j0] += 1j*K[i1,j1]*np.sin(-2*3.1415*(i0*i1/np.size(K,0) + j0*j1/np.size(K,1)))
                j1 += 1
            i1 += 1
        j0 += 1
    i0 += 1

# APlicar kernel
imagen_filtrada_espectro = KF*0 + 0.0j
for i in range(np.size(KF,0)):
    for j in range(np.size(KF,1)):
        imagen_filtrada_espectro[i,j] = IF[i,j]*KF[i,j]
        
# transformada inversa
imagen_filtrada = IF*0.0
i0 = 0;
while i0<np.size(I,0):
    j0 = 0;
    while j0<np.size(I,1):
        i1 = 0;
        while i1<np.size(I,0):
            j1 = 0;
            while j1<np.size(I,1):
                imagen_filtrada[i0,j0] += imagen_filtrada_espectro[i1,j1]*np.cos(2*3.1415*(i0*i1/np.size(I,0) + j0*j1/np.size(I,1)))
                imagen_filtrada[i0,j0] += 1j*imagen_filtrada_espectro[i1,j1]*np.sin(2*3.1415*(i0*i1/np.size(I,0) + j0*j1/np.size(I,1)))
                j1 += 1
            i1 += 1
        j0 += 1
    i0 += 1


plt.figure()
plt.imshow(imagen_filtrada.real, cmap='gray')
plt.savefig('suave.png')