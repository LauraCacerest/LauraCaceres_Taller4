#include "iostream"
#include "fstream"
#include <cstdio>
#include "cstdlib"
#include <math.h>
#include <complex.h>
using namespace std;

//void display(float n[2][16])
//{
//
//    cout << "Displaying Values: " << endl;
//    for(int i = 0;  i < 2; ++i)
//    {
//        for(int j = 0; j < 2; ++j)
//        {
//            cout << n[i][j] << " ";
//        }
//    }
//}

const int n=16;
void interpolar(float OData[2][n], float NData[2][n], float fs)
{
	int n = 16;
	int i;
	float pj;

	i = 0;
	do
	{
		NData[0][i] = i*1/fs;
		i++;
	} while (i<n);	

	
	i = 0;
	do
	{
		NData[1][i] = 0.0;
		for (int i1 = 0; i1 < n; ++i1)
		{
			pj = 1.0;
			for (int i2 = 0; i2 < n; ++i2)
			{
				if (i2 != i1)
				{
					pj = pj*(NData[0][i]-OData[0][i2])/(OData[0][i1]-OData[0][i2]);
				}
				else
				{
					pj = pj*1.0;
				}
				
			}
			NData[1][i] = NData[1][i] + OData[1][i1]*pj;
		}

		i++;
	} while (i<n);	
}

void transformada_fourier(float NData[2][n], float DFT[2][n], float fs)
{	
	int n = 16;
	float complex z;
	int i;
	int j;
	printf("%f\n", fs/2);
	i = 0;
	do
	{
		j = 0;
		do
		{
			z = NData[0][i]*cexp(-6.18*i*j/n*I);
			j++;
		}while (j<n);
		DFT[1][i] = creal(z);
		DFT[2][i] = cimag(z);

		DFT[0][i] = fs/2*((float)i/(float)(n/2));
		if (i >= n/2)
		{
			DFT[0][i] = fs/2*((float)(i-n/2)/(float)(n/2)-1.0);
		}

		i++;
	} while (i<n);
}

void exportar_datos(char *fname, float Data[2][n], int n0, int n1)
{
	FILE *out;

	out = fopen(fname, "w");
	for (int i = 0; i < n0; ++i)
	{
		for (int j = 0; j < n1; ++j)
		{
			fprintf(out, "%f ", Data[j][i]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
}

int main(int argc, char *argv[])
{	
	// Declarar Variables -----------------------------
	int i;
	//const int n=16;
	int sz;
	float fs; // frecuencia de muestreo

	FILE *in;
	

	// Inicializar Variables -----------------------------
	//n = 16;
	in = fopen(argv[1],"r");
	sz = n*sizeof(float);
	
	float OData[2][n]; // Datos originales
	float NData[2][n]; // Nuevos datos
	float DFT[3][n]; // Datos de la transformada
//	OData = malloc(2*sz);
//	NData = malloc(2*sz);
//	for (i = 0; i < 2; ++i)
//	{
//		OData[i]=malloc(sz);
//		NData[i]=malloc(sz);
//	}
//	DFT = malloc(3*sz);
//	for (i = 0; i < 3; ++i)	{DFT[i]=malloc(sz);}


	// Cargar Datos -------------------------------------
	i = 0;
	ifstream inn("datos.txt");

	do
	{

		fscanf(in,"%f %f\n",&OData[0][i],&OData[1][i]);
		//printf("%f %f\n", OData[0][i], OData[1][i]);
		i++;
	} while (i<n);

	// Generar senal muestreada uniformemente e interpolar -------------------------------------
	fs = 1.0/ ( (OData[0][n-1]-OData[0][0])/(n-1) ); // frecuencia de muestreo
	printf("\n\n%f\n\n",fs);

	float num[2][16];
	float copy_OData[2][n];
	int j;
	for (i = 0; i < 2; ++i)
	{
		for (j = 0; j < 2; ++j)
		{
			copy_OData[i][j] = i*j;
		}
	}

	interpolar(OData, NData, fs);

	// Calcular transformada de Fourier
	transformada_fourier(NData, DFT, fs);
	
	int nn = n;
	// Eportar datos
	exportar_datos("transformada.txt", DFT, nn, 3);
	
	return 0;
}