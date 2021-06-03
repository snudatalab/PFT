/*****************************************************************************
*    This is an implementation of Fast and Accurate                          *
*    Partial Fourier Transform for Time Series Data (submitted to KDD21).    *
*                                                                            *
*    Note that this code example is specifically designed for                *
*    real-valued input with a target range centered at zero                  *
*    for the best performance. A slight modification removes the             *
*    constraint; please refer to the original paper for more details.        * 
*                                                                            *
*    This code also contains the implementation of MKL DFTI.                 *
*                                                                            *
*    https://github.com/FPFTanonymous/FPFT                                   *
*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"
#include "ipps.h"
#include "ipp.h"
#include "mkl_dfti.h"
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

int main()
{
	// 2^10=1024 // 2^11=2048 // 2^12=4096 // 2^13=8192 // 2^14=16384 
	// 2^15=32768 // 2^16=65536 // 2^17=131072 // 2^18=262144 // 2^19=524288
	// 2^20=1048576 // 2^21=2097152 // 2^22=4194304 // 2^23=8388608 // 2^24=16777216
	int N = 4194304;
	int M = 131072;
	int p = 65536;
	string error = "e-7";
	// One of e-1, e-2, e-3, e-4, e-5, e-6, e-7 is allowed

	int num_D = ceil(M * (1.0 / p));
	int q = N / p;
	int r = 0;

	/* Pointers */
	double* W, * XI; 
	W = (double*)mkl_malloc(27 * sizeof(double), 64);
	XI = (double*)mkl_malloc(26 * sizeof(double), 64);

	/* Load precomputed xi */
	string file_name = "precomputed/" + error + ".csv";
	ifstream ip(file_name);
	string sd;
	int count = 0;

	while (count < 26)
	{
		getline(ip, sd, ',');
		XI[count] = stod(sd);
		count += 1;
	}

	/* Find r */
	while (XI[r] < ((1.0 * M) / p))
	{
		r += 1;
		if (r == 25)
			break;
	}
	r += 2;

	std::cout << "N, M = " << N << ", " << M << " // ";
	std::cout << "p, q, r = " << p << ", " << q << ", " << r << " // ";
	std::cout << "e = " << error << "\n";

	/* Load precomputed w */
	int w_start = (r * (r - 1)) / 2 - 1;

	for (int op = 0; op < r - 2; op++)
	{
		getline(ip, sd);
	}
	for (int op = 0; op < r; op++)
	{
		getline(ip, sd, ',');
		W[op] = stod(sd);
	}


	int ppt = p + 2;
	int ppth = ppt / 2;
	int pp = ppt * 2;
	int rr = r / 2;

	/* Pointers */
	float* A, * B, * C, * D, * Z, * F, * DS, * FS, * FSC, * FSS, * TCC, * TSS, * OUT;

	/* Initialize array */
	A = (float*)mkl_malloc(N * sizeof(float), 32);
	B = (float*)mkl_malloc(q * r * sizeof(float), 32);
	C = (float*)mkl_malloc(ppt * r * sizeof(float), 32);
	D = (float*)mkl_malloc(ppt * r * sizeof(float), 32);
	Z = (float*)mkl_malloc((N + 2) * sizeof(float), 32);
	F = (float*)mkl_malloc(pp * sizeof(float), 32);
	DS = (float*)mkl_malloc(ppt * r * num_D * sizeof(float), 32);
	FS = (float*)mkl_malloc(pp * num_D * sizeof(float), 32);
	FSC = (float*)mkl_malloc(ppt * num_D * sizeof(float), 32);
	FSS = (float*)mkl_malloc(ppt * num_D * sizeof(float), 32);
	TCC = (float*)mkl_malloc(ppt * num_D * sizeof(float), 32);
	TSS = (float*)mkl_malloc(ppt * num_D * sizeof(float), 32);
	OUT = (float*)mkl_malloc((M + 2) * sizeof(float), 32); 
	int n, k, l, j;

	/* Generate A */
	std::srand(std::time(0));
	for (n = 0; n != N; n++)
	{
		A[n] = (float)(((float)std::rand()) / RAND_MAX); // random 0 ~ 1
	}

	printf("\nShow A for %d,...,%d\n", 0, 7);
	for (int op = 0; op < 8; op++)
	{
		printf("%e  ", A[op]);
	}
	std::cout << "\n";

	/* Precompute B */
	for (j = 0; j != r; j++)
		for (l = 0; l != q; l++)
		{
			B[q * j + l] = (float)(W[j] * pow(-2.0 * (l - (q / 2.0)) * (1.0 / q), j));
		}

	/* Precompute D */
	for (j = 0; j != r; j++)
		for (k = 0; k != ppt; k++)
		{
			if (j % 4 == 0 || j % 4 == 3)
				D[ppt * j + k] = (float)(pow((k / 2) * (1.0 / p), j));
			else
				D[ppt * j + k] = (float)(-pow((k / 2) * (1.0 / p), j));
		}

	/* Precompute DS */
	for (int si = 0; si != num_D; si++)
	{
		for (j = 0; j != r; j++)
			for (k = 0; k != ppt; k++)
			{
				if (j % 4 == 0 || j % 4 == 3)
				{
					if (si % 2 == 0)
						DS[ppt * (si * r + j) + k] = (float)(pow(((si * p + k) / 2) * (1.0 / p), j));
					else
						DS[ppt * (si * r + j) + ppt - 1 - k] = (float)(pow(((si * p + k) / 2) * (1.0 / p), j));
				}
				else
				{
					if (si % 2 == 0)
						DS[ppt * (si * r + j) + k] = (float)(-pow(((si * p + k) / 2) * (1.0 / p), j));
					else 
						DS[ppt * (si * r + j) + ppt - 1 - k] = (float)(-pow(((si * p + k) / 2) * (1.0 / p), j));
				}
			}
	}

	/* Precompute twiddle factors */
	for (int si = 0; si != num_D; si++)
	{
		for (k = 0; k != ppt; k++)
		{
			TCC[ppt * si + k] = (float)(cos(((p * si + k) / 2) * M_PI * (1.0 / p)));
			TSS[ppt * si + k] = (float)(-sin(((p * si + k) / 2) * M_PI * (1.0 / p)));
		}

	}




	/* FFT */
	DFTI_DESCRIPTOR_HANDLE hand;
	DftiCreateDescriptor(&hand, DFTI_SINGLE, DFTI_REAL, 1, N);
	DftiSetValue(hand, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
	DftiCommitDescriptor(hand);
	/*********** DO FFT **********/
	DftiComputeForward(hand, A, Z);
	/*****************************/
	DftiFreeDescriptor(&hand);



	/* PFT */
	//DFTI_DESCRIPTOR_HANDLE hand;
	DftiCreateDescriptor(&hand, DFTI_SINGLE, DFTI_REAL, 1, p);
	DftiSetValue(hand, DFTI_NUMBER_OF_TRANSFORMS, r);
	DftiSetValue(hand, DFTI_INPUT_DISTANCE, ppt);
	DftiCommitDescriptor(hand);
	/*********** DO PFT **********/
	/* Matrix Multiplication */
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, r, p, q, 1.0, B, q, A, q, 0.0, C, ppt);
	/* FFTs */
	DftiComputeForward(hand, C);
	/* Post-processing */
	for (int si = 0; si != num_D; si++)
	{
		ippsMul_32f(C, DS + r * ppt * si, FS + pp * si, pp);
		for (int i = 1; i != (r / 2); ++i) {
			ippsAddProduct_32f(C + pp * i, DS + r * ppt * si + pp * i, FS + pp * si, pp);
		}
		if (r % 2 == 1)
			ippsAddProduct_32f(C + pp * rr, DS + r * ppt * si + pp * rr, FS + pp * si, ppt); 
		if (si % 2 == 0)
		{
			cblas_saxpy(ppth, 1.0, FS + pp * si + ppt + 1, 2, FS + pp * si, 2);
			cblas_saxpy(ppth, -1.0, FS + pp * si + ppt, 2, FS + pp * si + 1, 2);
		}
		else
		{
			cblas_saxpy(ppth, -1.0, FS + pp * si + ppt + 1, 2, FS + pp * si, 2);
			cblas_saxpy(ppth, 1.0, FS + pp * si + ppt, 2, FS + pp * si + 1, 2);
		}
		ippsMul_32f(FS + pp * si, TCC + ppt * si, FSC + ppt * si, ppt);
		ippsMul_32f(FS + pp * si, TSS + ppt * si, FSS + ppt * si, ppt);
		if (si % 2 == 0)
		{
			cblas_saxpy(ppth, -1.0, FSS + ppt * si + 1, 2, FSC + ppt * si, 2);
			cblas_saxpy(ppth, 1.0, FSS + ppt * si, 2, FSC + ppt * si + 1, 2);
		}
		else
		{
			cblas_saxpy(ppth, -1.0, FSS + ppt * si + 1, 2, FSC + ppt * si, 2);
			cblas_saxpy(ppth, 1.0, FSS + ppt * si, 2, FSC + ppt * si + 1, 2);
		}
	}
	/*****************************/
	DftiFreeDescriptor(&hand);


	/* Rearrange output */
	for (int si = 0; si != num_D; si++)
	{
		if (si % 2 == 0)
		{
			for (k = 0; k != ppt; k++)
			{
				OUT[p * si + k] = FSC[p * si + 2 * si + k];
			}
		}
		else
		{
			for (k = 0; k != ppt; k++)
			{
				OUT[p * si + k] = FSC[ppt * (si + 1) - k - 1];
			}
		}

	}

	printf("\nShow the first few Fourier coefficients\n");

	for (int op = 0; op < 0 + 6; op++)
	{
		printf("%e  ", Z[op]);
	}
	std::cout << "<-- FFT\n";

	for (int op = 0; op < 0 + 6; op++)
	{
		printf("%e  ", OUT[op]);
	}
	std::cout << "<-- PFT\n";

	printf("\nShow the last few Fourier coefficients\n");

	for (int op = M - 4; op < M + 2; op++)
	{
		printf("%e  ", Z[op]);
	}
	std::cout << "<-- FFT\n";

	for (int op = M - 4; op < M + 2; op++)
	{
		printf("%e  ", OUT[op]);
	}
	std::cout << "<-- PFT\n\n";



	Ipp32f pNor;
	Ipp32f fNor;

	for (int x = 0; x < M; x += 2)
	{
		DS[x / 2] = sqrt(pow(OUT[x] - Z[x], 2) + pow(OUT[x + 1] - Z[x + 1], 2));
	}
	ippsNorm_L2_32f(DS, M / 2, &pNor);

	for (int x = 0; x < M; x += 2)
	{
		DS[x / 2] = sqrt(pow(Z[x], 2) + pow(Z[x + 1], 2));
	}
	ippsNorm_L2_32f(DS, M / 2, &fNor);

	printf("%e  ", pNor);
	std::cout << "<-- diff_norm_L2 \n";
	printf("%e  ", fNor);
	std::cout << "<-- norm_L2 \n";
	printf("%e  ", pNor / fNor);
	std::cout << "<-- diff_norm_L2/norm_L2 \n\n";

	mkl_free(W);
	mkl_free(XI);
	mkl_free(A);
	mkl_free(B);
	mkl_free(C);
	mkl_free(D);
	mkl_free(Z);
	mkl_free(F);
	mkl_free(DS);
	mkl_free(FS);
	mkl_free(FSC);
	mkl_free(FSS);
	mkl_free(TCC);
	mkl_free(TSS);
	mkl_free(OUT);
}


