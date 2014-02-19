// Serial version of Cooley-Tukey FFT in OCR framework.
// All logic is performed in a single EDT with no parallelism.
//

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void ditfft2(float *X_real, float *X_imag, float *x_in, int N, int step, int x_in_offset) {
	if(N == 1) {
		//printf("x_in_offset: %d\n",x_in_offset);
		X_real[0] = x_in[x_in_offset];
		X_imag[0] = 0;
	} else {
		// DFT even side
		ditfft2(X_real, X_imag, x_in, N/2, 2 * step, x_in_offset);
		ditfft2(X_real+N/2, X_imag+N/2, x_in, N/2, 2 * step, x_in_offset + step);
		int k;
		for(k=0;k<N/2;k++) {
			float t_real = X_real[k];
			float t_imag = X_imag[k];
			double twiddle_real;
			double twiddle_imag;
			twiddle_imag = sin(-2 * M_PI * k / N);
			twiddle_real = cos(-2 * M_PI * k / N);
			float xr = X_real[k+N/2];
			float xi = X_imag[k+N/2];

			// (a+bi)(c+di) = (ac - bd) + (bc + ad)i
			X_real[k] = t_real +
				(twiddle_real*xr - twiddle_imag*xi);
			X_imag[k] = t_imag +
				(twiddle_imag*xr + twiddle_real*xi);
			X_real[k+N/2] = t_real -
				(twiddle_real*xr - twiddle_imag*xi);
			X_imag[k+N/2] = t_imag -
				(twiddle_imag*xr + twiddle_real*xi);
		}
	}
}


#define __OCR__
#include "ocr.h"

ocrGuid_t fftEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	// Use macros from compat.h for fsim/ocr compatibility
	int N = paramv[0];
	int iterations = paramv[1];
	float *x_in = (float*)malloc(sizeof(float) * N);
	float *X_real = (float*)malloc(sizeof(float) * N);
	float *X_imag = (float*)malloc(sizeof(float) * N);
	int i,j;
	
	for(i=0;i<N;i++) {
		x_in[i] = 0;
	}
	x_in[1] = 1;
	//x_in[3] = 2;
	//x_in[6] = 5;

	for(i=0;i<iterations;i++) {
		for(j=0;j<N;j++) {
			X_real[i] = 0;
			X_imag[i] = 0;
		}
		ditfft2(X_real, X_imag, x_in, N, 1, 0);
	}

	//for(i=0;i<N;i++) {
	//	PRINTF("%d [%f + %fi]\n",i,X_real[i],X_imag[i]);
	//}
	

	ocrShutdown(); // This is the last EDT to execute
	return NULL_GUID;
}

ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 argc = getArgc(depv[0].ptr);
	int i;

	if(argc < 2) {
		PRINTF("Pass size N to use, where size of matrix = 2^N.\n");
		ocrShutdown();
		return NULL_GUID;
	}
	for(i=0;i<argc;i++) {
		char *argv = getArgv(depv[0].ptr,i);
		PRINTF("argv[%d]: %s\n",i,argv);
	}

	u64 N = pow(2, strtol(getArgv(depv[0].ptr,1),NULL,10) );
	u64 iterations = (argc > 2 ? strtol(getArgv(depv[0].ptr,2),NULL,10) : 1);

	ocrGuid_t tempGuid;
	ocrEdtTemplateCreate(&tempGuid, &fftEdt, 2, 0);
	
	PRINTF("Running %d iterations\n",iterations);
	PRINTF("sizeof(ocrGuid_t): %d\n",sizeof(ocrGuid_t));
	PRINTF("paramc: %u\n",paramc);
	for(i=0;i<paramc;i++) {
		PRINTF("%ul\n",paramv[i]);
	}

	u64 params[2] = { N, iterations };

	ocrGuid_t edtGuid;
	ocrEdtCreate(&edtGuid, tempGuid, EDT_PARAM_DEF, params, EDT_PARAM_DEF, NULL_GUID, EDT_PROP_NONE, NULL_GUID, NULL);
	//END-TODO
	return NULL_GUID;
}
