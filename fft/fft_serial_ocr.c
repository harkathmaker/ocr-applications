// Serial version of Cooley-Tukey FFT in OCR framework.
// All logic is performed in a single EDT with no parallelism.
//

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "ditfft2.h"

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
		ditfft2(X_real, X_imag, x_in, N, 1);
	}

	//for(i=0;i<N;i++) {
	//	PRINTF("%d [%f + %fi]\n",i,X_real[i],X_imag[i]);
	//}
	

	PRINTF("Finished successfully\n");
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
	PRINTF("N: %u\n",N);
	for(i=0;i<paramc;i++) {
		PRINTF("%ul\n",paramv[i]);
	}

	u64 params[2] = { N, iterations };

	ocrGuid_t edtGuid;
	ocrEdtCreate(&edtGuid, tempGuid, EDT_PARAM_DEF, params, EDT_PARAM_DEF, NULL_GUID, EDT_PROP_NONE, NULL_GUID, NULL);
	//END-TODO
	return NULL_GUID;
}
