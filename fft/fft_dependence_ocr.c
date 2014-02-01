// Serial version of Cooley-Tukey FFT in plain C.
//

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#define PI 3.1415926

void ditfft2(float *X_real, float *X_imag, float *x_in, int N, int step) {
	if(N == 1) {
		X_real[0] = x_in[0];
		X_imag[0] = 0;
	} else {
		// DFT even side
		ditfft2(X_real, X_imag, x_in, N/2, 2 * step);
		ditfft2(X_real+N/2, X_imag+N/2, x_in+step, N/2, 2 * step);
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

// ex1: Create an EDT

ocrGuid_t fftStartEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	// Use macros from compat.h for fsim/ocr compatibility
	//PRINTF("Hello from EDT\n");
	//PRINTF("paramc: %d\n",paramc);
	int i;
	//for(i=0;i<paramc;i++) {
	//	PRINTF("paramv[%d]: %lu\n",i,paramv[i]);
	//}
	ocrGuid_t startGuid = paramv[0];
	ocrGuid_t endGuid = paramv[1];
	float *data = depv[0].ptr;
	ocrGuid_t dataGuid = depv[0].guid;
	u64 N = paramv[2];
	u64 step = paramv[3];
	u64 offset = paramv[4];
	u64 x_in_offset = paramv[5];
	float *x_in = (float*)data;
	float *X_real = (float*)(data+offset + N*step);
	float *X_imag = (float*)(data+offset + 2*N*step);
	//PRINTF("Step %d offset: %d N*step: %d\n", step, offset, N*step);
	
	//if(step == 1) {
	//for(i=0;i<N;i++) {
	//	PRINTF("%d %f \n",i,x_in[i]);
	//}
	//}

	//ditfft2(X_real, X_imag, x_in, N, 1);
	if(N == 1) {
		X_real[0] = x_in[x_in_offset];
		X_imag[0] = 0;
	} else {
		// DFT even side
		u64 childParamv[6] = { startGuid, endGuid, N/2, 2 * step, 0 + offset, x_in_offset };
		u64 childParamv2[6] = { startGuid, endGuid, N/2, 2 * step, N/2 + offset, x_in_offset + step };
		
		//PRINTF("Creating children of size %d\n",N/2);
		ocrGuid_t edtGuid, edtGuid2, endEdtGuid, finishEventGuid, finishEventGuid2;

		ocrEdtCreate(&edtGuid, startGuid, EDT_PARAM_DEF, childParamv, EDT_PARAM_DEF, &dataGuid, EDT_PROP_FINISH, NULL_GUID, &finishEventGuid);
		ocrEdtCreate(&edtGuid2, startGuid, EDT_PARAM_DEF, childParamv2, EDT_PARAM_DEF, &dataGuid, EDT_PROP_FINISH, NULL_GUID, &finishEventGuid2);
		//PRINTF("finishEventGuid after create: %lu\n",finishEventGuid);
		//ditfft2(X_real, X_imag, x_in, N/2, 2 * step);
		//ditfft2(X_real+N/2, X_imag+N/2, x_in+step, N/2, 2 * step);

		ocrGuid_t endDependencies[3] = { dataGuid, finishEventGuid, finishEventGuid2 };
		// Do calculations after having divided and conquered
		ocrEdtCreate(&endEdtGuid, endGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, endDependencies, EDT_PROP_FINISH, NULL_GUID, NULL);
	}

	//PRINTF("Task with size %d completed\n",N);
	return NULL_GUID;
}

ocrGuid_t fftEndEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	int i;
	ocrGuid_t startGuid = paramv[0];
	ocrGuid_t endGuid = paramv[1];
	float *data = depv[0].ptr;
	ocrGuid_t dataGuid = depv[0].guid;
	u64 N = paramv[2];
	u64 step = paramv[3];
	u64 offset = paramv[4];
	float *x_in = (float*)data+offset;
	float *X_real = (float*)(data+offset + N*step);
	float *X_imag = (float*)(data+offset + 2*N*step);
	//PRINTF("Reached end phase for step %d\n",step);
	//PRINTF("paramc: %d\n",paramc);

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
	//sleep(1);
	return NULL_GUID;
}

ocrGuid_t finalPrintEdt(u32 paramc, u64 *paramv, u32 depc, ocrEdtDep_t depv[]) {
	PRINTF("Final print EDT\n");
	int i;
	//for(i=0;i<paramc;i++) {
	//	PRINTF("paramv[%d]: %lu\n",i,paramv[i]);
	//}
	ocrGuid_t startGuid = paramv[0];
	ocrGuid_t endGuid = paramv[1];
	float *data = depv[1].ptr;
	ocrGuid_t dataGuid = depv[1].guid;
	u64 N = paramv[2];
	u64 step = paramv[3];
	u64 offset = paramv[4];
	float *x_in = (float*)data+offset;
	float *X_real = (float*)(data+offset + N*step);
	float *X_imag = (float*)(data+offset + 2*N*step);

	//for(i=0;i<N;i++) {
	//	PRINTF("%d [%f , %fi]\n",i,X_real[i],X_imag[i]);
	//}

	ocrShutdown();
}

ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	//TODO Create an EDT template for 'appEdt', no argument, no dependence
	ocrGuid_t startTempGuid,endTempGuid,printTempGuid;
	ocrEdtTemplateCreate(&startTempGuid, &fftStartEdt, 6, 1);
	ocrEdtTemplateCreate(&endTempGuid, &fftEndEdt, 6, 3);
	ocrEdtTemplateCreate(&printTempGuid, &finalPrintEdt, 6, 2);
	
	//END-TODO

	int N = 524288;
	// x_in, X_real, and X_imag in a contiguous block
	float *x;
	ocrGuid_t dataGuid;
	DBCREATE(&dataGuid, (void **) &x, sizeof(float) * N * 3, 0, NULL_GUID, NO_ALLOC);
	
	int i;
	for(i=0;i<N;i++) {
		x[i] = 0;
	}
	x[1] = 1;
	x[3] = -3;
	x[4] = 8;
	x[5] = 9;
	x[6] = 1;


	u64 edtParamv[6] = { startTempGuid, endTempGuid, N, 1 /* step size */, 0 /* offset */, 0 /* x_in_offset */ };
	
	PRINTF("sizeof(ocrGuid_t): %d\n",sizeof(ocrGuid_t));
	PRINTF("paramc: %u\n",paramc);
	for(i=0;i<paramc;i++) {
		PRINTF("%lu\n",paramv[i]);
	}

	//sleep(1);

	// Create an EDT out of the EDT template
	ocrGuid_t edtGuid, printEdtGuid, edtEventGuid;
	ocrEdtCreate(&edtGuid, startTempGuid, EDT_PARAM_DEF, edtParamv, EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &edtEventGuid);

	ocrGuid_t finishDependencies[2] = { edtEventGuid, dataGuid };
	ocrEdtCreate(&printEdtGuid, printTempGuid, EDT_PARAM_DEF, edtParamv, EDT_PARAM_DEF, finishDependencies, EDT_PROP_NONE, NULL_GUID, NULL);
	
	ocrAddDependence(dataGuid, edtGuid, 0, DB_MODE_ITW);

	return NULL_GUID;
}
