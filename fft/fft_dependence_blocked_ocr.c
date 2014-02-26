// OCR implementation of the Cooley-Tukey algorithm. Same as
// fft_dependence_ocr.c, but recursive creation of StartEDTs stops once
// the matrix size reaches SERIAL_BLOCK_SIZE. For these small matrices ditfft2
// is called to compute the answer serially. This is meant to minimize the overhead
// of creating EDTs while still maximizing parallelism.
//
// EndEDTs are also changed to divide their work to a number of slave EDTs, such that
// each slave handles SERIAL_BLOCK_SIZE elements.
//

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SERIAL_BLOCK_SIZE (1024*16)
#define PRINT_RESULTS 0

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
	ocrGuid_t endSlaveGuid = paramv[2];
	float *data = depv[0].ptr;
	ocrGuid_t dataGuid = depv[0].guid;
	u64 N = paramv[3];
	u64 step = paramv[4];
	u64 offset = paramv[5];
	u64 x_in_offset = paramv[6];
	float *x_in = (float*)data;
	float *X_real = (float*)(data+offset + N*step);
	float *X_imag = (float*)(data+offset + 2*N*step);
	//PRINTF("Step %d offset: %d N*step: %d\n", step, offset, N*step);
	
	//if(step == 1) {
	//for(i=0;i<N;i++) {
	//	PRINTF("%d %f \n",i,x_in[i]);
	//}
	//}

	if(N <= SERIAL_BLOCK_SIZE) {
		ditfft2(X_real, X_imag, x_in+x_in_offset, N, step);
	} else {
		// DFT even side
		u64 childParamv[7] = { startGuid, endGuid, endSlaveGuid, N/2, 2 * step, 0 + offset, x_in_offset };
		u64 childParamv2[7] = { startGuid, endGuid, endSlaveGuid, N/2, 2 * step, N/2 + offset, x_in_offset + step };
		
		//PRINTF("Creating children of size %d\n",N/2);
		ocrGuid_t edtGuid, edtGuid2, endEdtGuid, finishEventGuid, finishEventGuid2;

		ocrEdtCreate(&edtGuid, startGuid, EDT_PARAM_DEF, childParamv, EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &finishEventGuid);
		ocrEdtCreate(&edtGuid2, startGuid, EDT_PARAM_DEF, childParamv2, EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &finishEventGuid2);
		//PRINTF("finishEventGuid after create: %lu\n",finishEventGuid);

		ocrGuid_t endDependencies[3] = { dataGuid, finishEventGuid, finishEventGuid2 };
		// Do calculations after having divided and conquered
		ocrEdtCreate(&endEdtGuid, endGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, endDependencies, EDT_PROP_FINISH, NULL_GUID, NULL);

		ocrAddDependence(dataGuid, edtGuid, 0, DB_MODE_ITW);
		ocrAddDependence(dataGuid, edtGuid2, 0, DB_MODE_ITW);
	}

	//PRINTF("Task with size %d completed\n",N);
	return NULL_GUID;
}

ocrGuid_t fftEndEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	int i;
	ocrGuid_t startGuid = paramv[0];
	ocrGuid_t endGuid = paramv[1];
	ocrGuid_t endSlaveGuid = paramv[2];
	float *data = depv[0].ptr;
	ocrGuid_t dataGuid = depv[0].guid;
	u64 N = paramv[3];
	u64 step = paramv[4];
	u64 offset = paramv[5];
	float *x_in = (float*)data+offset;
	float *X_real = (float*)(data+offset + N*step);
	float *X_imag = (float*)(data+offset + 2*N*step);
	//PRINTF("Reached end phase for step %d\n",step);
	//PRINTF("paramc: %d\n",paramc);
	//
	ocrGuid_t *slaveGuids;
	u64 *slaveParamv;

	if(N/2 > SERIAL_BLOCK_SIZE) {
		// TODO: malloc will leak memory
		slaveGuids = malloc(sizeof(ocrGuid_t) * (N/2)/SERIAL_BLOCK_SIZE);
		slaveParamv = malloc(sizeof(u64) * 5 * (N/2)/SERIAL_BLOCK_SIZE);

#if PRINT_RESULTS==1
		PRINTF("Creating %d slaves for N=%d\n",(N/2)/SERIAL_BLOCK_SIZE,N);
#endif
		
		for(i=0;i<(N/2)/SERIAL_BLOCK_SIZE;i++) {
			slaveParamv[i*5] = N;
			slaveParamv[i*5+1] = step;
			slaveParamv[i*5+2] = offset;
			slaveParamv[i*5+3] = i*SERIAL_BLOCK_SIZE;
			slaveParamv[i*5+4] = (i+1)*SERIAL_BLOCK_SIZE;

			ocrEdtCreate(slaveGuids+i, endSlaveGuid, EDT_PARAM_DEF, slaveParamv+i*5, EDT_PARAM_DEF, &dataGuid, EDT_PROP_NONE, NULL_GUID, NULL);
		}
	} else {
		slaveGuids = malloc(sizeof(ocrGuid_t));
		slaveParamv = malloc(sizeof(u64*) * 5);
		
		slaveParamv[0] = N;
		slaveParamv[1] = step;
		slaveParamv[2] = offset;
		slaveParamv[3] = 0;
		slaveParamv[4] = N/2;

		ocrEdtCreate(slaveGuids, endSlaveGuid, EDT_PARAM_DEF, slaveParamv, EDT_PARAM_DEF, &dataGuid, EDT_PROP_NONE, NULL_GUID, NULL);
	}
	return NULL_GUID;
}
ocrGuid_t fftEndSlaveEdt(u32 paramc, u64 *paramv, u32 depc, ocrEdtDep_t depv[]) {
	int i;
	float *data = depv[0].ptr;
	ocrGuid_t dataGuid = depv[0].guid;
	u64 N = paramv[0];
	u64 step = paramv[1];
	u64 offset = paramv[2];
	float *x_in = (float*)data+offset;
	float *X_real = (float*)(data+offset + N*step);
	float *X_imag = (float*)(data+offset + 2*N*step);
	u64 kStart = paramv[3];
	u64 kEnd = paramv[4];

	int k;
	for(k=kStart;k<kEnd;k++) {
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
	ocrGuid_t endSlaveGuid = paramv[2];
	float *data = depv[1].ptr;
	ocrGuid_t dataGuid = depv[1].guid;
	u64 N = paramv[3];
	u64 step = paramv[4];
	u64 offset = paramv[5];
	float *x_in = (float*)data+offset;
	float *X_real = (float*)(data+offset + N*step);
	float *X_imag = (float*)(data+offset + 2*N*step);

#if PRINT_RESULTS==1
	for(i=0;i<N;i++) {
		PRINTF("%d [%f , %fi]\n",i,X_real[i],X_imag[i]);
	}
#endif

	ocrShutdown();
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

	// Create an EDT template for 'appEdt', no argument, no dependence
	ocrGuid_t startTempGuid,endTempGuid,printTempGuid,endSlaveTempGuid;
	ocrEdtTemplateCreate(&startTempGuid, &fftStartEdt, 7, 1);
	ocrEdtTemplateCreate(&endTempGuid, &fftEndEdt, 7, 3);
	ocrEdtTemplateCreate(&endSlaveTempGuid, &fftEndSlaveEdt, 5, 1);
	ocrEdtTemplateCreate(&printTempGuid, &finalPrintEdt, 7, 2);
	
	int N = pow( 2, strtol(getArgv(depv[0].ptr,1), NULL, 10) );
	// x_in, X_real, and X_imag in a contiguous block
	float *x;
	ocrGuid_t dataGuid;
	// TODO: OCR cannot handle large datablocks
	DBCREATE(&dataGuid, (void **) &x, sizeof(float) * N * 3, 0, NULL_GUID, NO_ALLOC);
	printf("Datablock of size %lu (N=%lu) created\n",sizeof(float)*N*3,N);
	
	for(i=0;i<N;i++) {
		x[i] = 0;
	}
	x[1] = 1;
	x[3] = -3;
	x[4] = 8;
	x[5] = 9;
	x[6] = 1;


	u64 edtParamv[7] = { startTempGuid, endTempGuid, endSlaveTempGuid, N, 1 /* step size */, 0 /* offset */, 0 /* x_in_offset */ };
	
	// Create an EDT out of the EDT template
	ocrGuid_t edtGuid, printEdtGuid, edtEventGuid;
	ocrEdtCreate(&edtGuid, startTempGuid, EDT_PARAM_DEF, edtParamv, EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &edtEventGuid);

	ocrGuid_t finishDependencies[2] = { edtEventGuid, dataGuid };
	ocrEdtCreate(&printEdtGuid, printTempGuid, EDT_PARAM_DEF, edtParamv, EDT_PARAM_DEF, finishDependencies, EDT_PROP_NONE, NULL_GUID, NULL);
	
	ocrAddDependence(dataGuid, edtGuid, 0, DB_MODE_ITW);

	return NULL_GUID;
}
