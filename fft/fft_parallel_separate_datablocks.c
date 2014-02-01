// Serial version of Cooley-Tukey FFT in plain C.
//

#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SERIAL_BLOCK_SIZE 1024
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
	float *data_in = depv[0].ptr;
	float *data_real = depv[1].ptr;
	float *data_imag = depv[2].ptr;
	ocrGuid_t dataInGuid = depv[0].guid;
	ocrGuid_t dataRealGuid = depv[1].guid;
	ocrGuid_t dataImagGuid = depv[2].guid;
	u64 N = paramv[3];
	u64 step = paramv[4];
	u64 offset = paramv[5];
	u64 x_in_offset = paramv[6];
	float *x_in = (float*)data_in;
	float *X_real = (float*)(data_real+offset);
	float *X_imag = (float*)(data_imag+offset);
	//PRINTF("Step %d offset: %d N*step: %d\n", step, offset, N*step);
	
	//if(step == 1) {
	//for(i=0;i<N;i++) {
	//	PRINTF("%d %f \n",i,x_in[i]);
	//}
	//}

	//ditfft2(X_real, X_imag, x_in, N, 1);
	if(N <= SERIAL_BLOCK_SIZE) {
		ditfft2(X_real, X_imag, x_in+x_in_offset, N, step);
		//X_real[0] = x_in[x_in_offset];
		//X_imag[0] = 0;
	} else {
		// DFT even side
		u64 childParamv[7] = { startGuid, endGuid, endSlaveGuid, N/2, 2 * step, 0 + offset, x_in_offset };
		u64 childParamv2[7] = { startGuid, endGuid, endSlaveGuid, N/2, 2 * step, N/2 + offset, x_in_offset + step };
		
		//PRINTF("Creating children of size %d\n",N/2);
		ocrGuid_t edtGuid, edtGuid2, endEdtGuid, finishEventGuid, finishEventGuid2;

		// Divide-and-conquer child EDTs will not need finish event guids
		ocrGuid_t childDependencies[3] = { dataInGuid, dataRealGuid, dataImagGuid };

		ocrEdtCreate(&edtGuid, startGuid, EDT_PARAM_DEF, childParamv, EDT_PARAM_DEF, childDependencies, EDT_PROP_FINISH, NULL_GUID, &finishEventGuid);
		ocrEdtCreate(&edtGuid2, startGuid, EDT_PARAM_DEF, childParamv2, EDT_PARAM_DEF, childDependencies, EDT_PROP_FINISH, NULL_GUID, &finishEventGuid2);
		//PRINTF("finishEventGuid after create: %lu\n",finishEventGuid);
		//ditfft2(X_real, X_imag, x_in, N/2, 2 * step);
		//ditfft2(X_real+N/2, X_imag+N/2, x_in+step, N/2, 2 * step);

		ocrGuid_t endDependencies[5] = { dataInGuid, dataRealGuid, dataImagGuid, finishEventGuid, finishEventGuid2 };

		// Do calculations after having divided and conquered
		ocrEdtCreate(&endEdtGuid, endGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, childDependencies, EDT_PROP_FINISH, NULL_GUID, NULL);
	}

	//PRINTF("Task with size %d completed\n",N);
	return NULL_GUID;
}

ocrGuid_t fftEndEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	int i;
	ocrGuid_t startGuid = paramv[0];
	ocrGuid_t endGuid = paramv[1];
	ocrGuid_t endSlaveGuid = paramv[2];
	float *data_in = depv[0].ptr;
	float *data_real = depv[1].ptr;
	float *data_imag = depv[2].ptr;
	ocrGuid_t dataGuids[3] = { depv[0].guid, depv[1].guid, depv[2].guid };
	u64 N = paramv[3];
	u64 step = paramv[4];
	u64 offset = paramv[5];
	float *x_in = (float*)data_in+offset;
	float *X_real = (float*)(data_real+offset);
	float *X_imag = (float*)(data_imag+offset);
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

			ocrEdtCreate(slaveGuids+i, endSlaveGuid, EDT_PARAM_DEF, slaveParamv+i*5, EDT_PARAM_DEF, dataGuids, EDT_PROP_NONE, NULL_GUID, NULL);
		}
	} else {
		slaveGuids = malloc(sizeof(ocrGuid_t));
		slaveParamv = malloc(sizeof(u64*) * 5);
		
		slaveParamv[0] = N;
		slaveParamv[1] = step;
		slaveParamv[2] = offset;
		slaveParamv[3] = 0;
		slaveParamv[4] = N/2;

		ocrEdtCreate(slaveGuids, endSlaveGuid, EDT_PARAM_DEF, slaveParamv, EDT_PARAM_DEF, dataGuids, EDT_PROP_NONE, NULL_GUID, NULL);
	}
	//int k;
	//for(k=0;k<N/2;k++) {
	//	float t_real = X_real[k];
	//	float t_imag = X_imag[k];
	//	double twiddle_real;
	//	double twiddle_imag;
	//	twiddle_imag = sin(-2 * M_PI * k / N);
	//	twiddle_real = cos(-2 * M_PI * k / N);
	//	float xr = X_real[k+N/2];
	//	float xi = X_imag[k+N/2];

	//	// (a+bi)(c+di) = (ac - bd) + (bc + ad)i
	//	X_real[k] = t_real +
	//		(twiddle_real*xr - twiddle_imag*xi);
	//	X_imag[k] = t_imag +
	//		(twiddle_imag*xr + twiddle_real*xi);
	//	X_real[k+N/2] = t_real -
	//		(twiddle_real*xr - twiddle_imag*xi);
	//	X_imag[k+N/2] = t_imag -
	//		(twiddle_imag*xr + twiddle_real*xi);
	//}
	//sleep(1);
	return NULL_GUID;
}
ocrGuid_t fftEndSlaveEdt(u32 paramc, u64 *paramv, u32 depc, ocrEdtDep_t depv[]) {
	int i;
	float *data_in = depv[0].ptr;
	float *data_real = depv[1].ptr;
	float *data_imag = depv[2].ptr;
	u64 N = paramv[0];
	u64 step = paramv[1];
	u64 offset = paramv[2];
	float *x_in = (float*)data_in+offset;
	float *X_real = (float*)(data_real+offset);
	float *X_imag = (float*)(data_imag+offset);
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
//		X_real[k] = t_real +
//			(twiddle_real*xr - twiddle_imag*xi);
//		X_imag[k] = t_imag +
//			(twiddle_imag*xr + twiddle_real*xi);
//		X_real[k+N/2] = t_real -
//			(twiddle_real*xr - twiddle_imag*xi);
//		X_imag[k+N/2] = t_imag -
//			(twiddle_imag*xr + twiddle_real*xi);
	}
//	sleep(1);

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
	float *data_in = depv[1].ptr;
	float *data_real = depv[2].ptr;
	float *data_imag = depv[3].ptr;
	u64 N = paramv[3];
	u64 step = paramv[4];
	u64 offset = paramv[5];
	float *x_in = (float*)data_in+offset;
	float *X_real = (float*)(data_real+offset);
	float *X_imag = (float*)(data_imag+offset);

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

	ocrGuid_t startTempGuid,endTempGuid,printTempGuid,endSlaveTempGuid;
	ocrEdtTemplateCreate(&startTempGuid, &fftStartEdt, 7, 3);
	ocrEdtTemplateCreate(&endTempGuid, &fftEndEdt, 7, 5);
	ocrEdtTemplateCreate(&endSlaveTempGuid, &fftEndSlaveEdt, 5, 3);
	ocrEdtTemplateCreate(&printTempGuid, &finalPrintEdt, 7, 4);
	
	int N = pow(2,strtol(getArgv(depv[0].ptr,1), NULL, 10));
	int iterations = ( argc > 2 ? strtol(getArgv(depv[0].ptr,2), NULL, 10) : 1 );

	// x_in, X_real, and X_imag in a contiguous block
	float *x_in;
	float *X_real;
	float *X_imag;
	ocrGuid_t dataInGuid,dataRealGuid,dataImagGuid;
	// TODO: OCR cannot handle large datablocks
	DBCREATE(&dataInGuid, (void **) &x_in, sizeof(float) * N, 0, NULL_GUID, NO_ALLOC);
	DBCREATE(&dataRealGuid, (void **) &X_real, sizeof(float) * N, 0, NULL_GUID, NO_ALLOC);
	DBCREATE(&dataImagGuid, (void **) &X_imag, sizeof(float) * N, 0, NULL_GUID, NO_ALLOC);
	printf("3 Datablocks of size %lu (N=%lu) created\n",sizeof(float)*N,N);
	
	for(i=0;i<N;i++) {
		x_in[i] = 0;
		X_real[i] = 0;
		X_imag[i] = 0;
	}
	x_in[1] = 1;
	x_in[3] = -3;
	x_in[4] = 8;
	x_in[5] = 9;
	x_in[6] = 1;


	u64 edtParamv[7] = { startTempGuid, endTempGuid, endSlaveTempGuid, N, 1 /* step size */, 0 /* offset */, 0 /* x_in_offset */ };
	
	// Create an EDT out of the EDT template
	ocrGuid_t edtGuid, printEdtGuid, edtEventGuid;
	ocrEdtCreate(&edtGuid, startTempGuid, EDT_PARAM_DEF, edtParamv, EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &edtEventGuid);

	ocrGuid_t finishDependencies[4] = { edtEventGuid, dataInGuid, dataRealGuid, dataImagGuid };
	ocrEdtCreate(&printEdtGuid, printTempGuid, EDT_PARAM_DEF, edtParamv, EDT_PARAM_DEF, finishDependencies, EDT_PROP_NONE, NULL_GUID, NULL);
	
	ocrAddDependence(dataInGuid, edtGuid, 0, DB_MODE_ITW);
	ocrAddDependence(dataRealGuid, edtGuid, 1, DB_MODE_ITW);
	ocrAddDependence(dataImagGuid, edtGuid, 2, DB_MODE_ITW);

	return NULL_GUID;
}
