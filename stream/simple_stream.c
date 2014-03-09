#define __OCR__
#include "ocr.h"

/* STREAM Headers */
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h> 

#ifndef STREAM_ARRAY_SIZE
#	define STREAM_ARRAY_SIZE	100
#endif

#ifndef NTIMES
#	define NTIMES	5
#endif

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	u64 i;
	STREAM_TYPE scalar = 3.0;

	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] = mysecond();

	for (i = 0; i < STREAM_ARRAY_SIZE; i++)
		data[2 * STREAM_ARRAY_SIZE + i] = data[i];
	// PRINTF("FINISHED COPY\n");

	for (i = 0; i < STREAM_ARRAY_SIZE; i++)
		data[STREAM_ARRAY_SIZE + i] = scalar * data[2 * STREAM_ARRAY_SIZE + i];
	// PRINTF("FINISHED SCALE\n");

	for (i = 0; i < STREAM_ARRAY_SIZE; i++)
		data[2 * STREAM_ARRAY_SIZE + i] = data[i] + data[STREAM_ARRAY_SIZE + i];
	// PRINTF("FINISHED ADD\n");

	for (i = 0; i < STREAM_ARRAY_SIZE; i++)
		data[i] = data[STREAM_ARRAY_SIZE + i] + scalar * data[2 * STREAM_ARRAY_SIZE + i];
	// PRINTF("FINISHED TRIAD\n");

	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] = mysecond() - data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1];

	// PRINTF("FINISHED ITER\n");

	if (paramv[0] < paramv[1]) {
		ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
		paramv[0] += 1;
		
		// EDT Template for next iteration
		ocrGuid_t iterTemplateGuid;
		ocrEdtTemplateCreate(&iterTemplateGuid, iterEdt, sizeof(paramv), 1);

		// EDT for next iteration
		ocrGuid_t iterGuid;
		ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
			EDT_PROP_FINISH, NULL_GUID, NULL);

		ocrAddDependence(dataGuid, iterGuid, 0, DB_MODE_EW);
	}
	return NULL_GUID;
}

ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[1].ptr;
	STREAM_TYPE a = data[0];
	STREAM_TYPE b = data[STREAM_ARRAY_SIZE];
	STREAM_TYPE c = data[2 * STREAM_ARRAY_SIZE];
	STREAM_TYPE ai, bi, ci, scalar;
	STREAM_TYPE timingsum = 0.0;

	// Reproduce initializations
	ai = 1.0;
	bi = 2.0;
	ci = 0.0;

	// Execute timing loop
	scalar = 3.0;
	for (i = 0; i < NTIMES; i++) {
		ci = ai;
		bi = scalar * ci;
		ci = ai + bi;
		ai = bi + scalar * ci;
	}

	PRINTF("Timing Results:\n");
	for (i = 0; i < NTIMES; i++) {
		STREAM_TYPE timing = data[3 * STREAM_ARRAY_SIZE + i];
		PRINTF("TRIAL %d: %f s\n", i + 1, timing);
		timingsum += timing;
	}
	PRINTF("AVERAGE: %f s\n", timingsum / NTIMES);

	PRINTF("After %d Iterations:\n", NTIMES);
	if ((ai - a + bi - b + ci - c) == 0)
		PRINTF("No differences between expected and actual\n");
	else  
		PRINTF("Expected a: %f, Actual a: %f\n"
			   "Expected b: %f, Actual b: %f\n"
			   "Expected c: %f, Actual c: %f\n", ai, a, bi, b, ci, c);

	ocrShutdown();
	return NULL_GUID;
}
ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, nparamv[2] = {1, NTIMES};

	// Preparing datablock
	STREAM_TYPE * dataArray;
	ocrGuid_t dataGuid;
	DBCREATE(&dataGuid,(void **) &dataArray, sizeof(ocrGuid_t)*(3 * STREAM_ARRAY_SIZE + NTIMES), 0, NULL_GUID, NO_ALLOC);
	for (i = 0; i < STREAM_ARRAY_SIZE; i++){
		dataArray[i] = 1.0;
		dataArray[STREAM_ARRAY_SIZE + i] = 2.0;
		dataArray[2 * STREAM_ARRAY_SIZE + i] = 0.0;
	}
	for (i = 3 * STREAM_ARRAY_SIZE; i < NTIMES; i++)
		dataArray[i] = 0.0;

	// EDT Template for iterator
	ocrGuid_t iterTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, iterEdt, sizeof(nparamv), 1);

	// EDT for iterator
	ocrGuid_t iterGuid, iterOutput;
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
		EDT_PROP_FINISH, NULL_GUID, &iterOutput);

	// Dependencies for iterator
	ocrAddDependence(dataGuid, iterGuid, 0, DB_MODE_EW);

	// EDT Template for results
	ocrGuid_t resultsTemplateGuid;
	ocrEdtTemplateCreate(&resultsTemplateGuid, resultsEdt, 0, 2);

	// EDT for results
	ocrGuid_t resultsGuid;
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, NULL_GUID, EDT_PARAM_DEF, NULL_GUID,
		EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for results
	ocrAddDependence(iterOutput, resultsGuid, 0, DB_MODE_RO);
	ocrAddDependence(dataGuid, resultsGuid, 1, DB_MODE_RO);

	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}