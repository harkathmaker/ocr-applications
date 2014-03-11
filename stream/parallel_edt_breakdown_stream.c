#define __OCR__
#include "ocr.h"

/* STREAM Headers */
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h> 

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

#ifndef NTIMES
#	define NTIMES	5
#endif

#ifndef STREAM_ARRAY_SIZE
#	define STREAM_ARRAY_SIZE	100
#endif

#ifndef SCALAR
#	define SCALAR	3.0
#endif

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

ocrGuid_t copyEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] = mysecond();
	for (i = 0; i < STREAM_ARRAY_SIZE; i++)
		data[2 * STREAM_ARRAY_SIZE + i] = data[i];
	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] = mysecond() - data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1];
	// PRINTF("FINISHED COPY\n");
	return NULL_GUID;
}

ocrGuid_t scaleEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < STREAM_ARRAY_SIZE; i++)
		data[STREAM_ARRAY_SIZE + i] = SCALAR * data[2 * STREAM_ARRAY_SIZE + i];
	time = mysecond() - time;
	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] += time;
	// PRINTF("FINISHED SCALE\n");
	return NULL_GUID;
}

ocrGuid_t addEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < STREAM_ARRAY_SIZE; i++)
		data[2 * STREAM_ARRAY_SIZE + i] = data[i] + data[STREAM_ARRAY_SIZE + i];
	time = mysecond() - time;
	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] += time;
	// PRINTF("FINISHED ADD\n");
	return NULL_GUID;
}

ocrGuid_t triadEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < STREAM_ARRAY_SIZE; i++)
		data[i] = data[STREAM_ARRAY_SIZE + i] + SCALAR * data[2 * STREAM_ARRAY_SIZE + i];
	time = mysecond() - time;
	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] += time;
	// PRINTF("FINISHED TRIAD\n");
	return NULL_GUID;
}

ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
	ocrGuid_t copyTemplateGuid = paramv[1];
	ocrGuid_t scaleTemplateGuid = paramv[2];
	ocrGuid_t addTemplateGuid = paramv[3];
	ocrGuid_t triadTemplateGuid = paramv[4];
	ocrGuid_t nextIterTemplateGuid = paramv[5];

	// EDTs for vector operations
	ocrGuid_t copyGuid, copyDone, scaleGuid, scaleDone, addGuid, addDone, triadGuid, triadDone, nextIterGuid;
	ocrEdtCreate(&copyGuid, copyTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &copyDone);
	ocrEdtCreate(&scaleGuid, scaleTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &scaleDone);
	ocrEdtCreate(&addGuid, addTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &addDone);
	ocrEdtCreate(&triadGuid, triadTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &triadDone);

	// PRINTF("FINISHED ITER\n");
	if (paramv[0] < NTIMES) {
		paramv[0] += 1;
		// EDT for next iteration
		ocrEdtCreate(&nextIterGuid, nextIterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
					 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence(dataGuid, nextIterGuid, 0, DB_MODE_EW);
		ocrAddDependence(triadDone, nextIterGuid, 1, DB_MODE_RO);
	}

	// Dependencies for vector operations
	ocrAddDependence(dataGuid, triadGuid, 0, DB_MODE_EW);
	ocrAddDependence(addDone, triadGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, addGuid, 0, DB_MODE_EW);
	ocrAddDependence(scaleDone, addGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, scaleGuid, 0, DB_MODE_EW);
	ocrAddDependence(copyDone, scaleGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, copyGuid, 0, DB_MODE_EW);
	return NULL_GUID;
}

ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE a = data[0];
	STREAM_TYPE b = data[STREAM_ARRAY_SIZE];
	STREAM_TYPE c = data[2 * STREAM_ARRAY_SIZE];
	STREAM_TYPE ai, bi, ci;
	STREAM_TYPE timingsum = 0.0;

	// Reproduce initializations
	ai = 1.0;
	bi = 2.0;
	ci = 0.0;

	// Execute timing loop
	for (i = 0; i < NTIMES; i++) {
		ci = ai;
		bi = SCALAR * ci;
		ci = ai + bi;
		ai = bi + SCALAR * ci;
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
	u64 i;
	u64 nparamc = 6;
	ocrGuid_t iterTemplateGuid, resultsTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, 0, 2);
	ocrEdtTemplateCreate(&copyTemplateGuid, &copyEdt, nparamc, 1);
	ocrEdtTemplateCreate(&scaleTemplateGuid, &scaleEdt, nparamc, 2);
	ocrEdtTemplateCreate(&addTemplateGuid, &addEdt, nparamc, 2);
	ocrEdtTemplateCreate(&triadTemplateGuid, &triadEdt, nparamc, 2);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, 2);

	u64 nparamv[6]= {1, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid}; 

	// Preparing datablock
	ocrGuid_t dataGuid;
	STREAM_TYPE * dataArray;
	DBCREATE(&dataGuid,(void **) &dataArray, sizeof(STREAM_TYPE)*(3 * STREAM_ARRAY_SIZE + NTIMES), 0, NULL_GUID, NO_ALLOC);
	for (i = 0; i < STREAM_ARRAY_SIZE; i++){
		dataArray[i] = 1.0;
		dataArray[STREAM_ARRAY_SIZE + i] = 2.0;
		dataArray[2 * STREAM_ARRAY_SIZE + i] = 0.0;
	}
	for (i = 3 * STREAM_ARRAY_SIZE; i < NTIMES; i++)
		dataArray[i] = 0.0;

	// Create iterator and results EDTs
	ocrGuid_t iterGuid, iterOutput, resultsGuid;
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &iterOutput);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, NULL_GUID, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for iterator and results
	ocrAddDependence(dataGuid, resultsGuid, 0, DB_MODE_RO);
	ocrAddDependence(iterOutput, resultsGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, iterGuid, 0, DB_MODE_EW);
	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}