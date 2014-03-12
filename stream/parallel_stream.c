#define __OCR__
#include "ocr.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h> 

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

#ifndef STREAM_ARRAY_SIZE
#	define STREAM_ARRAY_SIZE	100
#endif

#ifndef NSPLIT
#	define NSPLIT	10
#endif

#ifndef CHUNK
#	define CHUNK	STREAM_ARRAY_SIZE / NSPLIT
#endif

#ifndef NTIMES
#	define NTIMES	5
#endif

#ifndef SCALAR
#	define SCALAR	3.0
#endif

/*
	NOTES: 
	One large datablock is used with ranges of indexes corresponding to the following:
			  a --> [0 to (STREAM_ARRAY_SIZE - 1)]
			  b --> [STREAM_ARRAY_SIZE to (2 * STREAM_ARRAY_SIZE - 1)]
			  c --> [2 * STREAM_ARRAY_SIZE to (3 * STREAM_ARRAY_SIZE - 1)]
		timings --> [3 * STREAM_ARRAY_SIZE to (3 * STREAM_ARRAY_SIZE + NTIMES - 1)]
*/

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

ocrGuid_t pipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 begin = paramv[paramc - 2];
	u64 end = paramv[paramc - 1];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] = mysecond();

	for (i = begin; i < end; i++)
		data[2 * STREAM_ARRAY_SIZE + i] = data[i];
	// PRINTF("FINISHED COPY\n");

	for (i = begin; i < end; i++)
		data[STREAM_ARRAY_SIZE + i] = SCALAR * data[2 * STREAM_ARRAY_SIZE + i];
	// PRINTF("FINISHED SCALE\n");

	for (i = begin; i < end; i++)
		data[2 * STREAM_ARRAY_SIZE + i] = data[i] + data[STREAM_ARRAY_SIZE + i];
	// PRINTF("FINISHED ADD\n");

	for (i = begin; i < end; i++)
		data[i] = data[STREAM_ARRAY_SIZE + i] + SCALAR * data[2 * STREAM_ARRAY_SIZE + i];
	// PRINTF("FINISHED TRIAD\n");

	data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1] = mysecond() - data[3 * STREAM_ARRAY_SIZE + paramv[0] - 1];

	// PRINTF("FINISHED PIPELINE\n");
	return NULL_GUID;
}

ocrGuid_t pipeExecEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, nparamv[6] = {0};
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
	ocrGuid_t pipelineTemplateGuid = paramv[3];
	ocrGuid_t pipelineGuid;

	for (i = 0; i < paramc; i++)
		nparamv[i] = paramv[i];

	// Spawn pipeline children operating on CHUNK amounts of data
	for (i = 0; i < NSPLIT; i++) {
		nparamv[paramc] = i * CHUNK;
		nparamv[paramc + 1] = (i + 1) * CHUNK;
		ocrEdtCreate(&pipelineGuid, pipelineTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence(dataGuid, pipelineGuid, 0, DB_MODE_ITW);
	}
	return NULL_GUID;
}

ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
	ocrGuid_t pipeExecTemplateGuid = paramv[1];
	ocrGuid_t nextIterTemplateGuid = paramv[2];
	ocrGuid_t pipeExecGuid, pipeExecDone, nextIterGuid;

	ocrEdtCreate(&pipeExecGuid, pipeExecTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &pipeExecDone);

	// Setup next iteration
	if (paramv[0] < NTIMES) {
		paramv[0] += 1;
		ocrEdtCreate(&nextIterGuid, nextIterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
					 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence(dataGuid, nextIterGuid, 0, DB_MODE_EW);
		ocrAddDependence(pipeExecDone, nextIterGuid, 1, DB_MODE_RO);
	} 

	// Dependencies for pipeline
	ocrAddDependence(dataGuid, pipeExecGuid, 0, DB_MODE_EW);
	return NULL_GUID;
}

ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE a, ai, b, bi, c, ci;
	STREAM_TYPE totalsum, timingsum[NTIMES];
	a = b = c = 0;

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

	for (i = 0; i < NSPLIT; i++) {
		a += data[i];
		b += data[STREAM_ARRAY_SIZE + i];
		c += data[2 * STREAM_ARRAY_SIZE + i];
		for (j = 0; j < NTIMES; j++) {
			timingsum[j] += data[3 * STREAM_ARRAY_SIZE + j];
		}
	}

	PRINTF("Timing Results:\n");
	for (i = 0; i < NTIMES; i++) {
		PRINTF("TRIAL %d: %f s\n", i + 1, timingsum[i]);
		totalsum += timingsum[i];
	}
	PRINTF("AVERAGE Time Per CHUNK: %f s\n", totalsum / NTIMES / NSPLIT);
	PRINTF("AVERAGE Time Per Trial: %f s\n", totalsum / NTIMES);

	PRINTF("After %d Iterations:\n", NTIMES);
	if ((NSPLIT * ai - a + NSPLIT * bi - b + NSPLIT * ci - c) == 0)
		PRINTF("No differences between expected and actual\n");
	else  
		PRINTF("Expected a: %f, Actual a: %f\n"
			   "Expected b: %f, Actual b: %f\n"
			   "Expected c: %f, Actual c: %f\n"
			   "Note: Expected values are multiplied by NSPLIT (i.e. a * NSPLIT)\n"
				"     Actual values are sum of NSPLIT number of first elements", 
				NSPLIT * ai, a, NSPLIT * bi, b, NSPLIT * ci, c);

	ocrShutdown();
	return NULL_GUID;
}

ocrGuid_t mainEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, nparamc = 4, pparamc = 2;
	STREAM_TYPE * dataArray;
	ocrGuid_t dataGuid, iterTemplateGuid, iterGuid, iterDone, resultsTemplateGuid, resultsGuid, 
			  pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, 0, 2);
	ocrEdtTemplateCreate(&pipeExecTemplateGuid, &pipeExecEdt, nparamc, 1);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, 2);
	ocrEdtTemplateCreate(&pipelineTemplateGuid, &pipelineEdt, nparamc + pparamc, 1);

	u64 nparamv[4] = {1, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};

	// Formatting datablock
	DBCREATE(&dataGuid,(void **) &dataArray, sizeof(STREAM_TYPE) * (3 * STREAM_ARRAY_SIZE + NTIMES), 0, NULL_GUID, NO_ALLOC);
	for (i = 0; i < STREAM_ARRAY_SIZE; i++){
		dataArray[i] = 1.0;
		dataArray[STREAM_ARRAY_SIZE + i] = 2.0;
		dataArray[2 * STREAM_ARRAY_SIZE + i] = 0.0;
	}
	for (i = 3 * STREAM_ARRAY_SIZE; i < NTIMES; i++)
		dataArray[i] = 0.0;

	// Create iterator and results EDTs
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &iterDone);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, NULL_GUID, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for iterator and results
	ocrAddDependence(dataGuid, resultsGuid, 0, DB_MODE_RO);
	ocrAddDependence(iterDone, resultsGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, iterGuid, 0, DB_MODE_EW);
	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}