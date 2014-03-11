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
#	define NTIMES	10
#endif

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

#ifndef NSPLIT
#	define NSPLIT	10
#endif

#ifndef CHUNK
#	define CHUNK	STREAM_ARRAY_SIZE / NSPLIT
#endif

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	STREAM_TYPE a, b, c, timingsum[NTIMES];
	STREAM_TYPE ai, bi, ci, scalar;
	STREAM_TYPE totalsum = 0.0;
	a = b = c = 0;

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

	for (i = 0; i < NSPLIT; i++) {
		STREAM_TYPE * cur = (STREAM_TYPE *) depv[i].ptr;
		a += cur[0];
		b += cur[CHUNK];
		c += cur[2 * CHUNK];
		for (j = 0; j < NTIMES; j++)
			timingsum[j] += cur[3 * CHUNK + j];
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
			   "Expected c: %f, Actual c: %f\n", NSPLIT* ai, a, NSPLIT * bi, b, NSPLIT * ci, c);

	ocrShutdown();
	return NULL_GUID;
}

ocrGuid_t pipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE scalar = 3.0;

	data[3 * CHUNK + paramv[0] - 1] = mysecond();

	for (i = 0; i < CHUNK; i++)
		data[2 * CHUNK + i] = data[i];
	// PRINTF("FINISHED COPY\n");
	//for (i = 0; i < 3 * CHUNK; i++)
		//PRINTF("datablock[] = %f ", data[i]);

	for (i = 0; i < CHUNK; i++)
		data[CHUNK + i] = scalar * data[2 * CHUNK + i];
	// PRINTF("FINISHED SCALE\n");

	for (i = 0; i < CHUNK; i++)
		data[2 * CHUNK + i] = data[i] + data[CHUNK + i];
	// PRINTF("FINISHED ADD\n");

	for (i = 0; i < CHUNK; i++)
		data[i] = data[CHUNK + i] + scalar * data[2 * CHUNK + i];
	// PRINTF("FINISHED TRIAD\n");

	data[3 * CHUNK + paramv[0] - 1] = mysecond() - data[3 * CHUNK + paramv[0] - 1];

	// PRINTF("FINISHED PIPELINE\n");
	return NULL_GUID;
}

ocrGuid_t pipeExecEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	//ocrGuid_t * data = (ocrGuid_t *) depv[0].ptr;
	ocrGuid_t pipelineTemplateGuid = paramv[3];
	ocrGuid_t pipelineGuid;

	//                                    a                        b                                 c                          partial timings
	// Each element in dataGuid = [ (0, ... , CHUNK - 1), (CHUNK, ... , 2 * CHUNK - 1), (2 * CHUNK, ... ,3 * CHUNK - 1), (3 * CHUNK , ... , 3 * CHUNK + NTIMES - 1) ]
	// Size = 3 * CHUNK + NTIMES
	// Spawn pipeline children operating on CHUNK amounts of data
	for (i = 0; i < NSPLIT; i++) {
		ocrEdtCreate(&pipelineGuid, pipelineTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence((ocrGuid_t) depv[i].guid, pipelineGuid, 0, DB_MODE_ITW);
	}
	// PRINTF("FINISHED PIPE EXEC\n");
	return NULL_GUID;
}

ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	//ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
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
		for (i = 0; i < NSPLIT; i++)
			ocrAddDependence((ocrGuid_t) depv[i].guid, nextIterGuid, i, DB_MODE_EW);
		ocrAddDependence(pipeExecDone, nextIterGuid, NSPLIT, DB_MODE_RO);
	} 

	// Dependencies for pipeline
	for (i = 0; i < NSPLIT; i++)
		ocrAddDependence((ocrGuid_t) depv[i].guid, pipeExecGuid, i, DB_MODE_EW);
	// PRINTF("FINISHED ITER\n");
	return NULL_GUID;
}

ocrGuid_t mainEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	u64 nparamc = 4;

	ocrGuid_t iterTemplateGuid, resultsTemplateGuid, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, NSPLIT);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, 0, NSPLIT + 1);
	ocrEdtTemplateCreate(&pipeExecTemplateGuid, &pipeExecEdt, nparamc, NSPLIT);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, NSPLIT + 1);
	ocrEdtTemplateCreate(&pipelineTemplateGuid, &pipelineEdt, nparamc, 1);

	u64 nparamv[4] = {1, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};

	// Preparing datablock
	ocrGuid_t chunkGuid, dataGuids[NSPLIT] = {0};
	STREAM_TYPE * chunkArray;

	// 3 * STREAM_ARRAY_SIZE + NTIMES = (3 * CHUNK + NTIMES)* NSPLIT
	for (i = 0; i < NSPLIT; i++) {
		DBCREATE(&chunkGuid,(void **) &chunkArray, sizeof(STREAM_TYPE)*(3 * CHUNK + NTIMES), 0, NULL_GUID, NO_ALLOC);
		for (j = 0; j < CHUNK; j++) {
			chunkArray[j] = 1.0;
			chunkArray[CHUNK + j] = 2.0;
			chunkArray[2 * CHUNK + j] = 0.0;
		}
		dataGuids[i] = chunkGuid;
	}

	// Create iterator and results EDTs
	ocrGuid_t iterGuid, iterOutput, resultsGuid;
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &iterOutput);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, NULL_GUID, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for iterator and results

	for (i = 0; i < NSPLIT; i++)
		ocrAddDependence(dataGuids[i], resultsGuid, i, DB_MODE_RO);
	ocrAddDependence(iterOutput, resultsGuid, NSPLIT, DB_MODE_RO);
	for (i = 0; i < NSPLIT; i++)
		ocrAddDependence(dataGuids[i], iterGuid, i, DB_MODE_EW);

	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}