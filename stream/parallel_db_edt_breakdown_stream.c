#define __OCR__
#include "ocr.h"
#include <stdio.h>
#include <unistd.h>
#include <sys/time.h> 

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

#ifndef STREAM_ARRAY_SIZE
#	define STREAM_ARRAY_SIZE	10000000
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
	NSPLIT datablocks are used with ranges of indexes corresponding to the following:
			  a --> [0 to (CHUNK - 1)]
			  b --> [CHUNK to (2 * CHUNK - 1)]
			  c --> [2 * CHUNK to (3 * CHUNK - 1)]
		timings --> [3 * CHUNK to (3 * CHUNK + NTIMES - 1)]
*/

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

ocrGuid_t copyEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	data[3 * CHUNK + paramv[0] - 1] = mysecond();
	for (i = 0; i < CHUNK; i++)
		data[2 * CHUNK + i] = data[i];
	data[3 * CHUNK + paramv[0] - 1] = mysecond() - data[3 * CHUNK + paramv[0] - 1];
	// PRINTF("FINISHED COPY\n");
	return NULL_GUID;
}

ocrGuid_t scaleEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < CHUNK; i++)
		data[CHUNK + i] = SCALAR * data[2 * CHUNK + i];
	time = mysecond() - time;
	data[3 * CHUNK + paramv[0] - 1] += time;
	// PRINTF("FINISHED SCALE\n");
	return NULL_GUID;
}

ocrGuid_t addEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < CHUNK; i++)
		data[2 * CHUNK + i] = data[i] + data[CHUNK + i];
	time = mysecond() - time;
	data[3 * CHUNK + paramv[0] - 1] += time;
	// PRINTF("FINISHED ADD\n");
	return NULL_GUID;
}

ocrGuid_t triadEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < CHUNK; i++)
		data[i] = data[CHUNK + i] + SCALAR * data[2 * CHUNK + i];
	time = mysecond() - time;
	data[3 * CHUNK + paramv[0] - 1] += time;
	// PRINTF("FINISHED TRIAD\n");
	return NULL_GUID;
}

ocrGuid_t pipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	ocrGuid_t copyGuid, copyDone, scaleGuid, scaleDone, addGuid, addDone, triadGuid, triadDone, nextIterGuid;
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
	ocrGuid_t copyTemplateGuid = paramv[4];
	ocrGuid_t scaleTemplateGuid = paramv[5];
	ocrGuid_t addTemplateGuid = paramv[6];
	ocrGuid_t triadTemplateGuid = paramv[7];

	// EDTs for vector operations
	ocrEdtCreate(&copyGuid, copyTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &copyDone);
	ocrEdtCreate(&scaleGuid, scaleTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &scaleDone);
	ocrEdtCreate(&addGuid, addTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &addDone);
	ocrEdtCreate(&triadGuid, triadTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &triadDone);

	// Dependencies for vector operations
	ocrAddDependence(dataGuid, triadGuid, 0, DB_MODE_EW);
	ocrAddDependence(addDone, triadGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, addGuid, 0, DB_MODE_EW);
	ocrAddDependence(scaleDone, addGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, scaleGuid, 0, DB_MODE_EW);
	ocrAddDependence(copyDone, scaleGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, copyGuid, 0, DB_MODE_EW);

	// PRINTF("FINISHED PIPELINE\n");
	return NULL_GUID;
}

ocrGuid_t pipeExecEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	ocrGuid_t pipelineTemplateGuid = paramv[3];
	ocrGuid_t pipelineGuid;

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
	ocrGuid_t pipeExecTemplateGuid = paramv[1];
	ocrGuid_t nextIterTemplateGuid = paramv[2];
	ocrGuid_t pipeExecGuid, pipeExecDone, nextIterGuid;

	ocrEdtCreate(&pipeExecGuid, pipeExecTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &pipeExecDone);

	// PRINTF("FINISHED ITER\n");
	if (paramv[0] < NTIMES) {
		paramv[0] += 1;
		// EDT for next iteration
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

ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	STREAM_TYPE a, b, c, timingsum[NTIMES];
	STREAM_TYPE ai, bi, ci;
	STREAM_TYPE totalsum = 0.0;
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

ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j, nparamc = 8;
	STREAM_TYPE * chunkArray;
	ocrGuid_t dataGuids[NSPLIT] = {0};
	ocrGuid_t iterTemplateGuid, iterGuid, iterDone, resultsTemplateGuid, resultsGuid, pipeExecTemplateGuid, pipelineTemplateGuid, 
			  copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, NSPLIT);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, 0, NSPLIT + 1);

	ocrEdtTemplateCreate(&pipeExecTemplateGuid, &pipeExecEdt, nparamc, NSPLIT);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, NSPLIT + 1);
	ocrEdtTemplateCreate(&pipelineTemplateGuid, &pipelineEdt, nparamc, 1);

	ocrEdtTemplateCreate(&copyTemplateGuid, &copyEdt, nparamc, 1);
	ocrEdtTemplateCreate(&scaleTemplateGuid, &scaleEdt, nparamc, 2);
	ocrEdtTemplateCreate(&addTemplateGuid, &addEdt, nparamc, 2);
	ocrEdtTemplateCreate(&triadTemplateGuid, &triadEdt, nparamc, 2);

	u64 nparamv[8]= {1, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid}; 

	// Preparing datablock
	for (i = 0; i < NSPLIT; i++) {
		ocrGuid_t chunkGuid; 
		DBCREATE(&chunkGuid,(void **) &chunkArray, sizeof(STREAM_TYPE)*(3 * CHUNK + NTIMES), 0, NULL_GUID, NO_ALLOC);
		for (j = 0; j < CHUNK; j++) {
			chunkArray[j] = 1.0;
			chunkArray[CHUNK + j] = 2.0;
			chunkArray[2 * CHUNK + j] = 0.0;
		}
		dataGuids[i] = chunkGuid;
	}

	// Create iterator and results EDTs
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &iterDone);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, NULL_GUID, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for iterator and results
	for (i = 0; i < NSPLIT; i++)
		ocrAddDependence(dataGuids[i], resultsGuid, i, DB_MODE_RO);
	ocrAddDependence(iterDone, resultsGuid, NSPLIT, DB_MODE_RO);
	for (i = 0; i < NSPLIT; i++)
		ocrAddDependence(dataGuids[i], iterGuid, i, DB_MODE_EW);
	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}