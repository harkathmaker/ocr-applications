#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h> 

#ifndef __OCR__
#	define __OCR__
#	include "ocr.h"
#endif

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

/*
	NOTES: 
	split datablocks are used with ranges of indexes corresponding to the following:
			  a --> [0 to (chunk - 1)]
			  b --> [chunk to (2 * chunk - 1)]
			  c --> [2 * chunk to (3 * chunk - 1)]
		timings --> [3 * chunk to (3 * chunk + iterations - 1)]
*/

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void export_csv(char * name, u64 db_size, u64 iterations, STREAM_TYPE * trials, STREAM_TYPE avg) {
	char path[150];
	//strcpy(path, "./results/");
	strcpy(path, name);
	printf("path=%s\n", path);
	FILE * f = fopen(path, "a");
	if (f == NULL) {
		PRINTF("Error creating export file.");
		exit(1);
	}

	u64 i;
	for (i = 0; i < iterations; i++) 
		fprintf(f, "%llu %f\n", db_size, trials[i]);
	//fprintf(f, "%f\n", avg);

	fclose(f);
	return;
}

//                   0  1        2           3      4      5       6                     7 
// u64 nparamv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
//                   8                     9                 10                 11               12
//					 pipelineTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid}; 
ocrGuid_t copyEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	data[3 * chunk + paramv[0] - 1] = mysecond();
	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i];
	data[3 * chunk + paramv[0] - 1] = mysecond() - data[3 * chunk + paramv[0] - 1];
	// PRINTF("FINISHED COPY\n");
	return NULL_GUID;
}

ocrGuid_t scaleEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE scalar = paramv[5];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < chunk; i++)
		data[chunk + i] = scalar * data[2 * chunk + i];
	time = mysecond() - time;
	data[3 * chunk + paramv[0] - 1] += time;
	// PRINTF("FINISHED SCALE\n");
	return NULL_GUID;
}

ocrGuid_t addEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i] + data[chunk + i];
	time = mysecond() - time;
	data[3 * chunk + paramv[0] - 1] += time;
	// PRINTF("FINISHED ADD\n");
	return NULL_GUID;
}

ocrGuid_t triadEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE scalar = paramv[5];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE time = mysecond();
	for (i = 0; i < chunk; i++)
		data[i] = data[chunk + i] + scalar * data[2 * chunk + i];
	time = mysecond() - time;
	data[3 * chunk + paramv[0] - 1] += time;
	// PRINTF("FINISHED TRIAD\n");
	return NULL_GUID;
}

//                   0  1        2           3      4      5       6                     7 
// u64 nparamv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
//                   8                     9                 10                 11               12
//					 pipelineTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid}; 
ocrGuid_t pipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	ocrGuid_t copyTemplateGuid = paramv[9];
	ocrGuid_t scaleTemplateGuid = paramv[10];
	ocrGuid_t addTemplateGuid = paramv[11];
	ocrGuid_t triadTemplateGuid = paramv[12];
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
	ocrGuid_t copyGuid, copyDone, scaleGuid, scaleDone, addGuid, addDone, triadGuid, triadDone, nextIterGuid;

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
	ocrAddDependence(dataGuid, triadGuid, 0, DB_MODE_ITW);
	ocrAddDependence(addDone, triadGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, addGuid, 0, DB_MODE_ITW);
	ocrAddDependence(scaleDone, addGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, scaleGuid, 0, DB_MODE_ITW);
	ocrAddDependence(copyDone, scaleGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, copyGuid, 0, DB_MODE_ITW);

	// PRINTF("FINISHED PIPELINE\n");
	return NULL_GUID;
}

//                   0  1        2           3      4      5       6                     7 
// u64 nparamv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
//                   8                     9                 10                 11               12
//					 pipelineTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid}; 
ocrGuid_t pipeExecEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 split = paramv[3];
	ocrGuid_t pipelineTemplateGuid = paramv[8];
	ocrGuid_t pipelineGuid;

	// Spawn pipeline children operating on chunk amounts of data
	for (i = 0; i < split; i++) {
		ocrEdtCreate(&pipelineGuid, pipelineTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence((ocrGuid_t) depv[i].guid, pipelineGuid, 0, DB_MODE_ITW);
	}
	// PRINTF("FINISHED PIPE EXEC\n");
	return NULL_GUID;
}

//                   0  1        2           3      4      5       6                     7 
// u64 nparamv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
//                   8                     9                 10                 11               12
//					 pipelineTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid}; 
ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 iterations = paramv[2];
	u64 split = paramv[3];
	ocrGuid_t pipeExecTemplateGuid = paramv[6];
	ocrGuid_t nextIterTemplateGuid = paramv[7];
	ocrGuid_t pipeExecGuid, pipeExecDone, nextIterGuid;

	ocrEdtCreate(&pipeExecGuid, pipeExecTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &pipeExecDone);

	// PRINTF("FINISHED ITER %llu \n", iterations);
	if (paramv[0] < iterations) {
		paramv[0] += 1;;
		// EDT for next iteration
		ocrEdtCreate(&nextIterGuid, nextIterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
					 EDT_PROP_FINISH, NULL_GUID, NULL);
		for (i = 0; i < split; i++)
			ocrAddDependence((ocrGuid_t) depv[i].guid, nextIterGuid, i, DB_MODE_ITW);
		ocrAddDependence(pipeExecDone, nextIterGuid, split, DB_MODE_RO);
	}

	// Dependencies for pipeline
	for (i = 0; i < split; i++)
		ocrAddDependence((ocrGuid_t) depv[i].guid, pipeExecGuid, i, DB_MODE_ITW);
	return NULL_GUID;
}

//                   0        1           2      3      4       5       6
// u64 rparamv[7] = {db_size, iterations, split, chunk, verify, scalar, verbose}
ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	u64 db_size = paramv[0];
	u64 iterations = paramv[1];
	u64 split = paramv[2];
	u64 chunk = paramv[3];
	int verify = (int) paramv[4];
	int verbose = (int) paramv[6];
	STREAM_TYPE totalsum = 0.0, timingsum[iterations], avg;

	for (i = 0; i < split; i++) {
		STREAM_TYPE * cur = (STREAM_TYPE *) depv[i].ptr;
		for (j = 0; j < iterations; j++)
			timingsum[j] += cur[3 * chunk + j];
	}

	PRINTF("Timing Results:\n");
	for (i = 0; i < iterations; i++) {
		if (verbose)
			PRINTF("TRIAL %d: %f s\n", i + 1, timingsum[i]);
		totalsum += timingsum[i];
	}
	avg = totalsum / iterations;
	PRINTF("AVERAGE Time Per chunk: %f s\n", avg / split);
	PRINTF("AVERAGE Time Per Trial: %f s\n", avg);

	if (strcmp((char *) depv[split].ptr, "") != 0) 
		export_csv((char *) depv[split].ptr, db_size, iterations, timingsum, avg);

	if (verify) {
		STREAM_TYPE a = 0, ai, b = 0, bi, c = 0, ci;
		STREAM_TYPE scalar = paramv[5];

		// Sum actual values
		for (i = 0; i < split; i++) {
			STREAM_TYPE * cur = (STREAM_TYPE *) depv[i].ptr;
			a += cur[0];
			b += cur[chunk];
			c += cur[2 * chunk];
		}
		// Reproduce initializations
		ai = 1.0;
		bi = 2.0;
		ci = 0.0;

		// Execute timing loop
		for (i = 0; i < iterations; i++) {
			ci = ai;
			bi = scalar * ci;
			ci = ai + bi;
			ai = bi + scalar * ci;
		}
		
		PRINTF("After %d Iterations:\n", iterations);
		if ((split * ai - a + split * bi - b + split * ci - c) == 0)
			PRINTF("No differences between expected and actual\n");
		else  
			PRINTF("Expected a: %f, Actual a: %f\n"
				   "Expected b: %f, Actual b: %f\n"
				   "Expected c: %f, Actual c: %f\n", split* ai, a, split * bi, b, split * ci, c);
	}
	ocrShutdown();
	return NULL_GUID;
}

ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 c, i, j, argc = getArgc(depv[0].ptr);
	char * argv[argc];
	for(i = 0; i < argc; i++)
		argv[i] = getArgv(depv[0].ptr, i);

	// Getopt arguments
	u64 db_size;
	char efile[100] = "";
	u64 iterations;
	u64 split;
	u64 chunk;
	int verify;
	STREAM_TYPE scalar;
	int verbose;

	// Parse getopt commands, shutdown and exit if help is selected
	if (parseParallelOptions(argc, argv,  &db_size, efile, &iterations, &split, &chunk, &verify, &scalar, &verbose)) {
		ocrShutdown();
		return NULL_GUID;
	}

	u64 nparamc = 13, rparamc = 7;
	STREAM_TYPE * chunkArray;
	char * efileArray;
	ocrGuid_t dataGuids[split];
	ocrGuid_t iterTemplateGuid, iterGuid, iterDone, efileGuid, resultsTemplateGuid, resultsGuid, pipeExecTemplateGuid, pipelineTemplateGuid, 
			  copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, split);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, rparamc, split + 2);

	ocrEdtTemplateCreate(&pipeExecTemplateGuid, &pipeExecEdt, nparamc, split);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, split + 1);
	ocrEdtTemplateCreate(&pipelineTemplateGuid, &pipelineEdt, nparamc, 1);

	ocrEdtTemplateCreate(&copyTemplateGuid, &copyEdt, nparamc, 1);
	ocrEdtTemplateCreate(&scaleTemplateGuid, &scaleEdt, nparamc, 2);
	ocrEdtTemplateCreate(&addTemplateGuid, &addEdt, nparamc, 2);
	ocrEdtTemplateCreate(&triadTemplateGuid, &triadEdt, nparamc, 2);

	u64 nparamv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
					 pipelineTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid}; 
	u64 rparamv[7] = {db_size, iterations, split, chunk, verify, scalar, verbose};

	// Preparing datablock
	for (i = 0; i < split; i++) {
		ocrGuid_t chunkGuid; 
		DBCREATE(&chunkGuid,(void **) &chunkArray, sizeof(STREAM_TYPE)*(3 * chunk + iterations), 0, NULL_GUID, NO_ALLOC);
		for (j = 0; j < chunk; j++) {
			chunkArray[j] = 1.0;
			chunkArray[chunk + j] = 2.0;
			chunkArray[2 * chunk + j] = 0.0;
		}
		dataGuids[i] = chunkGuid;
	}

	// Create datablock for export file name
	DBCREATE(&efileGuid, (void **) &efileArray, sizeof(char) * sizeof(efile), 0, NULL_GUID, NO_ALLOC);
	strcpy(efileArray, efile);

	// Create iterator and results EDTs
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &iterDone);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, rparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for iterator and results
	for (i = 0; i < split; i++)
		ocrAddDependence(dataGuids[i], resultsGuid, i, DB_MODE_RO);
	ocrAddDependence(efileGuid, resultsGuid, split, DB_MODE_RO);
	ocrAddDependence(iterDone, resultsGuid, split + 1, DB_MODE_RO);
	for (i = 0; i < split; i++)
		ocrAddDependence(dataGuids[i], iterGuid, i, DB_MODE_ITW);
	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}