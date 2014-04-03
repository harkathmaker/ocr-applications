#include "options.h"
#include <stdio.h>
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
	One large datablock is used with ranges of indexes corresponding to the following:
			  a --> [0 to (db_size - 1)]
			  b --> [db_size to (2 * db_size - 1)]
			  c --> [2 * db_size to (3 * db_size - 1)]
		timings --> [3 * db_size to (3 * db_size + iterations - 1)]
*/

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void export_csv(char * name, u64 iterations, STREAM_TYPE * trials, STREAM_TYPE avg) {
	char path[150];
	strcpy(path, "./results/");
	strcat(path, name);
	FILE * f = fopen(path, "a");
	if (f == NULL) {
		PRINTF("Error creating export file.");
		exit(1);
	}

	u64 i;
	for (i = 0; i < iterations; i++) 
		fprintf(f, "%f, ", trials[i]);
	fprintf(f, "%f\n", avg);

	fclose(f);
	return;
}

ocrGuid_t pipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[1];
	u64 begin = paramv[paramc - 2];
	u64 end = paramv[paramc - 1];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[5];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	data[3 * db_size + paramv[0] - 1] = mysecond();

	for (i = begin; i < end; i++)
		data[2 * db_size + i] = data[i];
	// PRINTF("FINISHED COPY\n");

	for (i = begin; i < end; i++)
		data[db_size + i] = scalar * data[2 * db_size + i];
	// PRINTF("FINISHED SCALE\n");

	for (i = begin; i < end; i++)
		data[2 * db_size + i] = data[i] + data[db_size + i];
	// PRINTF("FINISHED ADD\n");

	for (i = begin; i < end; i++)
		data[i] = data[db_size + i] + scalar * data[2 * db_size + i];
	// PRINTF("FINISHED TRIAD\n");

	data[3 * db_size + paramv[0] - 1] = mysecond() - data[3 * db_size + paramv[0] - 1];

	// PRINTF("FINISHED PIPELINE\n");
	return NULL_GUID;
}

//                   0  1        2           3      4      5       6                     7                     8
// u64 nparamv[9] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};
ocrGuid_t pipeExecEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, nparamv[11] = {0};
	u64 split = paramv[3];
	u64 chunk = paramv[4];
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
	ocrGuid_t pipelineTemplateGuid = paramv[8];
	ocrGuid_t pipelineGuid;

	for (i = 0; i < paramc; i++)
		nparamv[i] = paramv[i];

	// Spawn pipeline children operating on chunk amounts of data
	for (i = 0; i < split; i++) {
		nparamv[paramc] = i * chunk;
		nparamv[paramc + 1] = (i + 1) * chunk;
		ocrEdtCreate(&pipelineGuid, pipelineTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence(dataGuid, pipelineGuid, 0, DB_MODE_ITW);
	}
	return NULL_GUID;
}

//                   0  1        2           3      4      5       6                     7                     8
// u64 nparamv[6] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};
ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 iterations = paramv[2];
	ocrGuid_t pipeExecTemplateGuid = paramv[6];
	ocrGuid_t nextIterTemplateGuid = paramv[7];
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
	ocrGuid_t pipeExecGuid, pipeExecDone, nextIterGuid;

	ocrEdtCreate(&pipeExecGuid, pipeExecTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &pipeExecDone);

	// Setup next iteration
	if (paramv[0] < iterations) {
		paramv[0] += 1;
		ocrEdtCreate(&nextIterGuid, nextIterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
					 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence(dataGuid, nextIterGuid, 0, DB_MODE_ITW);
		ocrAddDependence(pipeExecDone, nextIterGuid, 1, DB_MODE_RO);
	} 

	// Dependencies for pipeline
	ocrAddDependence(dataGuid, pipeExecGuid, 0, DB_MODE_ITW);
	// PRINTF("FINISHED ITER\n");
	return NULL_GUID;
}

//                   0        1           2      3       4       5
// u64 rparamv[6] = {db_size, iterations, split, verify, scalar, verbose};
ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	u64 db_size = paramv[0];
	u64 iterations = paramv[1];
	u64 split = paramv[2];
	int verify = (int) paramv[3];
	int verbose = (int) paramv[5];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE totalsum, timingsum[iterations], avg;

	for (i = 0; i < split; i++)
		for (j = 0; j < iterations; j++)
			timingsum[j] += data[3 * db_size + j];

	PRINTF("Timing Results:\n");
	for (i = 0; i < iterations; i++) {
		if (verbose)
			PRINTF("TRIAL %d: %f s\n", i + 1, timingsum[i]);
		totalsum += timingsum[i];
	}
	avg = totalsum / iterations;
	if (verbose)
		PRINTF("AVERAGE Time Per chunk: %f s\n", avg / split);
	PRINTF("AVERAGE Time Per Trial: %f s\n", avg);

	if (strcmp((char *) depv[1].ptr, "") != 0) 
		export_csv((char *) depv[1].ptr, iterations, timingsum, avg);

	if (verify) {
		STREAM_TYPE a = 0, ai, b = 0, bi, c = 0, ci;
		STREAM_TYPE scalar = (STREAM_TYPE) paramv[4];

		// Sum actual values
		for (i = 0; i < split; i++) {
			a += data[i];
			b += data[db_size + i];
			c += data[2 * db_size + i];
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
				   "Expected c: %f, Actual c: %f\n"
				   "Note: Expected values are multiplied by split (i.e. a * split)\n"
					"     Actual values are sum of split number of first elements", 
					split * ai, a, split * bi, b, split * ci, c);
	}
	ocrShutdown();
	return NULL_GUID;
}

ocrGuid_t mainEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 c, i, argc = getArgc(depv[0].ptr);
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

	u64 nparamc = 9, rparamc = 6, pparamc = 2;
	STREAM_TYPE * dataArray;
	char * efileArray;
	ocrGuid_t dataGuid, efileGuid, iterTemplateGuid, iterGuid, iterDone, resultsTemplateGuid, resultsGuid, 
			  pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, rparamc, 3);
	ocrEdtTemplateCreate(&pipeExecTemplateGuid, &pipeExecEdt, nparamc, 1);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, 2);
	ocrEdtTemplateCreate(&pipelineTemplateGuid, &pipelineEdt, nparamc + pparamc, 1);

	u64 nparamv[9] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};
	u64 rparamv[6] = {db_size, iterations, split, verify, scalar, verbose};

	// Formatting datablock
	DBCREATE(&dataGuid,(void **) &dataArray, sizeof(STREAM_TYPE) * (3 * db_size + iterations), 0, NULL_GUID, NO_ALLOC);
	for (i = 0; i < db_size; i++){
		dataArray[i] = 1.0;
		dataArray[db_size + i] = 2.0;
		dataArray[2 * db_size + i] = 0.0;
	}
	for (i = 3 * db_size; i < iterations; i++)
		dataArray[i] = 0.0;

	// Create datablock for export file name
	DBCREATE(&efileGuid, (void **) &efileArray, sizeof(char) * sizeof(efile), 0, NULL_GUID, NO_ALLOC);
	strcpy(efileArray, efile);

	// Create iterator and results EDTs
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &iterDone);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, rparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for iterator and results
	ocrAddDependence(dataGuid, resultsGuid, 0, DB_MODE_RO);
	ocrAddDependence(efileGuid, resultsGuid, 1, DB_MODE_RO);
	ocrAddDependence(iterDone, resultsGuid, 2, DB_MODE_RO);
	ocrAddDependence(dataGuid, iterGuid, 0, DB_MODE_ITW);
	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}