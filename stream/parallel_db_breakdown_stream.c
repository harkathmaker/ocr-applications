#include "helper.h"
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> 

#ifndef __OCR__
#	define __OCR__
#	include "ocr.h"
#endif

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

/*
	NOTES: 
	-- All results are in MB/s --> 1 MB = 10^6 B, NOT 2^20 B
	-- split datablocks are used with ranges of indexes corresponding to the following:
			  			 a --> [0 to (chunk - 1)]
			  			 b --> [chunk to (2 * chunk - 1)]
			  			 c --> [2 * chunk to (3 * chunk - 1)]
			  copy timing  --> [3 * chunk]
			  scale timing --> [3 * chunk + 1]
			  add timing   --> [3 * chunk + 2]
			  triad timing --> [3 * chunk + 3]
			  ...
			  ...
			  ...
			  copy timing  --> [3 * chunk + 4 * iterations]
			  scale timing --> [3 * chunk + 4 * iterations + 1]
			  add timing   --> [3 * chunk + 4 * iterations + 2]
			  triad timing --> [3 * chunk + 4 * iterations + 3]
*/

//                  0  1        2           3      4      5       6                     7                     8
// u64 paramv[9] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};
ocrGuid_t pipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[5];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start, stop;

	// COPY
	start = mysecond();
	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i];
	stop = mysecond();
	data[3 * chunk + 4 * (paramv[0] - 1)] = stop - start;

	// SCALE
	start = mysecond();
	for (i = 0; i < chunk; i++)
		data[chunk + i] = scalar * data[2 * chunk + i];
	stop = mysecond();
	data[3 * chunk + 4 * (paramv[0] - 1) + 1] = stop - start;

	// ADD
	start = mysecond();
	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i] + data[chunk + i];
	stop = mysecond();
	data[3 * chunk + 4 * (paramv[0] - 1) + 2] = stop - start;

	// TRIAD
	start = mysecond();
	for (i = 0; i < chunk; i++)
		data[i] = data[chunk + i] + scalar * data[2 * chunk + i];
	stop = mysecond();
	data[3 * chunk + 4 * (paramv[0] - 1) + 3] = stop - start;

	//PRINTF("FINISHED PIPELINE\n");
	return NULL_GUID;
}

//                  0  1        2           3      4      5       6                     7                     8
// u64 paramv[9] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};
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

	//PRINTF("FINISHED PIPE EXEC\n");
	return NULL_GUID;
}

//                  0  1        2           3      4      5       6                     7                     8
// u64 paramv[9] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};
ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 iterations = paramv[2];
	u64 split = paramv[3];
	ocrGuid_t pipeExecTemplateGuid = paramv[6];
	ocrGuid_t nextIterTemplateGuid = paramv[7];
	ocrGuid_t pipeExecGuid, pipeExecDone, nextIterGuid;

	ocrEdtCreate(&pipeExecGuid, pipeExecTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &pipeExecDone);

	// Setup next iteration
	if (paramv[0] < iterations) {
		paramv[0] += 1;
		ocrEdtCreate(&nextIterGuid, nextIterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
					 EDT_PROP_FINISH, NULL_GUID, NULL);
		for (i = 0; i < split; i++)
			ocrAddDependence((ocrGuid_t) depv[i].guid, nextIterGuid, i, DB_MODE_ITW);
		ocrAddDependence(pipeExecDone, nextIterGuid, split, DB_MODE_RO);
	} 

	// Dependencies for pipeline
	for (i = 0; i < split; i++)
		ocrAddDependence((ocrGuid_t) depv[i].guid, pipeExecGuid, i, DB_MODE_ITW);

	//PRINTF("FINISHED ITER\n");
	return NULL_GUID;
}

//                   0        1           2      3      4       5       6
// u64 rparamv[6] = {db_size, iterations, split, chunk, verify, scalar, verbose};
ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	u64 db_size = paramv[0];
	u64 iterations = paramv[1];
	u64 split = paramv[2];
	u64 chunk = paramv[3];
	int verify = (int) paramv[4];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[5];
	int verbose = (int) paramv[6];
	STREAM_TYPE totaltiming[4], timings[iterations][4], avg[4];

	// Bytes operated on for each vector operation
	// bytes[0] = copy, bytes[1] = scale, bytes[2] = add, bytes[3] = triad
	double bytes[4] = {
    	2 * sizeof(STREAM_TYPE) * db_size,
    	2 * sizeof(STREAM_TYPE) * db_size,
    	3 * sizeof(STREAM_TYPE) * db_size,
    	3 * sizeof(STREAM_TYPE) * db_size
   	};
	
	// Sum partial timings
	for (i = 0; i < split; i++) {
		STREAM_TYPE * cur = (STREAM_TYPE *) depv[i].ptr;
		for (j = 0; j < iterations; j++) {
			timings[j][0] += cur[3 * chunk + 4 * j];
			timings[j][1] += cur[3 * chunk + 4 * j + 1];
			timings[j][2] += cur[3 * chunk + 4 * j + 2];
			timings[j][3] += cur[3 * chunk + 4 * j + 3];
		}
	}

	// Set initial min and max values for each vector operation to first iteration
	STREAM_TYPE min[4] = {timings[0][0], timings[0][1], timings[0][2], timings[0][3]};
	STREAM_TYPE max[4] = {timings[0][0], timings[0][1], timings[0][2], timings[0][3]};

	// Sum timings from each iteration
	for (i = 0; i < iterations; i++) {
		// Print results from each iteration if verbose is specified
		if (verbose) {
			PRINTF(HLINE);
			PRINTF("ITERATION %d:\n", i + 1);
			PRINTF("Function       Rate MB/s     Time\n");
			for (j = 0; j < 4; j++)
				PRINTF("%s%12.1f %11.6f\n", label[j], 1.0E-06 * bytes[i] / timings[i][j], timings[i][j]);
		}
		totaltiming[0] += timings[i][0];
		totaltiming[1] += timings[i][1];
		totaltiming[2] += timings[i][2];
		totaltiming[3] += timings[i][3];
		for (j = 0; j < 4; j++) {
			if (timings[i][j] > max[j])
				max[j] = timings[i][j];
			if (timings[i][j] < min[j])
				min[j] = timings[i][j];
		}
	}

	// Compute averages
	for (i = 0; i < 4; i++)
		avg[i] = totaltiming[i] / iterations;

	// Print overall results from iterations
	PRINTF(HLINE);
	PRINTF("OVERALL:\n");
	PRINTF("Function    Best Rate MB/s  Avg time     Min time     Max time\n");
	for (i = 0; i < 4; i++)
		PRINTF("%s%12.1f  %11.6f  %11.6f  %11.6f\n", label[i], 1.0E-06 * bytes[i] / avg[i],
													 avg[i], min[i], max[i]);
	PRINTF(HLINE);

	// Export to CSV
	if (strcmp((char *) depv[split].ptr, "") != 0) 
		export_csv((char *) depv[split].ptr, db_size, iterations, split, scalar, timings,
				   1.0E-06 * bytes[0] / avg[0], 1.0E-06 * bytes[1] / avg[1],
				   1.0E-06 * bytes[2] / avg[2], 1.0E-06 * bytes[3] / avg[3]);

	// Verify results
	if (verify) {
		STREAM_TYPE ai, bi, ci;
		STREAM_TYPE scalar = (STREAM_TYPE) paramv[5];
		int diff = 0;

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

		// Compare against first elements of each "split" of actual
		PRINTF("After %d Iterations:\n", iterations);
		for (i = 0; i < split; i++) {
			STREAM_TYPE * cur = (STREAM_TYPE *) depv[i].ptr;
			if (cur[0] != ai) {
				diff += 1;
				PRINTF("Expected a: %f, Actual a: %f\n", ai, cur[0]);
			}
			if (cur[chunk] != bi) {
				diff += 1;
				PRINTF("Expected b: %f, Actual b: %f\n", bi, cur[chunk]);
			}
			if (cur[2 * chunk] != ci) {
				diff += 1;
				PRINTF("Expected c: %f, Actual c: %f\n", ci, cur[2 * chunk]);
			}
		}
		PRINTF("%d differences between expected and actual\n", diff);
	}

	ocrShutdown();
	//PRINTF("FINISHED RESULTS\n");
	return NULL_GUID;
}

ocrGuid_t mainEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
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

	u64 nparamc = 9, rparamc = 7;
	STREAM_TYPE * chunkArray;
	char * efileArray;
	ocrGuid_t dataGuids[split];
	ocrGuid_t iterTemplateGuid, iterGuid, iterDone, efileGuid, resultsTemplateGuid, resultsGuid, pipeExecTemplateGuid, 
			  nextIterTemplateGuid, pipelineTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, split);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, rparamc, split + 2);
	ocrEdtTemplateCreate(&pipeExecTemplateGuid, &pipeExecEdt, nparamc, split);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, split + 1);
	ocrEdtTemplateCreate(&pipelineTemplateGuid, &pipelineEdt, nparamc, 1);

	u64 nparamv[9] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};
	u64 rparamv[7] = {db_size, iterations, split, chunk, verify, scalar, verbose};

	// Preparing data block
	// 3 * db_size + 4 * iterations = (3 * chunk + 4 * iterations) * split
	for (i = 0; i < split; i++) {
		ocrGuid_t chunkGuid;
		DBCREATE(&chunkGuid,(void **) &chunkArray, sizeof(STREAM_TYPE)*(3 * chunk + 4 * iterations), 0, NULL_GUID, NO_ALLOC);
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

	//PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}
