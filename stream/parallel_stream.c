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
	-- One large datablock is used with ranges of indexes corresponding to the following:
			  			 a --> [0 to (db_size - 1)]
			  			 b --> [db_size to (2 * db_size - 1)]
			 			 c --> [2 * db_size to (3 * db_size - 1)]
		  	  copy timing  --> [3 * db_size]
			  scale timing --> [3 * db_size + 1]
			  add timing   --> [3 * db_size + 2]
			  triad timing --> [3 * db_size + 3]
			  ...
			  ...
			  ...
			  copy timing  --> [3 * db_size + 4 * iterations]
			  scale timing --> [3 * db_size + 4 * iterations + 1]
			  add timing   --> [3 * db_size + 4 * iterations + 2]
			  triad timing --> [3 * db_size + 4 * iterations + 3]
*/

ocrGuid_t pipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[1];
	u64 begin = paramv[paramc - 2];
	u64 end = paramv[paramc - 1];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[5];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start, stop;

	// COPY
	start = mysecond();
	for (i = begin; i < end; i++)
		data[2 * db_size + i] = data[i];
	stop = mysecond();
	data[3 * db_size + 4 * (paramv[0] - 1)] += stop - start;

	// SCALE
	start = mysecond();
	for (i = begin; i < end; i++)
		data[db_size + i] = scalar * data[2 * db_size + i];
	stop = mysecond();
	data[3 * db_size + 4 * (paramv[0] - 1) + 1] += stop - start;

	// ADD
	start = mysecond();
	for (i = begin; i < end; i++)
		data[2 * db_size + i] = data[i] + data[db_size + i];
	stop = mysecond();
	data[3 * db_size + 4 * (paramv[0] - 1) + 2] += stop - start;

	// TRIAD
	start = mysecond();
	for (i = begin; i < end; i++)
		data[i] = data[db_size + i] + scalar * data[2 * db_size + i];
	stop = mysecond();
	data[3 * db_size + 4 * (paramv[0] - 1) + 3] += stop - start;

	//PRINTF("FINISHED PIPELINE\n");
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

	//PRINTF("FINISHED PIPEEXEC\n");
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

	//PRINTF("FINISHED ITER\n");
	return NULL_GUID;
}

//                   0        1           2      3       4       5
// u64 rparamv[6] = {db_size, iterations, split, verify, scalar, verbose};
// ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
// 	u64 i, j;
// 	u64 db_size = paramv[0];
// 	u64 iterations = paramv[1];
// 	u64 split = paramv[2];
// 	int verify = (int) paramv[3];
// 	STREAM_TYPE scalar = (STREAM_TYPE) paramv[4];
// 	int verbose = (int) paramv[5];
// 	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
// 	STREAM_TYPE totaltiming[4], timings[iterations][4], min[4], avg[4], max[4];

// 	double bytes[4] = {
//     	2 * sizeof(STREAM_TYPE) * db_size,
//     	2 * sizeof(STREAM_TYPE) * db_size,
//     	3 * sizeof(STREAM_TYPE) * db_size,
//     	3 * sizeof(STREAM_TYPE) * db_size
//    	};

// 	// Sum timings from each iteration
// 	for (i = 0; i < iterations; i++) {
// 		// Copy iteration timings to local data structure
// 		timings[i][0] += data[3 * db_size + 4 * i];
// 		timings[i][1] += data[3 * db_size + 4 * i + 1];
// 		timings[i][2] += data[3 * db_size + 4 * i + 2];
// 		timings[i][3] += data[3 * db_size + 4 * i + 3];

// 		// Print results from each iteration if verbose is specified
// 		if (verbose) {
// 			PRINTF(HLINE);
// 			PRINTF("ITERATION %d:\n", i + 1);
// 			PRINTF("Function       Rate MB/s     Time\n");
// 			for (j = 0; j < 4; j++)
// 				PRINTF("%s%12.1f %11.6f\n", label[j], 1.0E-06 * bytes[i] / timings[i][j], timings[i][j]);
// 		}

// 		// Add iteration to running sums
// 		totaltiming[0] += timings[i][0];
// 		totaltiming[1] += timings[i][1];
// 		totaltiming[2] += timings[i][2];
// 		totaltiming[3] += timings[i][3];

// 		// Set initial min and max values for each vector operation to first iteration
// 		if (i == 0) {
// 			for (j = 0; j < 4; j++) {
// 				min[j] = timings[i][j];
// 				max[j] = timings[i][j];
// 			}
// 		}

// 		// Compare current max/min to current iteration
// 		for (j = 0; j < 4; j++) {
// 			if (timings[i][j] > max[j])
// 				max[j] = timings[i][j];
// 			if (timings[i][j] < min[j])
// 				min[j] = timings[i][j];
// 		}
// 	}

// 	// Compute averages
// 	for (i = 0; i < 4; i++)
// 		avg[i] = totaltiming[i] / iterations;

// 	// Print overall results from iterations
// 	PRINTF(HLINE);
// 	PRINTF("OVERALL:\n");
// 	PRINTF("Function    Best Rate MB/s  Avg time     Min time     Max time\n");
// 	for (i = 0; i < 4; i++)
// 		PRINTF("%s%12.1f  %11.6f  %11.6f  %11.6f\n", label[i], 1.0E-06 * bytes[i] / avg[i],
// 													 avg[i], min[i], max[i]);
// 	PRINTF(HLINE);

// 	// Export to CSV
// 	if (strcmp((char *) depv[split].ptr, "") != 0) 
// 		export_csv((char *) depv[split].ptr, db_size, iterations, split, scalar, timings,
// 				   1.0E-06 * bytes[0] / avg[0], 1.0E-06 * bytes[1] / avg[1],
// 				   1.0E-06 * bytes[2] / avg[2], 1.0E-06 * bytes[3] / avg[3]);

// 	// Verify results
// 	if (verify) {
// 		STREAM_TYPE ai, bi, ci;
// 		STREAM_TYPE scalar = (STREAM_TYPE) paramv[4];
// 		int diff = 0;

// 		// Reproduce initializations
// 		ai = 1.0;
// 		bi = 2.0;
// 		ci = 0.0;

// 		// Execute timing loop
// 		for (i = 0; i < iterations; i++) {
// 			ci = ai;
// 			bi = scalar * ci;
// 			ci = ai + bi;
// 			ai = bi + scalar * ci;
// 		}

// 		// Compare against actual
// 		PRINTF("After %d Iterations:\n", iterations);
// 		for (i = 0; i < split; i++) {
// 			if (data[i] != ai) {
// 				diff += 1;
// 				PRINTF("Expected a: %f, Actual a: %f\n", ai, data[i]);
// 			}
// 			if (data[db_size + i] != bi) {
// 				diff += 1;
// 				PRINTF("Expected b: %f, Actual b: %f\n", bi, data[db_size + i]);
// 			}
// 			if (data[2 * db_size + i] != ci) {
// 				diff += 1;
// 				PRINTF("Expected c: %f, Actual c: %f\n", ci, data[2 * db_size + i]);
// 			}
// 		}
// 		PRINTF("%d differences between expected and actual\n", diff);
// 	}

// 	ocrShutdown();
// 	//PRINTF("FINISHED RESULTS\n");
// 	return NULL_GUID;
// }

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
	if (parseOptions(argc, argv,  &db_size, efile, &iterations, &split, &chunk, &verify, &scalar, &verbose)) {
		ocrShutdown();
		return NULL_GUID;
	}

	u64 nparamc = 9, pparamc = 2;
	STREAM_TYPE * dataArray;
	char * efileArray;
	ocrGuid_t dataGuid, efileGuid, iterTemplateGuid, iterGuid, iterDone, resultsTemplateGuid, resultsGuid, 
			  pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, RPARAMC, 3);
	ocrEdtTemplateCreate(&pipeExecTemplateGuid, &pipeExecEdt, nparamc, 1);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, 2);
	ocrEdtTemplateCreate(&pipelineTemplateGuid, &pipelineEdt, nparamc + pparamc, 1);

	u64 nparamv[9] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid};
	u64 rparamv[RPARAMC] = {db_size, iterations, split, chunk, verify, scalar, verbose, 1};

	// Formatting datablock
	DBCREATE(&dataGuid,(void **) &dataArray, sizeof(STREAM_TYPE) * (3 * db_size + 4 * iterations), 0, NULL_GUID, NO_ALLOC);
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

	//PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}
