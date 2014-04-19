#include "helper.h"
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
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
	One large datablock is used with ranges of indexes corresponding to the following:
			  a --> [0 to (db_size - 1)]
			  b --> [db_size to (2 * db_size - 1)]
			  c --> [2 * db_size to (3 * db_size - 1)]
		timings --> [3 * db_size to (3 * db_size + iterations - 1)]
*/

//                   0  1        2           3       4
// u64 nparamv[5] = {1, db_size, iterations, scalar, iterTemplateGuid};

ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[1];
	u64 iterations = paramv[2];
	STREAM_TYPE scalar = paramv[3];
	ocrGuid_t iterGuid;
	ocrGuid_t iterTemplateGuid = paramv[4];
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	data[3 * db_size + paramv[0] - 1] = mysecond();

	for (i = 0; i < db_size; i++)
		data[2 * db_size + i] = data[i];
	// PRINTF("FINISHED COPY\n");

	for (i = 0; i < db_size; i++)
		data[db_size + i] = scalar * data[2 * db_size + i];
	// PRINTF("FINISHED SCALE\n");

	for (i = 0; i < db_size; i++)
		data[2 * db_size + i] = data[i] + data[db_size + i];
	// PRINTF("FINISHED ADD\n");

	for (i = 0; i < db_size; i++)
		data[i] = data[db_size + i] + scalar * data[2 * db_size + i];
	// PRINTF("FINISHED TRIAD\n");

	data[3 * db_size + paramv[0] - 1] = mysecond() - data[3 * db_size + paramv[0] - 1];
	// PRINTF("FINISHED ITER\n");

	if (paramv[0] < iterations) {
		paramv[0] += 1;
		// Create next iterator EDT
		ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
					 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence(dataGuid, iterGuid, 0, DB_MODE_ITW);
	}
	return NULL_GUID;
}

//                   0        1           2       3       4
// u64 rparamv[5] = {db_size, iterations, verify, scalar, verbose};
ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[0];
	u64 iterations = paramv[1];
	int verify = (int) paramv[2];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[3];
	int verbose = (int) paramv[4];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE timing[iterations];
	STREAM_TYPE timingsum = 0.0;
	STREAM_TYPE avg = 0.0;

	PRINTF("Timing Results:\n");
	for (i = 0; i < iterations; i++) {
		timing[i] = data[3 * db_size + i];
		if (verbose)
			PRINTF("TRIAL %d: %f s\n", i + 1, timing[i]);
		timingsum += timing[i];
	}
	avg = timingsum / iterations;
	PRINTF("AVERAGE Time Per Trial: %f s\n", avg);

	if (strcmp((char *) depv[1].ptr, "") != 0) 
		export_csv((char *) depv[1].ptr, db_size, iterations, 0, scalar, timing, avg);

	if (verify) {
		STREAM_TYPE a = data[0];
		STREAM_TYPE b = data[db_size];
		STREAM_TYPE c = data[2 * db_size];
		STREAM_TYPE ai, bi, ci;
		STREAM_TYPE scalar = (STREAM_TYPE) paramv[3];

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
		if ((ai - a + bi - b + ci - c) == 0)
			PRINTF("No differences between expected and actual\n");
		else  
			PRINTF("Expected a: %f, Actual a: %f\n"
				   "Expected b: %f, Actual b: %f\n"
				   "Expected c: %f, Actual c: %f\n", ai, a, bi, b, ci, c);
	}

	ocrShutdown();
	return NULL_GUID;
}

ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 c, i, argc = getArgc(depv[0].ptr);
	char * argv[argc];
	for(i = 0; i < argc; i++)
		argv[i] = getArgv(depv[0].ptr, i);

	// Getopt arguments
	u64 db_size;
	char efile[100] = "";
	u64 iterations;
	int verify;
	STREAM_TYPE scalar;
	int verbose;

	// Parse getopt commands, shutdown and exit if help is selected
	if (parseOptions(argc, argv,  &db_size, efile, &iterations, &verify, &scalar, &verbose)) {
		ocrShutdown();
		return NULL_GUID;
	}

	// Initialize parameters, variables, and templates
	u64 nparamc = 5, rparamc = 5;
	STREAM_TYPE * dataArray;
	char * efileArray;
	ocrGuid_t dataGuid, efileGuid, iterTemplateGuid, iterGuid, iterDone, resultsTemplateGuid, resultsGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, rparamc, 3);

	u64 nparamv[5] = {1, db_size, iterations, scalar, iterTemplateGuid};
	u64 rparamv[5] = {db_size, iterations, verify, scalar, verbose};

	// Format datablock
	DBCREATE(&dataGuid, (void **) &dataArray, sizeof(STREAM_TYPE) * (3 * db_size + iterations), 0, NULL_GUID, NO_ALLOC);
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
