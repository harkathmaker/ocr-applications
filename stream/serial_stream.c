#define __OCR__
#include "ocr.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h> 
#include <unistd.h>

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


// scalar
// db_size
//	u64 nparamv[2=5] = {1, db_size, scalar, iterations, iterTemplateGuid};
ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[1];
	u64 scalar = paramv[2];
	u64 iterations = paramv[3];
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
		ocrAddDependence(dataGuid, iterGuid, 0, DB_MODE_EW);
	}
	return NULL_GUID;
}

ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 verify = paramv[0];
	u64 verbose = paramv[1];
	u64 db_size = paramv[2];
	u64 iterations = paramv[3];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE timing;
	STREAM_TYPE timingsum = 0.0;

	PRINTF("Timing Results:\n");
	for (i = 0; i < iterations; i++) {
		timing = data[3 * db_size + i];
		PRINTF("TRIAL %d: %f s\n", i + 1, timing);
		timingsum += timing;
	}
	PRINTF("AVERAGE: %f s\n", timingsum / iterations);

	if (verify) {
		STREAM_TYPE a = data[0];
		STREAM_TYPE b = data[db_size];
		STREAM_TYPE c = data[2 * db_size];
		STREAM_TYPE ai, bi, ci;
		STREAM_TYPE scalar = paramv[4];

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
	char *argv[argc];
	for(i = 0; i < argc; i++)
		argv[i] = getArgv(depv[0].ptr,i);

	u64 db_size = 100;
	STREAM_TYPE scalar = 3.0;
	u64 iterations = 1;
	char efile[100];
	static int verify = 0;
	static int verbose = 0;

	while (1) {
		int option_index = 0;
		static struct option long_options[] = {
			{"export_to",   required_argument,  0,              'e'},
			{"help",        no_argument,        0,              'h'},
			{"iterations",  required_argument,  0,              'i'},
			{"db_size",     required_argument,  0,              'n'},
			{"scalar",      required_argument,  0,              's'},
			{"verify",      no_argument,        &verify,   1},
			{"verbose",     no_argument,        &verbose,  1},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, "e:hi:n:rs:v", long_options, &option_index);
		if (c == -1)
			break;
		switch (c) {
			case 'e':
				sscanf(optarg, "%s", &efile);
				break;
			case 'h':
				printf("serial_stream [-e file_name] [-h] [-i num_iter] [-n db_size] "
					   "[-r] [-s scalar_value] [-v]\n"
					   "List of Options\n" 
					   "-e = exports results to csv file\n"
					   "-h = list of options\n"
					   "-i = number of iterations (1 by default)\n"
					   "-n = size of data blocks a, b, c (10 million by default)\n"
					   "-r = verify results\n"
					   "-s = value of scalar (3.0 by default)\n"
					   "-v = verbose output\n");
					   ocrShutdown();
					   return NULL_GUID;
			case 'i':
				sscanf(optarg, "%d", &iterations);
				break;
			case 'n':
				sscanf(optarg, "%llu", &db_size);
				break;
			case 'r':
				verify = 1;
				break;
			case 's':
				sscanf(optarg, "%llu", &scalar);
				break;
			case 'v':
				verbose = 1;
				break;
			case '?':
				PRINTF("Unknown option -%c.\n",optopt);
		}
	}

	u64 nparamc = 5, rparamc = 5;
	STREAM_TYPE * dataArray;
	ocrGuid_t dataGuid, iterTemplateGuid, iterGuid, iterDone, resultsTemplateGuid, resultsGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, rparamc, 2);

	u64 nparamv[5] = {1, db_size, scalar, iterations, iterTemplateGuid};
	u64 rparamv[5] = {verify, verbose, db_size, iterations, scalar};

	// Formatting datablock
	DBCREATE(&dataGuid, (void **) &dataArray, sizeof(STREAM_TYPE) * (3 * db_size + iterations), 0, NULL_GUID, NO_ALLOC);
	for (i = 0; i < db_size; i++){
		dataArray[i] = 1.0;
		dataArray[db_size + i] = 2.0;
		dataArray[2 * db_size + i] = 0.0;
	}
	for (i = 3 * db_size; i < iterations; i++)
		dataArray[i] = 0.0;

	// Create iterator and results EDTs
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &iterDone);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, rparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for iterator and results
	ocrAddDependence(dataGuid, resultsGuid, 0, DB_MODE_RO);
	ocrAddDependence(iterDone, resultsGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, iterGuid, 0, DB_MODE_EW);
	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}