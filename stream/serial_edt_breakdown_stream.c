#include "helper.h"
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

//                 0  1        2           3       4                 5                  6                7                  8
// u64 paramv[9]= {1, db_size, iterations, scalar, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid};
ocrGuid_t copyEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[1];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start = mysecond();
	for (i = 0; i < db_size; i++)
		data[2 * db_size + i] = data[i];
	STREAM_TYPE stop = mysecond();
	data[3 * db_size + 4 * (paramv[0] - 1)] = stop - start;
	
	//PRINTF("FINISHED COPY\n");
	return NULL_GUID;
}

ocrGuid_t scaleEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[1];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[3];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start = mysecond();
	for (i = 0; i < db_size; i++)
		data[db_size + i] = scalar * data[2 * db_size + i];
	STREAM_TYPE stop = mysecond();
	data[3 * db_size + 4 * (paramv[0] - 1) + 1] = stop - start;
	
	//PRINTF("FINISHED SCALE\n");
	return NULL_GUID;
}

ocrGuid_t addEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[1];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start = mysecond();
	for (i = 0; i < db_size; i++)
		data[2 * db_size + i] = data[i] + data[db_size + i];
	STREAM_TYPE stop = mysecond();
	data[3 * db_size + 4 * (paramv[0] - 1) + 2] = stop - start;
	
	//PRINTF("FINISHED ADD\n");
	return NULL_GUID;
}

ocrGuid_t triadEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 db_size = paramv[1];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[3];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start = mysecond();
	for (i = 0; i < db_size; i++)
		data[i] = data[db_size + i] + scalar * data[2 * db_size + i];
	STREAM_TYPE stop = mysecond();
	data[3 * db_size + 4 * (paramv[0] - 1) + 3] = stop - start;
	
	//PRINTF("FINISHED TRIAD\n");
	return NULL_GUID;
}

//                 0  1        2           3       4                 5                  6                7                  8
// u64 paramv[9]= {1, db_size, iterations, scalar, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid}; 
ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 iterations = paramv[2];
	ocrGuid_t copyTemplateGuid = paramv[4];
	ocrGuid_t scaleTemplateGuid = paramv[5];
	ocrGuid_t addTemplateGuid = paramv[6];
	ocrGuid_t triadTemplateGuid = paramv[7];
	ocrGuid_t nextIterTemplateGuid = paramv[8];
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;

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

	if (paramv[0] < iterations) {
		paramv[0] += 1;
		// EDT for next iteration
		ocrEdtCreate(&nextIterGuid, nextIterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
					 EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence(dataGuid, nextIterGuid, 0, DB_MODE_ITW);
		ocrAddDependence(triadDone, nextIterGuid, 1, DB_MODE_RO);
	}

	// Dependencies for vector operations
	ocrAddDependence(dataGuid, triadGuid, 0, DB_MODE_ITW);
	ocrAddDependence(addDone, triadGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, addGuid, 0, DB_MODE_ITW);
	ocrAddDependence(scaleDone, addGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, scaleGuid, 0, DB_MODE_ITW);
	ocrAddDependence(copyDone, scaleGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid, copyGuid, 0, DB_MODE_ITW);
	
	//PRINTF("FINISHED ITER\n");
	return NULL_GUID;
}
//                   0        1           2       3       4
// u64 paramv[5] = {db_size, iterations, verify, scalar, verbose};
ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	u64 db_size = paramv[0];
	u64 iterations = paramv[1];
	int verify = (int) paramv[2];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[3];
	int verbose = (int) paramv[4];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE totaltiming[4], timings[iterations][4], min[4], avg[4], max[4];

	// Bytes operated on for each vector operation
	// bytes[0] = copy, bytes[1] = scale, bytes[2] = add, bytes[3] = triad
	double bytes[4] = {
    	2 * sizeof(STREAM_TYPE) * db_size,
    	2 * sizeof(STREAM_TYPE) * db_size,
    	3 * sizeof(STREAM_TYPE) * db_size,
    	3 * sizeof(STREAM_TYPE) * db_size
   	};

	// Sum timings from each iteration
	for (i = 0; i < iterations; i++) {
		// Copy iterations to local data structure
		timings[i][0] += data[3 * db_size + 4 * i];
		timings[i][1] += data[3 * db_size + 4 * i + 1];
		timings[i][2] += data[3 * db_size + 4 * i + 2];
		timings[i][3] += data[3 * db_size + 4 * i + 3];

		// Print results from each iteration if verbose is specified
		if (verbose) {
			PRINTF(HLINE);
			PRINTF("ITERATION %d:\n", i + 1);
			PRINTF("Function       Rate MB/s     Time\n");
			for (j = 0; j < 4; j++)
				PRINTF("%s%12.1f %11.6f\n", label[j], 1.0E-06 * bytes[i] / timings[i][j], timings[i][j]);
		}

		// Add iteration to running sums
		totaltiming[0] += timings[i][0];
		totaltiming[1] += timings[i][1];
		totaltiming[2] += timings[i][2];
		totaltiming[3] += timings[i][3];

		// Set initial min and max values for each vector operation to first iteration
		if (i == 0) {
			for (j = 0; j < 4; j++) {
				min[j] = timings[i][j];
				max[j] = timings[i][j];
			}
		}

		// Compare current max/min to current iteration
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
	if (strcmp((char *) depv[1].ptr, "") != 0) 
		export_csv((char *) depv[1].ptr, db_size, iterations, 0, scalar, timings,
				   1.0E-06 * bytes[0] / avg[0], 1.0E-06 * bytes[1] / avg[1],
				   1.0E-06 * bytes[2] / avg[2], 1.0E-06 * bytes[3] / avg[3]);

	// Verify results
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
	if (parseOptions(argc, argv,  &db_size, efile, &iterations, NULL, NULL, &verify, &scalar, &verbose)) {
		ocrShutdown();
		return NULL_GUID;
	}

	u64 nparamc = 9, rparamc = 5;
	STREAM_TYPE * dataArray;
	char * efileArray;
	ocrGuid_t dataGuid, efileGuid, iterGuid, iterDone, iterTemplateGuid, resultsGuid, resultsTemplateGuid, copyTemplateGuid, 
			  scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, rparamc, 3);
	ocrEdtTemplateCreate(&copyTemplateGuid, &copyEdt, nparamc, 1);
	ocrEdtTemplateCreate(&scaleTemplateGuid, &scaleEdt, nparamc, 2);
	ocrEdtTemplateCreate(&addTemplateGuid, &addEdt, nparamc, 2);
	ocrEdtTemplateCreate(&triadTemplateGuid, &triadEdt, nparamc, 2);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, 2);

	u64 nparamv[9]= {1, db_size, iterations, scalar, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid}; 
	u64 rparamv[5] = {db_size, iterations, verify, scalar, verbose};

	// Preparing datablock
	DBCREATE(&dataGuid,(void **) &dataArray, sizeof(STREAM_TYPE)*(3 * db_size + 4 * iterations), 0, NULL_GUID, NO_ALLOC);
	for (i = 0; i < db_size; i++){
		dataArray[i] = 1.0;
		dataArray[db_size + i] = 2.0;
		dataArray[2 * db_size + i] = 0.0;
	}
	for (i = 3 * db_size; i < 4 * iterations; i++)
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
