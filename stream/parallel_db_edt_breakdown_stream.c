#include "helper.h"
#include "options.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h> 

#ifndef __OCR__
    #define __OCR__
    #include "ocr.h"
#endif

#ifndef STREAM_TYPE
	#define STREAM_TYPE double
#endif

#define PARALLEL_STREAM

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

//                  0  1        2           3      4      5       6                     7 
// u64 paramv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
//                   8                     9                 10                 11               12
//					 pipelineTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid}; 
ocrGuid_t copyEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	//START_PROFILE("copy");
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start = mysecond();
	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i];
	STREAM_TYPE stop = mysecond();
	data[3 * chunk + 4 * (paramv[0] - 1)] = stop - start;

	// PRINTF("FINISHED COPY\n");
	return NULL_GUID;
}

ocrGuid_t scaleEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE scalar = paramv[5];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start = mysecond();
	for (i = 0; i < chunk; i++)
		data[chunk + i] = scalar * data[2 * chunk + i];
	STREAM_TYPE stop = mysecond();
	data[3 * chunk + 4 * (paramv[0] - 1) + 1] = stop - start;

	// PRINTF("FINISHED SCALE\n");
	return NULL_GUID;
}

ocrGuid_t addEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start = mysecond();
	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i] + data[chunk + i];
	STREAM_TYPE stop = mysecond();
	data[3 * chunk + 4 * (paramv[0] - 1) + 2] = stop - start;
	// PRINTF("FINISHED ADD\n");
	return NULL_GUID;
}

ocrGuid_t triadEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[4];
	STREAM_TYPE scalar = paramv[5];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE start = mysecond();
	for (i = 0; i < chunk; i++)
		data[i] = data[chunk + i] + scalar * data[2 * chunk + i];
	STREAM_TYPE stop = mysecond();
	data[3 * chunk + 4 * (paramv[0] - 1) + 3] = stop - start;

	// PRINTF("FINISHED TRIAD\n");
	return NULL_GUID;
}

//                  0  1        2           3      4      5       6                     7 
// u64 paramv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
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
	ocrAddDependence(dataGuid,  triadGuid, 0, DB_MODE_ITW);
	ocrAddDependence(addDone,   triadGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid,  addGuid,   0, DB_MODE_ITW);
	ocrAddDependence(scaleDone, addGuid,   1, DB_MODE_RO);
	ocrAddDependence(dataGuid,  scaleGuid, 0, DB_MODE_ITW);
	ocrAddDependence(copyDone,  scaleGuid, 1, DB_MODE_RO);
	ocrAddDependence(dataGuid,  copyGuid,  0, DB_MODE_ITW);

	// PRINTF("FINISHED PIPELINE\n");
	return NULL_GUID;
}

//                  0  1        2           3      4      5       6                     7 
// u64 paramv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
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

//                  0  1        2           3      4      5       6                     7 
// u64 paramv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
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
	if (parseOptions(argc, argv,  &db_size, efile, &iterations, &split, &chunk, &verify, &scalar, &verbose)) {
		ocrShutdown();
		return NULL_GUID;
	}

	u64 nparamc = 13;
	STREAM_TYPE * chunkArray;
	char * efileArray;
	ocrGuid_t dataGuids[split];
	ocrGuid_t iterTemplateGuid, iterGuid, iterDone, efileGuid, resultsTemplateGuid, resultsGuid, pipeExecTemplateGuid, pipelineTemplateGuid, 
			  copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid, nextIterTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, split);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, RPARAMC, split + 2);

	ocrEdtTemplateCreate(&pipeExecTemplateGuid, &pipeExecEdt, nparamc, split);
	ocrEdtTemplateCreate(&nextIterTemplateGuid, &iterEdt, nparamc, split + 1);
	ocrEdtTemplateCreate(&pipelineTemplateGuid, &pipelineEdt, nparamc, 1);

	ocrEdtTemplateCreate(&copyTemplateGuid, &copyEdt, nparamc, 1);
	ocrEdtTemplateCreate(&scaleTemplateGuid, &scaleEdt, nparamc, 2);
	ocrEdtTemplateCreate(&addTemplateGuid, &addEdt, nparamc, 2);
	ocrEdtTemplateCreate(&triadTemplateGuid, &triadEdt, nparamc, 2);

	u64 nparamv[13]= {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, 
					 pipelineTemplateGuid, copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid}; 
	u64 rparamv[RPARAMC] = {db_size, iterations, split, chunk, verify, scalar, verbose, 0};

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

	// PRINTF("FINISHED MAIN\n");
	return NULL_GUID;
}
