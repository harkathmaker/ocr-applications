#include <string.h>
#include <sys/time.h> 

#ifndef __OCR__
#	define __OCR__
#	include "ocr.h"
#endif

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

#ifndef RPARAMC
#	define RPARAMC 8
#endif

#define HLINE "---------------------------------------------------------------\n"

char * label[4] = {"Copy:      ", "Scale:     ", "Add:       ", "Triad:     "};

struct args {
	u64 db_size;
	char efile[100];
	u64 iterations;
	u64 split;
	u64 chunk;
	int verify;
	STREAM_TYPE scalar;
	int verbose;
};

struct templates {
	ocrGuid_t nextIterTemplateGuid;
	ocrGuid_t pipeExecTemplateGuid;
	ocrGuid_t finePipelineTemplateGuid;
	ocrGuid_t coarsePipelineTemplateGuid;
	ocrGuid_t copyTemplateGuid;
	ocrGuid_t scaleTemplateGuid;
	ocrGuid_t addTemplateGuid;
	ocrGuid_t triadTemplateGuid;
};

int init_args(u64 argc, char ** argv, struct args * a) {
	if (parseOptions(argc, argv, a)) {
		ocrShutdown();
		return 0;
	}
	return 1;
}

void export_csv(char * name, u64 db_size, u64 iterations, u64 split, STREAM_TYPE scalar, STREAM_TYPE trials[][4],
				STREAM_TYPE cavg, STREAM_TYPE savg, STREAM_TYPE aavg, STREAM_TYPE tavg) {
	char path[150];
	//strcpy(path, "./results/");
	//strcat(path, name);
	strcpy(path, name);
	FILE * f = fopen(path, "a");
	if (f == NULL) {
		PRINTF("Error creating export file.");
		exit(1);
	}

	u64 i;
		//for (i = 0; i < iterations; i++) 
			//fprintf(f, "%llu %f, ", db_size, trials[i]);
	fprintf(f, "%llu, %llu, %llu, %.2f, %f, %f, %f, %f\n", db_size, iterations, split, scalar, cavg, savg, aavg, tavg); 

	fclose(f);
	return;
}

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

ocrGuid_t copyEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
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

ocrGuid_t finePipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
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

	PRINTF("FINISHED PIPELINE SETUP\n");
	return NULL_GUID;
}

ocrGuid_t coarsePipelineEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
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

ocrGuid_t pipeExecEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 split = paramv[3];
	u64 chunk = paramv[4];
	ocrGuid_t pipelineTemplateGuid = paramv[8];
	u64 parallel = paramv[paramc - 3];
	ocrGuid_t pipelineGuid;
	ocrGuid_t dataGuid = (ocrGuid_t) depv[0].guid;

	// Spawn pipeline children operating on chunk amounts of data
	for (i = 0; i < split; i++) {
		if (parallel == 1) {
			paramv[paramc - 2] = i * chunk;
			paramv[paramc - 1] = (i + 1) * chunk;
		}
		ocrEdtCreate(&pipelineGuid, pipelineTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, NULL);
		if (parallel == 1)
			ocrAddDependence(dataGuid, pipelineGuid, 0, DB_MODE_ITW);
		else 
			ocrAddDependence((ocrGuid_t) depv[i].guid, pipelineGuid, 0, DB_MODE_ITW);
	}

	PRINTF("FINISHED PIPE EXEC\n");
	return NULL_GUID;
}

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

	PRINTF("FINISHED ITER\n");
	return NULL_GUID;
}

//                                 0 	         1               2         3         4                5            6               7
// u64 rparamv[8] = {db_size, iterations, verify, scalar, verbose, <split>, <chunk>, <parallel>}
ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j;
	u64 db_size = paramv[0];
	u64 iterations = paramv[1];
	u64 split = paramv[2];
	u64 chunk = paramv[3];
	int verify = (int) paramv[4];
	STREAM_TYPE scalar = (STREAM_TYPE) paramv[5];
	int verbose = (int) paramv[6];
	int parallel = (int) paramv[7];
	STREAM_TYPE totaltiming[4], timings[iterations][4], avg[4];

	u64 actual_split = 0;
	u64 actual_chunk = 0;
	if (parallel == 1) {
		actual_split = split;
		split = 1;
		actual_chunk = chunk;
		chunk = db_size;
	}

	// Initially set values to zero
	memset(totaltiming, 0, 4 * sizeof(totaltiming[0]));
	memset(timings, 0, iterations * 4 * sizeof(timings[0][0]));
	memset(avg, 0, 4 * sizeof(avg[0]));

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
	if (parallel == 1) {
		split = actual_split;
		chunk = actual_chunk;
	}
	if (strcmp((char *) depv[depc - 2].ptr, "") != 0) 
		export_csv((char *) depv[split].ptr, db_size, iterations, split, scalar, timings, 
				   1.0E-06 * bytes[0] / avg[0], 1.0E-06 * bytes[1] / avg[1],
				   1.0E-06 * bytes[2] / avg[2], 1.0E-06 * bytes[3] / avg[3]);

	// Verify results
	if (verify) {
		STREAM_TYPE ai, bi, ci;
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

		// Compare against actual
		PRINTF("After %d Iterations:\n", iterations);
		STREAM_TYPE * cur = (STREAM_TYPE *) depv[0].ptr;
		for (i = 0; i < split; i++) {
			if (parallel != 1) {
				STREAM_TYPE * cur = (STREAM_TYPE *) depv[i].ptr;
			}
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
	return NULL_GUID;
}


void create_templates(struct templates * t, u64 split, u64 nparamc){//(u64 argc, char * argv) {

	ocrEdtTemplateCreate(&t->pipeExecTemplateGuid, &pipeExecEdt, nparamc, split);
	ocrEdtTemplateCreate(&t->nextIterTemplateGuid, &iterEdt, nparamc, split + 1);
	ocrEdtTemplateCreate(&t->finePipelineTemplateGuid, &finePipelineEdt, nparamc, 1);
	ocrEdtTemplateCreate(&t->coarsePipelineTemplateGuid, &coarsePipelineEdt, nparamc, 1);

	ocrEdtTemplateCreate(&t->copyTemplateGuid, &copyEdt, nparamc, 1);
	ocrEdtTemplateCreate(&t->scaleTemplateGuid, &scaleEdt, nparamc, 2);
	ocrEdtTemplateCreate(&t->addTemplateGuid, &addEdt, nparamc, 2);
	ocrEdtTemplateCreate(&t->triadTemplateGuid, &triadEdt, nparamc, 2);
	return;
}

// void init_single_datablock(ocrGuid_t * dataGuid) {
// 	u64 i;
// 	STREAM_TYPE * dataArray;
// 	char * efileArray;
// 	DBCREATE(dataGuid,(void **) &dataArray, sizeof(STREAM_TYPE) * (3 * db_size + 4 * iterations), 0, NULL_GUID, NO_ALLOC);
// 	for (i = 0; i < db_size; i++){
// 		dataArray[i] = 1.0;
// 		dataArray[db_size + i] = 2.0;
// 		dataArray[2 * db_size + i] = 0.0;
// 	}
// 	for (i = 3 * db_size; i < iterations; i++)
// 		dataArray[i] = 0.0;

// 	DBCREATE(&efileGuid, (void **) &efileArray, sizeof(char) * sizeof(efile), 0, NULL_GUID, NO_ALLOC);
// 	strcpy(efileArray, efile);
// }

void create_multi_dbs(u64 split, u64 chunk, u64 iterations, ocrGuid_t *dataGuids) {
	//printf("split = %d\n", split);
	//printf("db_size = %d\n", db_size);
	// Preparing data block
	// 3 * db_size + 4 * iterations = (3 * chunk + 4 * iterations) * split
	u64 i, j;
	STREAM_TYPE * chunkArray;
	for (i = 0; i < split; i++) {
		ocrGuid_t chunkGuid; 
		DBCREATE(&chunkGuid,(void **) &chunkArray, 
		sizeof(STREAM_TYPE)*(3 * chunk + 4 * iterations), 0, NULL_GUID, NO_ALLOC);
		for (j = 0; j < chunk; j++) {
			chunkArray[j] = 1.0;
			chunkArray[chunk + j] = 2.0;
			chunkArray[2 * chunk + j] = 0.0;
		}
		dataGuids[i] = chunkGuid;
	}
	return;
}



// void init_single_edt(ocrGuid_t * dataGuid) {
// 	u64 nparamv[12] = {1, db_size, iterations, split, chunk, scalar, pipeExecTemplateGuid, nextIterTemplateGuid, pipelineTemplateGuid, 1, 0, chunk};
// 	u64 rparamv[RPARAMC] = {db_size, iterations, split, chunk, verify, scalar, verbose, 1};
// 	//u64 rparamv[RPARAMC] = {db_size, iterations, split, chunk, verify, scalar, verbose, 1};
// 	// Create iterator and results EDTs
// 	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
// 				 EDT_PROP_FINISH, NULL_GUID, &iterDone);
// 	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, rparamv, EDT_PARAM_DEF, NULL_GUID,
// 				 EDT_PROP_NONE, NULL_GUID, NULL);

// 	// Dependencies for iterator and results
// 	ocrAddDependence(*dataGuid, resultsGuid, 0, DB_MODE_RO);
// 	ocrAddDependence(efileGuid, resultsGuid, 1, DB_MODE_RO);
// 	ocrAddDependence(iterDone, resultsGuid, 2, DB_MODE_RO);
// 	ocrAddDependence(*dataGuid, iterGuid, 0, DB_MODE_ITW);
// }

void create_multi_db_edts(u64 nparamc, u64 * nparamv, ocrGuid_t * dataGuids, struct args * a) {
	u64 i;

	u64 rparamv[RPARAMC] = {a->db_size, a->iterations, a->split, a->chunk, a->verify, a->scalar, a->verbose, 0};

	// // Create datablock for export file name
	ocrGuid_t efileGuid;
	char * efileArray;
	DBCREATE(&efileGuid, (void **) &efileArray, sizeof(char) * sizeof(a->efile), 0, NULL_GUID, NO_ALLOC);
	strcpy(efileArray, a->efile);

	// Create iter and results templates
	ocrGuid_t iterGuid, iterDone, iterTemplateGuid, resultsGuid, resultsTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, nparamc, a->split);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, RPARAMC, a->split + 2);

	// Create iterator and results EDTs
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_FINISH, NULL_GUID, &iterDone);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, rparamv, EDT_PARAM_DEF, NULL_GUID,
				 EDT_PROP_NONE, NULL_GUID, NULL);

	// Dependencies for iterator and results
	for(i = 0; i < a->split; i++)
		ocrAddDependence(dataGuids[i], resultsGuid, i, DB_MODE_RO);
	ocrAddDependence(efileGuid, resultsGuid, a->split, DB_MODE_RO);
	ocrAddDependence(iterDone, resultsGuid, a->split + 1, DB_MODE_RO);
	for (i = 0; i < a->split; i++)
		ocrAddDependence(dataGuids[i], iterGuid, i, DB_MODE_ITW);
}

// printf("db_size = %llu\n", a->db_size);