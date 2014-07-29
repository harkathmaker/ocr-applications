#include <string.h>
#include <sys/time.h> 

#ifndef __OCR__
	#define __OCR__
	#include "ocr.h"
#endif

#ifndef STREAM_TYPE
	#define STREAM_TYPE double
#endif

#ifndef RPARAMC
	#define RPARAMC 8
#endif

#define NUM_VECOP 4

#define COPY_IDX    0
#define SCALE_IDX  1
#define ADD_IDX     2
#define TRIAD_IDX   3

#define PARAMC 1

//#define DEBUG

#define HLINE "---------------------------------------------------------------\n"

char * label[4] = {"Copy:      ", "Scale:     ", "Add:       ", "Triad:     "};

// ----------------------------------
// GLOBAL STRUCTS
// ----------------------------------
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
	ocrGuid_t streamTemplateGuid;
	//ocrGuid_t serialStreamTemplateGuid;
};

struct params{
	u64 cur_itr;
	struct templates tmplt;
	struct args args;
};// params;

// ----------------------------------
// TIMERS
// ----------------------------------

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void export_csv(char * name, u64 db_size, u64 iterations, u64 split, STREAM_TYPE scalar, STREAM_TYPE trials[][4],
	STREAM_TYPE cavg, STREAM_TYPE savg, STREAM_TYPE aavg, STREAM_TYPE tavg) {
	u64 i;
	FILE * f = fopen(name, "a");

	if (f == NULL) {
		PRINTF("Error creating export file.");
		exit(1);
	}

	fprintf(f, "%llu, %llu, %llu, %.2f, %f, %f, %f, %f\n", db_size, iterations, split, scalar, cavg, savg, aavg, tavg); 

	fclose(f);
	return;
}

ocrGuid_t endTimerEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	STREAM_TYPE stop = mysecond();
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	PRINTF("INSIDE END TIMER EDT\n");
	int time_idx = NUM_VECOP * (paramv[0] - 1) + paramv[1];
	data[time_idx] = stop - data[time_idx];
	//printf("paramv[0] =  %llu\n", paramv[0]);
	//printf("time_idx = %d\n", time_idx);
	//PRINTF("time=%f\n", data[time_idx]);

	return NULL_GUID;
}

// ----------------------------------
// VECTOR OPERATION EDTS
// ----------------------------------
ocrGuid_t copyEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[0];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i];

	PRINTF("FINISHED COPY\n");
	return NULL_GUID;
}

ocrGuid_t scaleEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[0];
	STREAM_TYPE scalar = paramv[1];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	for (i = 0; i < chunk; i++)
		data[chunk + i] = scalar * data[2 * chunk + i];

	PRINTF("FINISHED SCALE\n");
	return NULL_GUID;
}

ocrGuid_t addEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[0];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i] + data[chunk + i];

	PRINTF("FINISHED ADD\n");
	return NULL_GUID;
}

ocrGuid_t triadEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[0];
	STREAM_TYPE scalar = paramv[1];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	for (i = 0; i < chunk; i++)
		data[i] = data[chunk + i] + scalar * data[2 * chunk + i];

	PRINTF("FINISHED TRIAD\n");
	return NULL_GUID;
}

// ----------------------------------
// SETUP + RESULTS EDTS
// ----------------------------------

ocrGuid_t runVecOpEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	//printf("iterations = %llu\n", params->args.iterations);
	u64 i;
	u64 cur_itr= paramv[0];
	//printf("run_vecop cur_itr = %llu\n", cur_itr);
	struct params * params = (struct params *) depv[depc - 1].ptr;
	u64 split = params->args.split;
	u64 chunk = params->args.chunk;
	u64 scalar = params->args.scalar;
	u64 vecOp = paramv[1];
	ocrGuid_t endTimerGuid, endTimerTemplateGuid;
	ocrGuid_t vecOpGuid, vecOpDone;
	ocrGuid_t allVecOpDone, allVecOpReady;
	ocrGuid_t vecOpTemplateGuid = (ocrGuid_t) paramv[2];

	// Create End Timer EDT
	u64 tparamc = 2;
	u64 tparamv[2] ={cur_itr, vecOp};

	ocrEdtTemplateCreate(&endTimerTemplateGuid, &endTimerEdt, tparamc, split + 1);
	ocrEdtCreate(&endTimerGuid, endTimerTemplateGuid, EDT_PARAM_DEF, tparamv, EDT_PARAM_DEF, NULL_GUID,
		EDT_PROP_FINISH, NULL_GUID, NULL);
	ocrAddDependence((ocrGuid_t) depv[split].guid, endTimerGuid, 0, DB_MODE_ITW);

	// Create "Start" Event signalling all instances of the vector operation to run
	ocrEventCreate(&allVecOpReady, OCR_EVENT_STICKY_T, true);

	// Create "split" number of selected vector operation EDTs
	u64 nparamv[2] = {chunk, scalar};

	for (i = 0; i < split; i++) {
		ocrEdtCreate(&vecOpGuid, vecOpTemplateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF, NULL_GUID,
			EDT_PROP_FINISH, NULL_GUID, &vecOpDone);
		// Run the vector operation EDTs once all of the vector operation EDTs have been created
		ocrAddDependence((ocrGuid_t) depv[i].guid, vecOpGuid, 0, DB_MODE_ITW);
		ocrAddDependence(allVecOpReady, vecOpGuid, 1, DB_MODE_RO);
		ocrAddDependence(vecOpDone, endTimerGuid, i + 1, DB_MODE_RO);
	}

	// Start the timer and run "split" number of the vector operation EDTs
	PRINTF("FINISHED CREATING ALL VECTOR OPERATIONS\n");
	u64 timer_idx =  NUM_VECOP * (cur_itr- 1) + vecOp;
	//printf("cur_itr = %llu\n", cur_itr);
	//printf("vecop = %llu\n", vecOp);
	//printf("index = %llu\n", timer_idx);
	STREAM_TYPE * data = (STREAM_TYPE *) depv[split].ptr;
	data[timer_idx] = mysecond();

	// Start all selected vector operation EDTs
	ocrEventSatisfy(allVecOpReady, NULL_GUID);

	PRINTF("FINISHED RUN VECOP\n");
	return NULL_GUID;
}

ocrGuid_t parallelStreamEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	// Create local variables of parameters

	//u64 chunk = params->args.chunk;
	//STREAM_TYPE scalar = params->args.scalar;

	struct params * params = (struct params *) depv[0].ptr;
	u64 split = params->args.split;
	u64 cur_itr = params->cur_itr;
	//printf("parallel_cur_itr = %llu\n", cur_itr);

	// Create templates for each vector operation
	u64 rparamc = 3;
	ocrGuid_t copyTemplateGuid, scaleTemplateGuid, addTemplateGuid, triadTemplateGuid;
	ocrEdtTemplateCreate(&copyTemplateGuid, &copyEdt, rparamc, 2);
	ocrEdtTemplateCreate(&scaleTemplateGuid, &scaleEdt, rparamc, 2);
	ocrEdtTemplateCreate(&addTemplateGuid, &addEdt, rparamc, 2);
	ocrEdtTemplateCreate(&triadTemplateGuid, &triadEdt, rparamc, 2);

	// Create run vector operation template guids
	// Note: otherVecOps waits for the previous vector operation to complete hence extra dependency
	ocrGuid_t firstVecOpTemplateGuid, otherVecOpTemplateGuid;
	ocrEdtTemplateCreate(&firstVecOpTemplateGuid, &runVecOpEdt, rparamc, split + 2);
	ocrEdtTemplateCreate(&otherVecOpTemplateGuid, &runVecOpEdt, rparamc, split + 3);

	// Parameter vectors for the vector operations
	u64 cparamv[3] = {cur_itr, COPY_IDX, copyTemplateGuid};
	u64 sparamv[3] = {cur_itr, SCALE_IDX, scaleTemplateGuid};
	u64 aparamv[3] = {cur_itr, ADD_IDX, addTemplateGuid};
	u64 tparamv[3] = {cur_itr, TRIAD_IDX, triadTemplateGuid};

	// Create each vector operation EDT
	ocrGuid_t copyGuid, scaleGuid, addGuid, triadGuid, copyDone, scaleDone, addDone;
	ocrEdtCreate(&copyGuid, firstVecOpTemplateGuid, EDT_PARAM_DEF, cparamv, EDT_PARAM_DEF, NULL_GUID,
	EDT_PROP_FINISH, NULL_GUID, &copyDone);

	ocrEdtCreate(&scaleGuid, otherVecOpTemplateGuid, EDT_PARAM_DEF, sparamv, EDT_PARAM_DEF, NULL_GUID,
	EDT_PROP_FINISH, NULL_GUID, &scaleDone);

	ocrEdtCreate(&addGuid, otherVecOpTemplateGuid, EDT_PARAM_DEF, aparamv, EDT_PARAM_DEF, NULL_GUID,
	EDT_PROP_FINISH, NULL_GUID, &addDone);

	ocrEdtCreate(&triadGuid, otherVecOpTemplateGuid, EDT_PARAM_DEF, tparamv, EDT_PARAM_DEF, NULL_GUID,
	EDT_PROP_FINISH, NULL_GUID, NULL);

	// Add datablock dependencies for each vector operation
	u64 i;
	for (i = 0; i < split + 1; i++) {
		ocrAddDependence((ocrGuid_t) paramv[i], triadGuid, i, DB_MODE_ITW);
		ocrAddDependence((ocrGuid_t) paramv[i], addGuid, i, DB_MODE_ITW);
		ocrAddDependence((ocrGuid_t) paramv[i], scaleGuid, i, DB_MODE_ITW);
	}



	ocrAddDependence(depv[0].guid, triadGuid, split + 2, DB_MODE_RO);
	ocrAddDependence(depv[0].guid, addGuid, split + 2, DB_MODE_RO);
	ocrAddDependence(depv[0].guid, scaleGuid, split + 2, DB_MODE_RO);
	ocrAddDependence(addDone, triadGuid, split + 1, DB_MODE_RO);
	ocrAddDependence(scaleDone, addGuid, split + 1, DB_MODE_RO);
	ocrAddDependence(copyDone, scaleGuid, split + 1, DB_MODE_RO);

	// The dependencies for the first vector operation to run are added last to prevent
	// the other vector operation EDTs from running before the first vector operation
	for (i = 0; i < split + 1; i++)
		ocrAddDependence((ocrGuid_t) paramv[i], copyGuid, i, DB_MODE_ITW);
	ocrAddDependence(depv[0].guid, copyGuid, split + 1, DB_MODE_RO);

	PRINTF("FINISHED PARALLEL STREAM\n");
	return NULL_GUID;
}

ocrGuid_t serialStreamEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
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

ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	// Create local variables of parameters
	struct params * params = (struct params *) depv[0].ptr;
	u64 iterations = params->args.iterations;

	// Initialize next Iteration and Stream templates
	ocrGuid_t nextIterTemplateGuid = params->tmplt.nextIterTemplateGuid;
	ocrGuid_t streamTemplateGuid = params->tmplt.streamTemplateGuid;

	params->cur_itr += 1;
	// Create Stream EDT
	ocrGuid_t streamGuid, streamDone;
	ocrEdtCreate(&streamGuid, streamTemplateGuid, paramc, paramv,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &streamDone);

	// Setup next iteration
	if (params->cur_itr < iterations) {
		ocrGuid_t nextIterGuid;
		ocrEdtCreate(&nextIterGuid, nextIterTemplateGuid, paramc, paramv,
			EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, NULL);
		ocrAddDependence(depv[0].guid, nextIterGuid, 0, DB_MODE_ITW);
		ocrAddDependence(streamDone, nextIterGuid, 1, DB_MODE_RO);
	}

	// Dependencies for Stream
	ocrAddDependence(depv[0].guid, streamGuid, 0, DB_MODE_ITW);

#ifdef DEBUG
	PRINTF("FINISHED ITER\n");
#endif

	return NULL_GUID;
}

ocrGuid_t resultsEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	struct params * params = (struct params *) depv[depc - 1].ptr;
	u64 i, j;
	printf("cur_itr = %llu\n", params->cur_itr);


	u64 db_size = params->args.db_size;
	u64 iterations = params->args.iterations;
	u64 split = params->args.split;
	u64 chunk = params->args.chunk;
	int verify = params->args.verify;
	STREAM_TYPE scalar = params->args.scalar;
	int verbose = params->args.verbose;
	
	// Only used in parallel_stream (dummy value for now)
	int parallel = 1;
	STREAM_TYPE totaltiming[4], timings[iterations][4], avg[4];

	// Used in parallel_stream (ignore for now)
	// u64 actual_split = 0;
	// u64 actual_chunk = 0;
	// if (parallel == 1) {
	// 	actual_split = split;
	// 	split = 1;
	// 	actual_chunk = chunk;
	// 	chunk = db_size;
	// }

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

	// Set initial min and max values for each vector operation to first iteration
	STREAM_TYPE * cur = (STREAM_TYPE *) depv[split].ptr;
	STREAM_TYPE min[4] = {cur[0], cur[1], cur[2], cur[3]};
	STREAM_TYPE max[4] = {cur[0], cur[1], cur[2], cur[3]};

	
	// Sum timings from each iteration
	for (i = 0; i < iterations; i++) {
		timings[i][0] = cur[4 * i];
		timings[i][1] = cur[4 * i + 1];
		timings[i][2] = cur[4 * i + 2];
		timings[i][3] = cur[4 * i + 3];
		// Print results from each iteration if verbose is specified
		if (verbose) {
			PRINTF(HLINE);
			PRINTF("ITERATION %d:\n", i + 1);
			PRINTF("Function       Rate MB/s     Time\n");
			for (j = 0; j < 4; j++)
				PRINTF("%s%12.1f %11.6f\n", label[j], 1.0E-06 * bytes[j] / timings[i][j], timings[i][j]);
		}
		totaltiming[0] += timings[i][0];
		totaltiming[1] += timings[i][1];
		totaltiming[2] += timings[i][2];
		totaltiming[3] += timings[i][3];
		for (j = 0; j < 4; j++) {
			if (timings[i][j] > max[j])
				 max[j] = timings[i][j];
			printf("timings[%llu][%llu] = %f\n", timings[i][j]);
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
	// Ignore for now
	// if (parallel == 1) {
	// 	split = actual_split;
	// 	chunk = actual_chunk;
	// }
	if (strcmp((char *) depv[split+1].ptr, "") != 0)
		export_csv((char *) depv[split + 1].ptr, db_size, iterations, split, scalar, timings, 
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
		// PRINTF("After %d Iterations:\n", iterations);
		STREAM_TYPE * cur = (STREAM_TYPE *) depv[0].ptr;
		for (i = 0; i < split; i++) {
			if (parallel == 1) {
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

	printf("FINISHED RESULTS\n");
	ocrShutdown();
	return NULL_GUID;
}

// ----------------------------------
// INITIALIZATION FUNCTIONS
// ----------------------------------
int initArgs(u64 argc, char ** argv, struct args * a) {
	if (parseOptions(argc, argv, a)) {
		ocrShutdown();
		return 0;
	}
	return 1;
}


void initTemplates(struct templates * t, u64 split, ocrGuid_t (* streamEdt) ()/*u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[])*/, u64 nparamc) {
	ocrEdtTemplateCreate(&t->streamTemplateGuid, streamEdt, split + 1, 1);
	ocrEdtTemplateCreate(&t->nextIterTemplateGuid, &iterEdt, split + 1, 2);
	//ocrEdtTemplateCreate(&t->serialStreamTemplateGuid, &serialStreamEdt, 0, split + 1);
	return;
}

// [a_chunk, b_chunk, c_chunk] [NUM_VECOP * iterations]
void initDbs(ocrGuid_t paramsGuid, struct params * paramsArray, ocrGuid_t * dataGuid) {
	// Create data block for "a, b, c" arrays
	u64 i, j;
	u64 split = paramsArray->args.split;
	u64 chunk = paramsArray->args.chunk;
	u64 iterations = paramsArray->args.iterations;
	ocrGuid_t chunkGuid;
	STREAM_TYPE * chunkArray;
	for (i = 0; i < split; i++) {
		DBCREATE(&chunkGuid,(void **) &chunkArray,
			sizeof(STREAM_TYPE)*(3 * chunk), 0, NULL_GUID, NO_ALLOC);
		for (j = 0; j < chunk; j++) {
			chunkArray[j] = 1.0;
			chunkArray[chunk + j] = 2.0;
			chunkArray[2 * chunk + j] = 0.0;
		}
		dataGuid[i] = chunkGuid;
	}

	// Create data block to store timings
	ocrGuid_t timingGuid;
	STREAM_TYPE * timingArray;
	DBCREATE(&timingGuid,(void **) &timingArray,
		sizeof(STREAM_TYPE)*(NUM_VECOP * iterations), 0, NULL_GUID, NO_ALLOC);
	for (j = 0; j < NUM_VECOP * iterations; j++)
		timingArray[j] = 0.0;
	dataGuid[split] = timingGuid;

	// Add paramsGuid containing arguments to data block
	//dataGuid[split + 1] = paramsGuid;
	return;
}

void initEdts (ocrGuid_t * dataGuids, struct params * paramsArray, ocrGuid_t paramsGuid) {
	// Create iterator and results templates
	u64 split = paramsArray->args.split;
	ocrGuid_t iterGuid, iterDone, iterTemplateGuid, nextIterTemplateGuid, resultsGuid, resultsTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, split+ 1, 1);//split + 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, 0, split + 4);

	// Add next iterator and stream template guids to paramv 
	//u64 paramv[1]={(u64) dataGuids[0]};//, (u64) dataGuids[1]};
	u64 paramv[split + 1];
	u64 i;
	for (i = 0; i < split + 1; i++)
		paramv[i] = (u64) dataGuids[i];

	// Create iterator and results EDTs
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, paramv, EDT_PARAM_DEF, NULL_GUID,
		EDT_PROP_FINISH, NULL_GUID, &iterDone);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, NULL_GUID, EDT_PARAM_DEF, NULL_GUID,
		EDT_PROP_NONE, NULL_GUID, NULL);

	// Create datablock for export file name
	ocrGuid_t efileGuid;
	//printf("export file_name = %s\n", paramsArray->args.efile);
	char * efileArray;
	DBCREATE(&efileGuid, (void **) &efileArray, sizeof(char) * sizeof(paramsArray->args.efile), 0, NULL_GUID, NO_ALLOC);
	strcpy(efileArray, paramsArray->args.efile);
	//printf("export file array = %s\n", efileArray);

	// Dependencies for iter and results
	for(i = 0; i < split + 1; i++)
		ocrAddDependence(dataGuids[i], resultsGuid, i, DB_MODE_RO);
	ocrAddDependence(efileGuid, resultsGuid, split + 1, DB_MODE_RO);
	ocrAddDependence(iterDone, resultsGuid, split + 2, DB_MODE_RO);
	ocrAddDependence(paramsGuid, resultsGuid, split + 3, DB_MODE_RO);
	
	ocrAddDependence(paramsGuid, iterGuid, 0, DB_MODE_ITW);
	//for (i = 0; i < split + 1; i++)
	//	ocrAddDependence(dataGuids[i], iterGuid, i, DB_MODE_ITW);

#ifdef DEBUG
	PRINTF("FINISHED CREATE EDTS\n");
#endif
	return;
}
