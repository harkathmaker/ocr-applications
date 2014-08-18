#include <string.h>
#include <sys/time.h>

#ifndef __OCR__
	#define __OCR__
	#include "ocr.h"
#endif

#ifndef STREAM_TYPE
	#define STREAM_TYPE double
#endif

#define NUM_VECOP 4

// #define DEBUG

#define HLINE "--------------------------------------------------------------\n"

char * label[4] = {"Copy:      ", "Scale:     ", "Add:       ", "Triad:     "};

// ----------------------------------
// GLOBAL STRUCTS
// ----------------------------------

struct args {
	u64  chunk;
	u64  db_size;
	char efile[100];
	u64  iterations;
	u64  split;
	int  verify;
	int  verbose;
	STREAM_TYPE scalar;
};

struct templates {
	ocrGuid_t nextIterTemplateGuid;
	ocrGuid_t streamTemplateGuid;
	ocrGuid_t copyTemplateGuid;
	ocrGuid_t scaleTemplateGuid;
	ocrGuid_t addTemplateGuid;
	ocrGuid_t triadTemplateGuid;
	ocrGuid_t firstVecOpTemplateGuid;
	ocrGuid_t otherVecOpTemplateGuid;
};

struct params{
	u64 cur_itr;
	struct templates tmplt;
	struct args args;
};

// ----------------------------------
// TIMERS
// ----------------------------------

// Original timer function used in native STREAM
double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);

	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void export_csv(char * name, u64 db_size, u64 iterations, u64 split,
	STREAM_TYPE scalar, STREAM_TYPE trials[][4], STREAM_TYPE cavg,
	STREAM_TYPE savg, STREAM_TYPE aavg, STREAM_TYPE tavg) {

	u64 i;
	FILE * f = fopen(name, "a");

	if (f == NULL) {
		PRINTF("Error creating export file 	.");
		exit(1);
	}

	fprintf(f, "%llu, %llu, %llu, %.2f, %12.1f, %12.1f, %12.1f, %12.1f\n",
		db_size, iterations, split, scalar, cavg, savg, aavg, tavg);

	fclose(f);

	return;
}

ocrGuid_t endTimerEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	STREAM_TYPE stop = mysecond();
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	int time_idx = NUM_VECOP * (paramv[0] - 1) + paramv[1];
	data[time_idx] = stop - data[time_idx];

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

#ifdef DEBUG
	PRINTF("FINISHED COPY\n");
#endif

	return NULL_GUID;
}

ocrGuid_t scaleEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[0];
	STREAM_TYPE scalar = paramv[1];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	for (i = 0; i < chunk; i++)
		data[chunk + i] = scalar * data[2 * chunk + i];

#ifdef DEBUG
	PRINTF("FINISHED SCALE\n");
#endif

	return NULL_GUID;
}

ocrGuid_t addEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[0];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	for (i = 0; i < chunk; i++)
		data[2 * chunk + i] = data[i] + data[chunk + i];

#ifdef DEBUG
	PRINTF("FINISHED ADD\n");
#endif

	return NULL_GUID;
}

ocrGuid_t triadEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i;
	u64 chunk = paramv[0];
	STREAM_TYPE scalar = paramv[1];
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;

	for (i = 0; i < chunk; i++)
		data[i] = data[chunk + i] + scalar * data[2 * chunk + i];

#ifdef DEBUG
	PRINTF("FINISHED TRIAD\n");
#endif

	return NULL_GUID;
}

// ----------------------------------
// SETUP + RESULTS EDTS
// ----------------------------------

ocrGuid_t runVecOpEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	// Create local variables of parameters
	u64 cur_itr= paramv[0];
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
	ocrEdtCreate(&endTimerGuid, endTimerTemplateGuid, EDT_PARAM_DEF, tparamv,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, NULL);
	ocrAddDependence((ocrGuid_t) depv[split].guid, endTimerGuid, 0, DB_MODE_ITW);

	// Create "Start" Event signalling all instances of the vector operation to run
	ocrEventCreate(&allVecOpReady, OCR_EVENT_STICKY_T, true);

	// Create "split" number of selected vector operation EDTs
	u64 i;
	u64 nparamv[2] = {chunk, scalar};
	for (i = 0; i < split; i++) {
		ocrEdtCreate(&vecOpGuid, vecOpTemplateGuid, EDT_PARAM_DEF, nparamv,
			EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &vecOpDone);
		// Run the vector operation EDTs once all of the vector operation EDTs have
		// been created
		ocrAddDependence((ocrGuid_t) depv[i].guid, vecOpGuid, 0, DB_MODE_ITW);
		ocrAddDependence(allVecOpReady, vecOpGuid, 1, DB_MODE_RO);
		ocrAddDependence(vecOpDone, endTimerGuid, i + 1, DB_MODE_RO);
	}

#ifdef DEBUG
	PRINTF("FINISHED CREATING ALL VECTOR OPERATIONS\n");
#endif

	// Start the timer and run "split" number of the vector operation EDTs
	u64 timer_idx =  NUM_VECOP * (cur_itr- 1) + vecOp;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[split].ptr;
	data[timer_idx] = mysecond();

	// Start all selected vector operation EDTs
	ocrEventSatisfy(allVecOpReady, NULL_GUID);

#ifdef DEBUG
	PRINTF("FINISHED RUN VECOP\n");
#endif

	return NULL_GUID;
}

ocrGuid_t parallelStreamEdt(u32 paramc, u64 * paramv, u32 depc,
	ocrEdtDep_t depv[]) {

	// Create local variables of parameters
	struct params * params = (struct params *) depv[0].ptr;
	u64 cur_itr = params->cur_itr;
	u64 split = params->args.split;

	ocrGuid_t copyTemplateGuid = params->tmplt.copyTemplateGuid;
	ocrGuid_t scaleTemplateGuid = params->tmplt.scaleTemplateGuid;
	ocrGuid_t addTemplateGuid = params->tmplt.addTemplateGuid;
	ocrGuid_t triadTemplateGuid = params->tmplt.triadTemplateGuid;

	ocrGuid_t firstVecOpTemplateGuid = params->tmplt.firstVecOpTemplateGuid;
	ocrGuid_t otherVecOpTemplateGuid = params->tmplt.otherVecOpTemplateGuid;

	// Parameter vectors for the vector operations
	u64 cparamv[3] = {cur_itr, 0, copyTemplateGuid};
	u64 sparamv[3] = {cur_itr, 1, scaleTemplateGuid};
	u64 aparamv[3] = {cur_itr, 2, addTemplateGuid};
	u64 tparamv[3] = {cur_itr, 3, triadTemplateGuid};

	// Create each vector operation EDT
	ocrGuid_t copyGuid, scaleGuid, addGuid, triadGuid, copyDone, scaleDone, addDone;
	ocrEdtCreate(&copyGuid, firstVecOpTemplateGuid, EDT_PARAM_DEF, cparamv,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &copyDone);

	ocrEdtCreate(&scaleGuid, otherVecOpTemplateGuid, EDT_PARAM_DEF, sparamv,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &scaleDone);

	ocrEdtCreate(&addGuid, otherVecOpTemplateGuid, EDT_PARAM_DEF, aparamv,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &addDone);

	ocrEdtCreate(&triadGuid, otherVecOpTemplateGuid, EDT_PARAM_DEF, tparamv,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, NULL);

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

#ifdef DEBUG
	PRINTF("FINISHED PARALLEL STREAM\n");
#endif

	return NULL_GUID;
}

ocrGuid_t serialStreamEdt(u32 paramc, u64 * paramv, u32 depc,
	ocrEdtDep_t depv[]) {

	ocrGuid_t serialGuid;
	struct params * params = (struct params *) depv[0].ptr;
	ocrGuid_t runSerialStreamTemplateGuid = params->tmplt.firstVecOpTemplateGuid;
	u64 split = params->args.split;
	ocrEdtCreate(&serialGuid, runSerialStreamTemplateGuid, 0, NULL_GUID,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, NULL);

	u64 i;
	for (i = 0; i < split + 1; i++)
		ocrAddDependence((ocrGuid_t) paramv[i], serialGuid, i, DB_MODE_RO);
	ocrAddDependence(depv[0].guid, serialGuid, split + 1, DB_MODE_RO);

#ifdef DEBUG
	PRINTF("FINISHED SERIAL STREAM\n");
#endif

	return NULL_GUID;
}

ocrGuid_t runSerialStreamEdt(u32 paramc, u64 * paramv, u32 depc,
	ocrEdtDep_t depv[]) {

	struct params * params = (struct params *) depv[depc - 1].ptr;
	u64 i;
	u64 db_size = params->args.db_size;
	u64 begin = 0;
	u64 end = params->args.db_size;
	STREAM_TYPE scalar = params->args.scalar;
	STREAM_TYPE * data = (STREAM_TYPE *) depv[0].ptr;
	STREAM_TYPE * timings = (STREAM_TYPE *) depv[1].ptr;
	STREAM_TYPE start, stop;
	u64 cur_itr = params->cur_itr;

	// COPY
	start = mysecond();
	for (i = begin; i < end; i++) {
		data[2 * db_size + i] = data[i];
		//printf("data[%d] = %f\n", i, data[i]);
	}
	stop = mysecond();
	timings[4 * (cur_itr - 1)] = stop - start;

	// SCALE
	start = mysecond();
	for (i = begin; i < end; i++)
		data[db_size + i] = scalar * data[2 * db_size + i];
	stop = mysecond();
	timings[4 * (cur_itr - 1) + 1] = stop - start;

	// ADD
	start = mysecond();
	for (i = begin; i < end; i++)
		data[2 * db_size + i] = data[i] + data[db_size + i];
	stop = mysecond();
	timings[4 * (cur_itr - 1) + 2] = stop - start;

	// TRIAD
	start = mysecond();
	for (i = begin; i < end; i++)
		data[i] = data[db_size + i] + scalar * data[2 * db_size + i];
	stop = mysecond();
	timings[4 * (cur_itr - 1) + 3] = stop - start;

#ifdef DEBUG
	PRINTF("FINISHED RUN SERIAL STREAM\n");
#endif

	return NULL_GUID;
}

ocrGuid_t iterEdt(u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[]) {
	// Create local variables of parameters
	struct params * params = (struct params *) depv[0].ptr;
	u64 iterations = params->args.iterations;

	// Initialize next Iteration and Stream templates
	ocrGuid_t nextIterTemplateGuid = params->tmplt.nextIterTemplateGuid;
	ocrGuid_t streamTemplateGuid = params->tmplt.streamTemplateGuid;

	// Increment current iteration number
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
	u64 db_size = params->args.db_size;
	u64 iterations = params->args.iterations;
	u64 split = params->args.split;
	u64 chunk = params->args.chunk;
	int verify = params->args.verify;
	STREAM_TYPE scalar = params->args.scalar;
	int verbose = params->args.verbose;

	STREAM_TYPE totaltiming[4], timings[iterations][4], avg[4];

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
	u64 i, j;
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
				PRINTF("%s%12.1f %11.6f\n", label[j],
					1.0E-06 * bytes[j] / timings[i][j], timings[i][j]);
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
		PRINTF("%s%12.1f  %11.6f  %11.6f  %11.6f\n", label[i],
			1.0E-06 * bytes[i] / avg[i], avg[i], min[i], max[i]);
	PRINTF(HLINE);

	// Export to CSV
	if (strcmp((char *) depv[split+1].ptr, "") != 0)
		export_csv((char *) depv[split + 1].ptr, db_size, iterations, split,
			scalar, timings,
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
		STREAM_TYPE * cur;
		for (i = 0; i < split; i++) {
			cur = (STREAM_TYPE *) depv[i].ptr;
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

#ifdef DEBUG
	printf("FINISHED RESULTS\n");
#endif

	ocrShutdown();

	return NULL_GUID;
}

// ----------------------------------
// INITIALIZATION FUNCTIONS
// ----------------------------------

int initArgs(u64 argc, char ** argv, struct args * a) {
	if (parseOptions(argc, argv, a))
		return 0;

	return 1;
}

void initTemplates(struct templates * t, u64 split,
	ocrGuid_t (* streamEdt) (u32 paramc, u64 * paramv, u32 depc, ocrEdtDep_t depv[])) {

	ocrEdtTemplateCreate(&t->nextIterTemplateGuid, &iterEdt, split + 1, 2);
	ocrEdtTemplateCreate(&t->streamTemplateGuid,   streamEdt, split + 1, 1);

	u64 rparamc = 3;
	ocrEdtTemplateCreate(&t->copyTemplateGuid,  &copyEdt,  rparamc, 2);
	ocrEdtTemplateCreate(&t->scaleTemplateGuid, &scaleEdt, rparamc, 2);
	ocrEdtTemplateCreate(&t->addTemplateGuid,   &addEdt,   rparamc, 2);
	ocrEdtTemplateCreate(&t->triadTemplateGuid, &triadEdt, rparamc, 2);

	ocrEdtTemplateCreate(&t->firstVecOpTemplateGuid, &runVecOpEdt, rparamc,
		split + 2);
	ocrEdtTemplateCreate(&t->otherVecOpTemplateGuid, &runVecOpEdt, rparamc,
		split + 3);

	return;
}

void createDbs(ocrGuid_t paramsGuid, struct params * paramsArray, ocrGuid_t * dataGuid) {

	// Create data block for "a, b, c" arrays
	u64 i, j;
	ocrGuid_t      chunkGuid;
	STREAM_TYPE *  chunkArray;
	u64 split      = paramsArray->args.split;
	u64 chunk      = paramsArray->args.chunk;
	u64 iterations = paramsArray->args.iterations;

	for (i = 0; i < split; i++) {
		DBCREATE(&chunkGuid,(void **) &chunkArray,
			sizeof(STREAM_TYPE) * (3 * chunk), 0, NULL_GUID, NO_ALLOC);
		for (j = 0; j < chunk; j++) {
			chunkArray[j]             = 1.0;
			chunkArray[chunk + j]     = 2.0;
			chunkArray[2 * chunk + j] = 0.0;
		}
		dataGuid[i] = chunkGuid;
	}

	// Create data block to store timings
	ocrGuid_t     timingGuid;
	STREAM_TYPE * timingArray;
	DBCREATE(&timingGuid,(void **) &timingArray,
		sizeof(STREAM_TYPE)*(NUM_VECOP * iterations), 0, NULL_GUID, NO_ALLOC);
	for (j = 0; j < NUM_VECOP * iterations; j++)
		timingArray[j] = 0.0;
	dataGuid[split] = timingGuid;

	return;
}

void startStream (ocrGuid_t * dataGuids, struct params * paramsArray,
	ocrGuid_t paramsGuid) {

	// Create iterator and results templates
	u64 split = paramsArray->args.split;
	ocrGuid_t iterGuid, iterDone, iterTemplateGuid, nextIterTemplateGuid,
		resultsGuid, resultsTemplateGuid;
	ocrEdtTemplateCreate(&iterTemplateGuid, &iterEdt, split+ 1, 1);
	ocrEdtTemplateCreate(&resultsTemplateGuid, &resultsEdt, 0, split + 4);

	// Add next iterator and stream template guids to paramv
	u64 i;
	u64 paramv[split + 1];
	for (i = 0; i < split + 1; i++)
		paramv[i] = (u64) dataGuids[i];

	// Create iterator and results EDTs
	ocrEdtCreate(&iterGuid, iterTemplateGuid, EDT_PARAM_DEF, paramv,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_FINISH, NULL_GUID, &iterDone);
	ocrEdtCreate(&resultsGuid, resultsTemplateGuid, EDT_PARAM_DEF, NULL_GUID,
		EDT_PARAM_DEF, NULL_GUID, EDT_PROP_NONE, NULL_GUID, NULL);

	// Create datablock for export file name
	ocrGuid_t efileGuid;
	char * efileArray;
	DBCREATE(&efileGuid, (void **) &efileArray,
		sizeof(char) * sizeof(paramsArray->args.efile), 0, NULL_GUID, NO_ALLOC);
	strcpy(efileArray, paramsArray->args.efile);

	// Dependencies for results
	for(i = 0; i < split + 1; i++)
		ocrAddDependence(dataGuids[i], resultsGuid, i, DB_MODE_RO);
	ocrAddDependence(efileGuid,  resultsGuid, split + 1, DB_MODE_RO);
	ocrAddDependence(iterDone,   resultsGuid, split + 2, DB_MODE_RO);
	ocrAddDependence(paramsGuid, resultsGuid, split + 3, DB_MODE_RO);

	// Dependency for iterator
	ocrAddDependence(paramsGuid, iterGuid, 0, DB_MODE_ITW);

#ifdef DEBUG
	PRINTF("FINISHED START STREAM\n");
#endif

	return;
}
