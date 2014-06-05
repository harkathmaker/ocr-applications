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

	// CHANGE: MAY NEED TO MODIFY SPLIT TO PARAMC - 1
	// Export to CSV
	if (parallel == 1) {
		split = actual_split;
		chunk = actual_chunk;
	}
	if (strcmp((char *) depv[depc - 2].ptr, "") != 0) 
		export_csv((char *) depv[split].ptr, db_size, iterations, split, scalar, timings, 
				   1.0E-06 * bytes[0] / avg[0], 1.0E-06 * bytes[1] / avg[1],
				   1.0E-06 * bytes[2] / avg[2], 1.0E-06 * bytes[3] / avg[3]);

	// CHANGE: MODIFY FOR SERIAL VERSIONS
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