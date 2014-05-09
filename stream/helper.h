#include <string.h>
#include <sys/time.h> 

#ifndef __OCR__
#	define __OCR__
#	include "ocr.h"
#endif

#ifndef STREAM_TYPE
#	define STREAM_TYPE double
#endif

double mysecond() {
	struct timeval tp;
	struct timezone tzp;
	int i = gettimeofday(&tp, &tzp);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}

void export_csv(char * name, u64 db_size, u64 iterations, u64 split, STREAM_TYPE scalar, STREAM_TYPE trials[][4], STREAM_TYPE avg[4]) {
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
	fprintf(f, "%llu, %llu, %llu, %.2f, %f, %f, %f, %f\n", db_size, iterations, split, scalar, avg[0], avg[1], avg[2], avg[3]);

	fclose(f);
	return;
}
