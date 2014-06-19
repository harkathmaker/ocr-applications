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

#ifndef __HELPER__
	#define __HELPER__
	#include "helper.h"
#endif

#define NPARAMC 13

ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j, argc = getArgc(depv[0].ptr);
	char * argv[argc];
	for (i = 0; i < argc; i++)
		argv[i] = getArgv(depv[0].ptr, i);

	// Parse getopt commands and fill values into args struct
	struct args a;
	if (init_args(argc, argv, &a) == 0) {
		return NULL_GUID;
	}

	// Create templates
	struct templates t;
	create_templates(&t, a.split, NPARAMC);

	// Initialize param vector to be passed around edts
	u64 nparamv[NPARAMC] = {1, a.db_size, a.iterations, a.split, a.chunk, a.scalar, t.pipeExecTemplateGuid, 
		t.nextIterTemplateGuid, t.finePipelineTemplateGuid, t.copyTemplateGuid, t.scaleTemplateGuid, 
		t.addTemplateGuid, t.triadTemplateGuid};

	// Create datablocks
	ocrGuid_t dataGuids[a.split];
	create_multi_dbs(a.split, a.chunk, a.iterations, dataGuids);

	// create iter and results edts
	create_multi_db_edts(NPARAMC, nparamv, dataGuids, &a);
	//ocrShutdown();
	return NULL_GUID; 
}


int test_getopt(u64 argc, char ** argv) {
	struct args a;
	init_args(argc, argv, &a);
	printf("db_size = %llu\n", a.db_size);
	printf("iterations = %llu\n", a.iterations);
	ocrShutdown();
}