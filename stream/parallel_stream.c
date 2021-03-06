#include "options.h"

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

ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	// Get argc and argv values from user input
	u64 i, j, argc = getArgc(depv[0].ptr);
	char * argv[argc];
	for (i = 0; i < argc; i++)
		argv[i] = getArgv(depv[0].ptr, i);

	// Create params data block struct to store argument data
	ocrGuid_t paramsGuid;
	struct params * paramsArray;
	DBCREATE(&paramsGuid,(void **) &paramsArray,
		sizeof(struct params), 0, NULL_GUID, NO_ALLOC);

	// Set current iteration to 0
	paramsArray->cur_itr = 0;

	// Parse getopt commands and add values into args sub-struct
	// Release and destroy data block if help selected or user gives unknown option
	if (initArgs(argc, argv, &paramsArray->args) == 0) {
		ocrDbRelease(paramsGuid);
		ocrDbDestroy(paramsGuid);
		ocrShutdown();
		return NULL_GUID;
	}

	// Initialize templates
	initTemplates(&paramsArray->tmplt, paramsArray->args.split, &parallelStreamEdt);

	// Create one datablock containing addresses to arrays, timings, export file name,
	// and parameters
	ocrGuid_t dataGuids[paramsArray->args.split + 1];
	createDbs(paramsGuid, paramsArray, dataGuids);

	// Create iter and results EDTs
	startStream(dataGuids, paramsArray, paramsGuid);

#ifdef DEBUG
	PRINTF("FINISHED MAIN\n");
#endif

	return NULL_GUID;
}
