#include "options.h"

#ifndef __OCR__
	#define __OCR__
	#include "ocr.h"
#endif

#ifndef __HELPER__
	#define __HELPER__
	#include "helper.h"
#endif

#ifndef STREAM_TYPE
	#define STREAM_TYPE double
#endif

#define DEBUG

/*
	NOTES:
	One large datablock is used with ranges of indexes corresponding to the following:
			  a --> [0 to (db_size - 1)]
			  b --> [db_size to (2 * db_size - 1)]
			  c --> [2 * db_size to (3 * db_size - 1)]
		timings --> [3 * db_size to (3 * db_size + iterations - 1)]
*/

ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	u64 i, j, argc = getArgc(depv[0].ptr);
	char * argv[argc];
	for (i = 0; i < argc; i++)
		argv[i] = getArgv(depv[0].ptr, i);

	// Create params data block struct containing argument data
	ocrGuid_t paramsGuid;
	struct params * paramsArray;
	DBCREATE(&paramsGuid,(void **) &paramsArray,
		sizeof(struct params), 0, NULL_GUID, NO_ALLOC);

	// Initialize params struct and set current iteration to 0
	paramsArray->cur_itr = 0;
	paramsArray->args.split = -1;

	// Parse getopt commands and add values into args sub-struct
	if (initArgs(argc, argv, &paramsArray->args) == 0) {
		ocrDbRelease(paramsGuid);
		ocrDbDestroy(paramsGuid);
		ocrShutdown();
		return NULL_GUID;
	}
	paramsArray->args.split = 1;

	// Create templates and set paramv
	ocrGuid_t nextIterTemplateGuid, streamTemplateGuid, firstVecOpTemplateGuid;
	u64 split = paramsArray->args.split;
	ocrEdtTemplateCreate(&paramsArray->tmplt.nextIterTemplateGuid, iterEdt, split + 1, 2);
	ocrEdtTemplateCreate(&paramsArray->tmplt.streamTemplateGuid, serialStreamEdt, split + 1, 1);
	ocrEdtTemplateCreate(&paramsArray->tmplt.firstVecOpTemplateGuid, runSerialStreamEdt, 0, split + 2);

	// Create one datablock containing addresses to arrays, timings, export file name, and parameters
	ocrGuid_t dataGuids[paramsArray->args.split + 1];
	createDbs(paramsGuid, paramsArray, dataGuids);

	// Create iter and results EDTs
	startStream(dataGuids, paramsArray, paramsGuid);

	return NULL_GUID;
}
