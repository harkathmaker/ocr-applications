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

#define NPARAMC 1

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

	// Parse getopt commands and add values into args sub-struct
	if (initArgs(argc, argv, &paramsArray->args) == 0)
		return NULL_GUID;

	// Create templates and set paramv
	initTemplates(&paramsArray->tmplt, paramsArray->args.split, &parallelStreamEdt, NPARAMC);

	// Create one datablock containing addresses to arrays, timings, export file name, and parameters
	ocrGuid_t dataGuids[paramsArray->args.split + 1];
	//ocrGuid_t * dataGuidsArray;
	//DBCREATE(&dataGuids,(void **) &dataGuidsArray,
	//	sizeof(dataGuids) * (paramsArray->args.split + 2), 0, NULL_GUID, NO_ALLOC);
	initDbs(paramsGuid, paramsArray, dataGuids);

	// Create Parallel Stream Template Guid
	 //ocrGuid_t parallelStreamTemplateGuid;
	// ocrEdtTemplateCreate(&parallelStreamTemplateGuid, &parallelStreamEdt, 0, paramsArray->args.split + 1);

	// Create iter and results EDTs
	initEdts(dataGuids, paramsArray, paramsGuid);

	PRINTF("FINISHED MAIN\n");
	return NULL_GUID; 
}