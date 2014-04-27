
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stack>

#include "OCRSparseMatrix.h"
#include "OCRDenseMatrix.h"

#define __OCR__
#include <ocr.h>

ocrGuid_t testEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	OCRDenseMatrix *dm = new OCRDenseMatrix(depv[0].ptr,depv[0].guid);
	dm->print();
	PRINTF("\n");
	OCRSparseMatrix *spm = new OCRSparseMatrix(dm,0.0);
	spm->print();
	PRINTF("\n");

	OCRSparseMatrix *sparseProduct = dm->multiplySparse(spm);
	sparseProduct->print();
	PRINTF("\n");
	OCRDenseMatrix *denseProduct = dm->multiplyDense(spm);
	denseProduct->print();
	PRINTF("\n");
	PRINTF("Finished successfully\n");
	
	delete spm;
	delete dm;
	delete sparseProduct;
	delete denseProduct;
	ocrShutdown(); // This is the last EDT to execute
	return NULL_GUID;
}

extern "C" ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	ocrGuid_t testTempEdt, testRunEdt;
	ocrEdtTemplateCreate(&testTempEdt, &testEdt, 0, 1);
	OCRDenseMatrix *dm = new OCRDenseMatrix(10, 10);
	dm->setElement(0,0,1.0);
	dm->setElement(1,2,1.0);
	dm->setElement(2,4,1.0);
	dm->setElement(3,6,1.0);
	dm->setElement(4,9,8.0);
	dm->setElement(4,1,1.0);
	dm->setElement(5,3,1.0);
	dm->setElement(6,2,2.0);
	dm->setElement(6,5,1.0);
	dm->setElement(8,8,1.0);
	dm->print();
	PRINTF("\n");

	ocrGuid_t matrixData = dm->getDatablock();

	ocrEdtCreate(&testRunEdt, testTempEdt, EDT_PARAM_DEF, NULL, EDT_PARAM_DEF, &matrixData, EDT_PROP_NONE, NULL_GUID, NULL);
	//delete dm;

	return NULL_GUID;
}
