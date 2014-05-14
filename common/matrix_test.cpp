
#define _USE_MATH_DEFINES

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stack>

#include "OCRSparseMatrix.h"
#include "OCRDenseMatrix.h"
#include "OCRDiagMatrix.h"

#define __OCR__
#include <ocr.h>

ocrGuid_t testEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	OCRDenseMatrix *dm = new OCRDenseMatrix(depv[0].ptr,depv[0].guid);
	OCRDiagMatrix *diag = new OCRDiagMatrix(depv[1].ptr,depv[1].guid);
	PRINTF("Dense:\n");
	dm->print();
	PRINTF("Diag:\n");
	diag->print();
	PRINTF("\n");
	OCRSparseMatrix *spm = new OCRSparseMatrix(dm,0.0);
	PRINTF("Sparse conversion:\n");
	spm->print();
	PRINTF("\n");

	OCRSparseMatrix *sparseProduct = dm->multiplySparse(spm);
	PRINTF("Sparse product:\n");
	sparseProduct->print();
	PRINTF("\n");
	PRINTF("Dense product:\n");
	OCRDenseMatrix *denseProduct = dm->multiplyDense(spm);
	denseProduct->print();
	PRINTF("\n");
	PRINTF("Finished successfully\n");

	OCRMatrix *inv = dm->getInverse();
	PRINTF("Diagonal Inverse:\n");
	inv->print();

	PRINTF("Determinant: %f\n",dm->getDeterminant());
	
	delete inv;
	delete spm;
	delete dm;
	delete sparseProduct;
	delete denseProduct;
	ocrShutdown(); // This is the last EDT to execute
	return NULL_GUID;
}

extern "C" ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]) {
	ocrGuid_t testTempEdt, testRunEdt;
	ocrEdtTemplateCreate(&testTempEdt, &testEdt, 0, 2);
	OCRDenseMatrix *dm = new OCRDenseMatrix(3, 3);
//	dm->setElement(0,0,1.0);
//	dm->setElement(1,2,1.0);
//	dm->setElement(2,4,1.0);
//	dm->setElement(3,6,1.0);
//	dm->setElement(4,9,8.0);
//	dm->setElement(4,1,1.0);
//	dm->setElement(5,3,1.0);
//	dm->setElement(6,2,2.0);
//	dm->setElement(6,5,1.0);
//	dm->setElement(7,9,3.2);
//	dm->setElement(8,8,1.0);
//	dm->setElement(9,7,2.25);
	dm->setElement(0,0,0); dm->setElement(0,1,2); dm->setElement(0,2,3);
	dm->setElement(1,0,0); dm->setElement(1,1,0); dm->setElement(1,2,3);
	dm->setElement(2,0,1); dm->setElement(2,1,0); dm->setElement(2,2,8);
	dm->print();
	PRINTF("\n");

	OCRDiagMatrix *diag = new OCRDiagMatrix(4);
	diag->setElement(0,0,1.0);
	diag->setElement(1,1,26);
	diag->setElement(2,2,3.14);
	diag->setElement(3,3,141);
	diag->print();

	ocrGuid_t matrixData[2] = { dm->getDatablock(), diag->getDatablock() };

	ocrEdtCreate(&testRunEdt, testTempEdt, EDT_PARAM_DEF, NULL, EDT_PARAM_DEF, matrixData, EDT_PROP_NONE, NULL_GUID, NULL);
	//delete dm;

	return NULL_GUID;
}
