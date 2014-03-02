/*
 * File:   main.cpp
 * Author: commlaptop
 *
 * Created on February 7, 2014, 1:38 PM
 */

#include <iostream>
#include <cmath>
#include <time.h>

extern "C" {
#include <pthread.h>
}

#define DEBUG_MESSAGES

#ifndef __OCR__
#define __OCR__
#endif

#include "ocr.h"

#include "matrix.h"

#define K_ITERATIONS 10

using namespace std;

//Parameters: result_event, A_rows, A_columns
//Dependencies: A_event, scalar_event (double)
//Output: NULL_GUID

extern "C" ocrGuid_t scaleEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "scaleEdt()" << endl;
#endif
	double *a = (double*) depv[0].ptr;
	double *b = (double*) depv[1].ptr;
	ocrGuid_t dataBlock = matrixScale(a, (int) paramv[1], (int) paramv[2], b);
	ocrEventSatisfy((ocrGuid_t) paramv[0], dataBlock);
#ifdef DEBUG_MESSAGES
	cout << "scaleEdt(): result satisfied" << endl;
#endif
	return NULL_GUID;
}

//Parameters: result_event, A_rows, A_columns, B_rows, B_columns
//Dependencies: A_event, B_event
//Output: NULL_GUID

extern "C" ocrGuid_t productEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "productEdt()" << endl;
#endif
	ocrGuid_t dataBlock = matrixProduct((double*) depv[0].ptr, (int) paramv[1], (int) paramv[2],
		(double*) depv[1].ptr, (int) paramv[3], (int) paramv[4]);
	ocrEventSatisfy((ocrGuid_t)paramv[0], dataBlock);
#ifdef DEBUG_MESSAGES
	cout << "productEdt(): result satisfied" << endl;
#endif
	return NULL_GUID;
}

//Parameters: result_event, A_rows, A_columns
//Dependencies: A_event
//Output: NULL_GUID

extern "C" ocrGuid_t transposeEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "transposeEdt()" << endl;
#endif
	ocrGuid_t dataBlock = matrixTranspose((double*) depv[0].ptr, (int) paramv[1], (int) paramv[2]);
	ocrEventSatisfy((ocrGuid_t) paramv[0], dataBlock);
#ifdef DEBUG_MESSAGES
	cout << "transposeEdt(): result satisfied" << endl;
#endif
	return NULL_GUID;
}

//Parameters: result_event, A_rows, A_columns
//Dependencies: A_event, B_event (A and B need same rows/columns)
//Output: NULL_GUID

extern "C" ocrGuid_t addEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "addEdt()" << endl;
#endif
	ocrGuid_t dataBlock = matrixAdd((double*) depv[0].ptr, (int) paramv[1], (int) paramv[2], (double*) depv[1].ptr);
	ocrEventSatisfy((ocrGuid_t) paramv[0], dataBlock);
#ifdef DEBUG_MESSAGES
	cout << "addEdt(): result satisfied" << endl;
#endif
	return NULL_GUID;
}

//Parameters: result_event, A_rows, A_columns
//Dependencies: A_event, B_event (A and B need same rows/columns)
//Output: NULL_GUID

extern "C" ocrGuid_t subtractEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "subtractEdt()" << endl;
#endif
	ocrGuid_t dataBlock = matrixSubtract((double*) depv[0].ptr, (int) paramv[1], (int) paramv[2],
		(double*) depv[1].ptr);
	ocrEventSatisfy((ocrGuid_t) paramv[0], dataBlock);
#ifdef DEBUG_MESSAGES
	cout << "subtractEdt(): result satisfied" << endl;
#endif
	return NULL_GUID;
}

//Parameters: result_event
//Dependencies: A_event (double), B_event (double)
//Output: NULL_GUID

extern "C" ocrGuid_t divideEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "printEdt()" << endl;
#endif
	double* r;
	ocrGuid_t dataBlock;
	DBCREATE(&dataBlock, (void**) &r, 1 * sizeof (double), DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	double* a = (double*) depv[0].ptr;
	double* b = (double*) depv[1].ptr;
	r[0] = a[0] / b[0];
	ocrEventSatisfy((ocrGuid_t) paramv[0], dataBlock);
#ifdef DEBUG_MESSAGES
	cout << "divideEdt(): result satisfied" << endl;
#endif
	return NULL_GUID;
}

//Parameters: x_rows
//Dependencies: A_event
//Output: NULL_GUID
//Not thread safe

extern "C" ocrGuid_t printEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "printEdt()" << endl;
#endif
	int rows = (int) paramv[0];
	double* r = (double*) depv[0].ptr;

	cout << "Solution:" << endl;
	for (int i = 0; i < rows; i++)
		cout << i << ",1: " << r[i] << endl;

	cout << "OCR Shutdown..." << endl;
	ocrShutdown();
	return NULL_GUID;
}

ocrGuid_t createEdt(ocrGuid_t templateGuid, ocrGuid_t outputGuid,
	int A_rows = 0, int A_columns = 0, int B_rows = 0, int B_columns = 0)
{
#ifdef DEBUG_MESSAGES
	cout << "createEdt()" << endl;
#endif
	u64 nparamv[5];
	nparamv[0] = outputGuid;
	nparamv[1] = A_rows;
	nparamv[2] = A_columns;
	nparamv[3] = B_rows;
	nparamv[4] = B_columns;
#ifdef DEBUG_MESSAGES
	cout << "createEdt(): initialized parameters" << endl;
#endif
	ocrGuid_t myEdt;
	ocrEdtCreate(&myEdt, templateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL_GUID, EDT_PROP_NONE, NULL_GUID, NULL);
	return myEdt;
}

//Parameters: A_rows, A_columns, X_old_rows, X_old_columns, k, doneEvent
//Dependencies: A_db, B_db, X_old_db, P_old_db, R_old_db
//Output: Guid data block with result of x_new after K_ITERATIONS

extern "C" ocrGuid_t CgEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t result = (ocrGuid_t) paramv[0];
	int A_rows = (int) paramv[1];
	int A_columns = (int) paramv[2];
	int X_old_rows = (int) paramv[3];
	int X_old_columns = (int) paramv[4];
	int k = (int) paramv[5];
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Retrieved parameters: CgEdt" << endl;
#endif

	ocrGuid_t A = depv[0].guid;
	ocrGuid_t B = depv[1].guid;
	ocrGuid_t x_old = depv[2].guid;
	ocrGuid_t p_old = depv[3].guid;
	ocrGuid_t r_old = depv[4].guid;
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Retrieved dependencies: CgEdt" << endl;
#endif

	if (k == K_ITERATIONS) {
		ocrEventSatisfy(result, x_old);
#ifdef DEBUG_MESSAGES
		cout << "k = " << k << ". Satisfied result of CgEdt" << endl;
#endif
		return NULL_GUID;
	}

	ocrGuid_t scaleEdtTemplate;
	ocrEdtTemplateCreate(&scaleEdtTemplate, scaleEdt, 3, 2);

	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 5, 2);

	ocrGuid_t transposeEdtTemplate;
	ocrEdtTemplateCreate(&transposeEdtTemplate, transposeEdt, 3, 1);

	ocrGuid_t addEdtTemplate;
	ocrEdtTemplateCreate(&addEdtTemplate, addEdt, 3, 2);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 3, 2);

	ocrGuid_t divideEdtTemplate;
	ocrEdtTemplateCreate(&divideEdtTemplate, divideEdt, 1, 2);

#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt templates: scale, product, transpose, add, subtract, divide" << endl;
#endif

	ocrGuid_t rT;
	ocrEventCreate(&rT, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtA = createEdt(transposeEdtTemplate, rT, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtA" << endl;
#endif

	ocrGuid_t rTr;
	ocrEventCreate(&rTr, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtB = createEdt(productEdtTemplate, rTr, X_old_columns,
		X_old_rows, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtB" << endl;
#endif

	ocrGuid_t pT;
	ocrEventCreate(&pT, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtC = createEdt(transposeEdtTemplate, pT, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtC" << endl;
#endif

	ocrGuid_t pTA;
	ocrEventCreate(&pTA, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtD = createEdt(productEdtTemplate, pTA, X_old_columns,
		X_old_rows, A_rows, A_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtD" << endl;
#endif

	ocrGuid_t pTAp;
	ocrEventCreate(&pTAp, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtE = createEdt(productEdtTemplate, pTAp, X_old_columns,
		X_old_rows, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtE" << endl;
#endif

	ocrGuid_t alpha;
	ocrEventCreate(&alpha, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtF = createEdt(divideEdtTemplate, alpha, 1, 1);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtF" << endl;
#endif

	//ocrDbDestroy(rT);
	//ocrDbDestroy(pT);
	//ocrDbDestroy(pTA);
	//ocrDbDestroy(pTAp);

	ocrGuid_t ap;
	ocrEventCreate(&ap, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtG = createEdt(scaleEdtTemplate, ap, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtG" << endl;
#endif

	ocrGuid_t x_new;
	ocrEventCreate(&x_new, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtH = createEdt(addEdtTemplate, x_new, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtH" << endl;
#endif

	//ocrDbDestroy(ap);

	ocrGuid_t aA;
	ocrEventCreate(&aA, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtI = createEdt(scaleEdtTemplate, aA, A_rows, A_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtI" << endl;
#endif

	ocrGuid_t aAp;
	ocrEventCreate(&aAp, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtJ = createEdt(productEdtTemplate, aAp, A_rows,
		A_columns, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtJ" << endl;
#endif

	ocrGuid_t r_new;
	ocrEventCreate(&r_new, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtK = createEdt(subtractEdtTemplate, r_new, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtK" << endl;
#endif

	//ocrDbDestroy(aA);
	//ocrDbDestroy(aAp);

	ocrGuid_t rT_new;
	ocrEventCreate(&rT_new, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtL = createEdt(transposeEdtTemplate, rT_new, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtL" << endl;
#endif

	ocrGuid_t rT_newr;
	ocrEventCreate(&rT_newr, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtM = createEdt(productEdtTemplate, rT_newr, X_old_columns,
		X_old_rows, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtM" << endl;
#endif

	ocrGuid_t beta;
	ocrEventCreate(&beta, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtN = createEdt(divideEdtTemplate, beta, rT_newr, 1, 1, rTr);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtN" << endl;
#endif

	//ocrDbDestroy(rT_new);
	//ocrDbDestroy(rT_newr);
	//ocrDbDestroy(rTr);

	ocrGuid_t bp;
	ocrEventCreate(&bp, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtO = createEdt(scaleEdtTemplate, bp, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtO" << endl;
#endif

	ocrGuid_t p_new;
	ocrEventCreate(&p_new, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtP = createEdt(addEdtTemplate, p_new, X_old_rows, X_old_columns);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtP" << endl;
#endif

	//ocrDbDestroy(bp);
	//k++;

	//ocrDbDestroy(r_old);
	//r_old = r_new;
	//		if (k != 0) {
	//			ocrDbDestroy(x_old);
	//			ocrDbDestroy(p_old);
	//		}
	//x_old = x_new;
	//p_old = p_new;

	ocrGuid_t CgEdtTemplate;
	ocrEdtTemplateCreate(&CgEdtTemplate, CgEdt, 6, 5);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt template: CgEdt" << endl;
#endif

	int z = k + 1;
	u64 nparamv[6];
	nparamv[0] = result;
	nparamv[1] = A_rows;
	nparamv[2] = A_columns;
	nparamv[3] = X_old_rows;
	nparamv[4] = X_old_columns;
	nparamv[5] = z;
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Initialized parameters: CgEdt" << endl;
#endif

	ocrGuid_t myEdt;
	//Satisfy with A, B, x_new, p_new, r_new
	ocrEdtCreate(&myEdt, CgEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: CgEdt" << endl;
#endif

	ocrAddDependence(r_old, edtA, 0, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtA" << endl;
#endif

	ocrAddDependence(rT, edtB, 0, DB_MODE_RO);
	ocrAddDependence(r_old, edtB, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtB" << endl;
#endif

	ocrAddDependence(p_old, edtC, 0, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtC" << endl;
#endif

	ocrAddDependence(pT, edtD, 0, DB_MODE_RO);
	ocrAddDependence(A, edtD, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtD" << endl;
#endif

	ocrAddDependence(pTA, edtE, 0, DB_MODE_RO);
	ocrAddDependence(p_old, edtE, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtE" << endl;
#endif

	ocrAddDependence(rTr, edtF, 0, DB_MODE_RO);
	ocrAddDependence(pTAp, edtF, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtF" << endl;
#endif

	ocrAddDependence(p_old, edtG, 0, DB_MODE_RO);
	ocrAddDependence(alpha, edtG, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtG" << endl;
#endif

	ocrAddDependence(x_old, edtH, 0, DB_MODE_RO);
	ocrAddDependence(ap, edtH, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtH" << endl;
#endif

	ocrAddDependence(A, edtI, 0, DB_MODE_RO);
	ocrAddDependence(alpha, edtI, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtI" << endl;
#endif

	ocrAddDependence(aA, edtJ, 0, DB_MODE_RO);
	ocrAddDependence(p_old, edtJ, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtJ" << endl;
#endif

	ocrAddDependence(r_old, edtK, 0, DB_MODE_RO);
	ocrAddDependence(aAp, edtK, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtK" << endl;
#endif

	ocrAddDependence(r_new, edtL, 0, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtL" << endl;
#endif

	ocrAddDependence(rT_new, edtM, 0, DB_MODE_RO);
	ocrAddDependence(r_new, edtM, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtM" << endl;
#endif

	ocrAddDependence(rT_newr, edtN, 0, DB_MODE_RO);
	ocrAddDependence(rTr, edtN, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtN" << endl;
#endif

	ocrAddDependence(p_old, edtO, 0, DB_MODE_RO);
	ocrAddDependence(beta, edtO, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtO" << endl;
#endif

	ocrAddDependence(r_new, edtP, 0, DB_MODE_RO);
	ocrAddDependence(bp, edtP, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtP" << endl;
#endif

	ocrAddDependence(A, myEdt, 0, DB_MODE_RO);
	ocrAddDependence(B, myEdt, 1, DB_MODE_RO);
	ocrAddDependence(x_new, myEdt, 2, DB_MODE_RO);
	ocrAddDependence(p_new, myEdt, 3, DB_MODE_RO);
	ocrAddDependence(r_new, myEdt, 4, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: CgEdt" << endl;
#endif

	return NULL_GUID;
}

//do error checks later

Matrix* conjugateGradient(Matrix *A, Matrix *x, Matrix *B)
{
	Matrix *x_new;
	Matrix bar(x);
	Matrix *x_old = &bar;
	Matrix *Ax = matrixProduct(A, x_old);
	Matrix *r_old = matrixSubtract(B, Ax);
	delete Ax;
	Matrix foo(r_old);
	Matrix *p_old = &foo;

	int k = 0;
	while (k < K_ITERATIONS) {
		Matrix *rT = matrixTranspose(r_old);
		Matrix *rTr = matrixProduct(rT, r_old);
		Matrix *pT = matrixTranspose(p_old);
		Matrix *pTA = matrixProduct(pT, A);
		Matrix *pTAp = matrixProduct(pTA, p_old);
		double alpha = rTr->getValue(0, 0) / pTAp->getValue(0, 0);
		delete rT;
		delete pT;
		delete pTA;
		delete pTAp;

		Matrix *ap = matrixScale(p_old, alpha);
		x_new = matrixAdd(x_old, ap);
		delete ap;

		Matrix *aA = matrixScale(A, alpha);
		Matrix *aAp = matrixProduct(aA, p_old);
		Matrix *r_new = matrixSubtract(r_old, aAp);
		delete aA;
		delete aAp;

		Matrix *rT_new = matrixTranspose(r_new);
		Matrix *rT_newr = matrixProduct(rT_new, r_new);
		double beta = rT_newr->getValue(0, 0) / rTr->getValue(0, 0);
		delete rT_new;
		delete rT_newr;
		delete rTr;

		Matrix *bp = matrixScale(p_old, beta);
		Matrix *p_new = matrixAdd(r_new, bp);
		delete bp;
		k++;

		delete r_old;
		r_old = r_new;
		//if (k != 0) {
		//delete x_old;
		//delete p_old;
		//}
		x_old = x_new;
		p_old = p_new;
	}

	return x_new;
}

//do error checks later

void conjugateGradient_OCR(Matrix *A_m, Matrix *x_m, Matrix *B_m, ocrGuid_t result)
{
	ocrGuid_t A;
	ocrGuid_t B;
	ocrGuid_t x;
#ifdef DEBUG_MESSAGES
	cout << "Init. Created guid: A, B, x" << endl;
#endif
	ocrEventCreate(&A, OCR_EVENT_STICKY_T, true);
	ocrEventCreate(&B, OCR_EVENT_STICKY_T, true);
	ocrEventCreate(&x, OCR_EVENT_STICKY_T, true);
#ifdef DEBUG_MESSAGES
	cout << "Init. Created events: A, B, x" << endl;
#endif
	ocrEventSatisfy(A, A_m->getDataBlock());
	ocrEventSatisfy(B, B_m->getDataBlock());
	ocrEventSatisfy(x, x_m->getDataBlock());
#ifdef DEBUG_MESSAGES
	cout << "Init. Satisfied events: A, B, x" << endl;
#endif

	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 5, 2);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 3, 2);

#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt templates: product, subtract" << endl;
#endif

	ocrGuid_t Ax;
	ocrEventCreate(&Ax, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtA = createEdt(productEdtTemplate, Ax, A_m->getRows(), A_m->getColumns(),
		x_m->getRows(), x_m->getColumns());
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt: Ax" << endl;
#endif

	ocrGuid_t r_old;
	ocrEventCreate(&r_old, OCR_EVENT_STICKY_T, true);
	ocrGuid_t edtB = createEdt(subtractEdtTemplate, r_old, B_m->getRows(), B_m->getColumns(),
		x_m->getRows(), x_m->getColumns());
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt: r_old" << endl;
#endif

	//ocrDbDestroy(Ax);

	int k = 0;

	//ocrGuid_t p_old = r_old;

	ocrGuid_t CgEdtTemplate;
	ocrEdtTemplateCreate(&CgEdtTemplate, CgEdt, 6, 5);
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt template: CgEdt" << endl;
#endif

	u64 nparamv[6];
	nparamv[0] = result;
	nparamv[1] = A_m->getRows();
	nparamv[2] = A_m->getColumns();
	nparamv[3] = x_m->getRows();
	nparamv[4] = x_m->getColumns();
	nparamv[5] = k;

#ifdef DEBUG_MESSAGES
	cout << "Init. Initialized parameters: CgEdt" << endl;
#endif

	ocrGuid_t myEdt;
	//Satisfy with A, B, x, p_old, r_old
	ocrEdtCreate(&myEdt, CgEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt: CgEdt" << endl;
#endif

	ocrAddDependence(A, edtA, 0, DB_MODE_RO);
	ocrAddDependence(x, edtA, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: edtA" << endl;
#endif

	ocrAddDependence(B, edtB, 0, DB_MODE_RO);
	ocrAddDependence(Ax, edtB, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: edtB" << endl;
#endif

	ocrAddDependence(A, myEdt, 0, DB_MODE_RO);
	ocrAddDependence(B, myEdt, 1, DB_MODE_RO);
	ocrAddDependence(x, myEdt, 2, DB_MODE_RO);
	ocrAddDependence(r_old, myEdt, 3, DB_MODE_RO); //p_old
	ocrAddDependence(r_old, myEdt, 4, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: CgEdt" << endl;
#endif
}

extern "C" ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	//Test values
	Matrix A(2, 2);
	A.setValue(0, 0, 4.0f);
	A.setValue(0, 1, 1.0f);
	A.setValue(1, 0, 1.0f);
	A.setValue(1, 1, 3.0f);
	Matrix x(2, 1);
	x.setValue(0, 0, 2.0f);
	x.setValue(1, 0, 1.0f);
	Matrix B(2, 1);
	B.setValue(0, 0, 1.0f);
	B.setValue(1, 0, 2.0f);

	cout << "Testing Conjugate Gradient..." << endl;
	Matrix *fishy = conjugateGradient(&A, &x, &B);
	cout << "Solution:" << endl;
	cout << fishy->getValue(0, 0) << endl;
	cout << fishy->getValue(1, 0) << endl << endl;

	cout << "Testing Conjugate Gradient OCR..." << endl;
	ocrGuid_t result;
	ocrEventCreate(&result, OCR_EVENT_STICKY_T, true);
	conjugateGradient_OCR(&A, &x, &B, result);

	u64 nparamv[1];
	nparamv[0] = x.getRows();

	ocrGuid_t printEdtTemplate;
	ocrEdtTemplateCreate(&printEdtTemplate, printEdt, 1, 1);

	ocrGuid_t printEdt;
	ocrEdtCreate(&printEdt, printEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_FINISH, NULL_GUID, NULL);
	ocrAddDependence(result, printEdt, 0, DB_MODE_RO);

	return NULL_GUID;
}

