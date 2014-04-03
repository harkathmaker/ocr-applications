/*
 * File:   main.cpp
 * Author: commlaptop
 *
 * Created on February 7, 2014, 1:38 PM
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <cstdlib>

extern "C" {
#include <pthread.h>
}

//#define DEBUG_MESSAGES

#ifndef __OCR__
#define __OCR__
#endif

#include "ocr.h"

#include "matrix.h"
#include "timer.h"

#define K_ITERATIONS 100
#define MATRIX_N 500

using namespace std;

// Generates a random double
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

//Parameters: A_rows, A_columns
//Dependencies: A_event, scalar_event (double)
//Output: result

extern "C" ocrGuid_t scaleEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	double *a = (double*) depv[0].ptr;
	double *b = (double*) depv[1].ptr;
	ocrGuid_t dataBlock = matrixScale(a, (int) paramv[0], (int) paramv[1], b);
#ifdef DEBUG_MESSAGES
	cout << "scaleEdt()" << endl;
#endif
	return dataBlock;
}

//Parameters: A_rows, A_columns, B_rows, B_columns
//Dependencies: A_event, B_event
//Output: result

extern "C" ocrGuid_t productEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixProduct((double*) depv[0].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[1].ptr, (int) paramv[2], (int) paramv[3]);
#ifdef DEBUG_MESSAGES
	cout << "productEdt()" << endl;
#endif
	return dataBlock;
}

//Parameters: A_rows, A_columns
//Dependencies: A_event
//Output: result

extern "C" ocrGuid_t transposeEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixTranspose((double*) depv[0].ptr, (int) paramv[0], (int) paramv[1]);
#ifdef DEBUG_MESSAGES
	cout << "transposeEdt()" << endl;
#endif
	return dataBlock;
}

//Parameters: A_rows, A_columns
//Dependencies: A_event, B_event (A and B need same rows/columns)
//Output: result

extern "C" ocrGuid_t addEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixAdd((double*) depv[0].ptr, (int) paramv[0], (int) paramv[1], (double*) depv[1].ptr);
#ifdef DEBUG_MESSAGES
	cout << "addEdt()" << endl;
#endif
	return dataBlock;
}

//Parameters: A_rows, A_columns
//Dependencies: A_event, B_event (A and B need same rows/columns)
//Output: result

extern "C" ocrGuid_t subtractEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixSubtract((double*) depv[0].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[1].ptr);
#ifdef DEBUG_MESSAGES
	cout << "subtractEdt()" << endl;
#endif
	return dataBlock;
}

//Parameters: None
//Dependencies: A_event (double), B_event (double)
//Output: result

extern "C" ocrGuid_t divideEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	double* r;
	ocrGuid_t dataBlock;
	DBCREATE(&dataBlock, (void**) &r, 1 * sizeof (double), DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	double* a = (double*) depv[0].ptr;
	double* b = (double*) depv[1].ptr;
	r[0] = a[0] / b[0];
#ifdef DEBUG_MESSAGES
	cout << "divideEdt()" << endl;
#endif
	return dataBlock;
}

//Parameters: x_rows, t1.sec, t1.usec
//Dependencies: A_event
//Output: NULL_GUID
//Not thread safe

extern "C" ocrGuid_t printEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "printEdt()" << endl;
#endif
	int rows = (int) paramv[0];
	time_t sec = (time_t) paramv[1];
	time_t usec = (suseconds_t) paramv[2];
	double* r = (double*) depv[0].ptr;

	cout << "Solution:" << endl;
	for (int i = 0; i < rows; i++)
		cout << r[i] << " ";
	cout << endl;

	timeval t2;
	t2.tv_sec = sec;
	t2.tv_usec = usec;
	double cgOcrTimeElapsed = tock(t2);

	cout << "Time to test Conjugate Gradient OCR: " << cgOcrTimeElapsed << " ms" << endl << endl;

#ifdef DEBUG_MESSAGES
	cout << "OCR Shutdown..." << endl;
#endif
	ocrShutdown();
	return NULL_GUID;
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

	u64 edtParams[4];

	ocrGuid_t scaleEdtTemplate;
	ocrEdtTemplateCreate(&scaleEdtTemplate, scaleEdt, 2, 2);

	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 4, 2);

	ocrGuid_t transposeEdtTemplate;
	ocrEdtTemplateCreate(&transposeEdtTemplate, transposeEdt, 2, 1);

	ocrGuid_t addEdtTemplate;
	ocrEdtTemplateCreate(&addEdtTemplate, addEdt, 2, 2);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 2, 2);

	ocrGuid_t divideEdtTemplate;
	ocrEdtTemplateCreate(&divideEdtTemplate, divideEdt, 0, 2);

#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt templates: scale, product, transpose, add, subtract, divide" << endl;
#endif

	ocrGuid_t rT;
	ocrGuid_t edtA;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtA, transposeEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &rT);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtA" << endl;
#endif

	ocrGuid_t rTr;
	ocrGuid_t edtB;
	edtParams[0] = X_old_columns;
	edtParams[1] = X_old_rows;
	edtParams[2] = X_old_rows;
	edtParams[3] = X_old_columns;
	ocrEdtCreate(&edtB, productEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &rTr);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtB" << endl;
#endif

	ocrGuid_t pT;
	ocrGuid_t edtC;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtC, transposeEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &pT);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtC" << endl;
#endif

	ocrGuid_t pTA;
	ocrGuid_t edtD;
	edtParams[0] = X_old_columns;
	edtParams[1] = X_old_rows;
	edtParams[2] = A_rows;
	edtParams[3] = A_columns;
	ocrEdtCreate(&edtD, productEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &pTA);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtD" << endl;
#endif

	ocrGuid_t pTAp;
	ocrGuid_t edtE;
	edtParams[0] = X_old_columns;
	edtParams[1] = X_old_rows;
	edtParams[2] = X_old_rows;
	edtParams[3] = X_old_columns;
	ocrEdtCreate(&edtE, productEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &pTAp);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtE" << endl;
#endif

	ocrGuid_t alpha;
	ocrGuid_t edtF;
	ocrEdtCreate(&edtF, divideEdtTemplate, EDT_PARAM_DEF, NULL, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &alpha);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtF" << endl;
#endif

	//ocrDbDestroy(rT);
	//ocrDbDestroy(pT);
	//ocrDbDestroy(pTA);
	//ocrDbDestroy(pTAp);

	ocrGuid_t ap;
	ocrGuid_t edtG;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtG, scaleEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &ap);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtG" << endl;
#endif

	ocrGuid_t x_new;
	ocrGuid_t edtH;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtH, addEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &x_new);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtH" << endl;
#endif

	//ocrDbDestroy(ap);

	ocrGuid_t Ap;
	ocrGuid_t edtI;
	edtParams[0] = A_rows;
	edtParams[1] = A_columns;
	edtParams[2] = X_old_rows;
	edtParams[3] = X_old_columns;
	ocrEdtCreate(&edtI, productEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &Ap);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtI" << endl;
#endif

	ocrGuid_t aAp;
	ocrGuid_t edtJ;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtJ, scaleEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &aAp);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtJ" << endl;
#endif

	ocrGuid_t r_new;
	ocrGuid_t edtK;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtK, subtractEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &r_new);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtK" << endl;
#endif

	//ocrDbDestroy(aA);
	//ocrDbDestroy(aAp);

	ocrGuid_t rT_new;
	ocrGuid_t edtL;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtL, transposeEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &rT_new);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtL" << endl;
#endif

	ocrGuid_t rT_newr;
	ocrGuid_t edtM;
	edtParams[0] = X_old_columns;
	edtParams[1] = X_old_rows;
	edtParams[2] = X_old_rows;
	edtParams[3] = X_old_columns;
	ocrEdtCreate(&edtM, productEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &rT_newr);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtM" << endl;
#endif

	ocrGuid_t beta;
	ocrGuid_t edtN;
	ocrEdtCreate(&edtN, divideEdtTemplate, EDT_PARAM_DEF, NULL, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &beta);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtN" << endl;
#endif

	//ocrDbDestroy(rT_new);
	//ocrDbDestroy(rT_newr);
	//ocrDbDestroy(rTr);

	ocrGuid_t bp;
	ocrGuid_t edtO;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtO, scaleEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &bp);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtO" << endl;
#endif

	ocrGuid_t p_new;
	ocrGuid_t edtP;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtP, addEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &p_new);
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

	ocrAddDependence(A, myEdt, 0, DB_MODE_RO);
	ocrAddDependence(B, myEdt, 1, DB_MODE_RO);
	ocrAddDependence(x_new, myEdt, 2, DB_MODE_RO);
	ocrAddDependence(p_new, myEdt, 3, DB_MODE_RO);
	ocrAddDependence(r_new, myEdt, 4, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: CgEdt" << endl;
#endif

	ocrAddDependence(r_new, edtP, 0, DB_MODE_RO);
	ocrAddDependence(bp, edtP, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtP" << endl;
#endif

	ocrAddDependence(p_old, edtO, 0, DB_MODE_RO);
	ocrAddDependence(beta, edtO, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtO" << endl;
#endif

	ocrAddDependence(rT_newr, edtN, 0, DB_MODE_RO);
	ocrAddDependence(rTr, edtN, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtN" << endl;
#endif

	ocrAddDependence(rT_new, edtM, 0, DB_MODE_RO);
	ocrAddDependence(r_new, edtM, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtM" << endl;
#endif

	ocrAddDependence(r_new, edtL, 0, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtL" << endl;
#endif

	ocrAddDependence(r_old, edtK, 0, DB_MODE_RO);
	ocrAddDependence(aAp, edtK, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtK" << endl;
#endif

	ocrAddDependence(Ap, edtJ, 0, DB_MODE_RO);
	ocrAddDependence(alpha, edtJ, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtJ" << endl;
#endif

	ocrAddDependence(A, edtI, 0, DB_MODE_RO);
	ocrAddDependence(p_old, edtI, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtI" << endl;
#endif

	ocrAddDependence(x_old, edtH, 0, DB_MODE_RO);
	ocrAddDependence(ap, edtH, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtH" << endl;
#endif

	ocrAddDependence(p_old, edtG, 0, DB_MODE_RO);
	ocrAddDependence(alpha, edtG, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtG" << endl;
#endif

	ocrAddDependence(rTr, edtF, 0, DB_MODE_RO);
	ocrAddDependence(pTAp, edtF, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtF" << endl;
#endif

	ocrAddDependence(pTA, edtE, 0, DB_MODE_RO);
	ocrAddDependence(p_old, edtE, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtE" << endl;
#endif

	ocrAddDependence(pT, edtD, 0, DB_MODE_RO);
	ocrAddDependence(A, edtD, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtD" << endl;
#endif

	ocrAddDependence(p_old, edtC, 0, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtC" << endl;
#endif

	ocrAddDependence(rT, edtB, 0, DB_MODE_RO);
	ocrAddDependence(r_old, edtB, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtB" << endl;
#endif

	ocrAddDependence(r_old, edtA, 0, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtA" << endl;
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

	cout << "derp" << endl;

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
		if (k >= 2) {
			if (x_old)
				delete x_old;
			if (p_old)
				delete p_old;
		}
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
	ocrEventCreate(&A, OCR_EVENT_STICKY_T, true);
	ocrEventCreate(&B, OCR_EVENT_STICKY_T, true);
	ocrEventCreate(&x, OCR_EVENT_STICKY_T, true);
	ocrEventSatisfy(A, A_m->getDataBlock());
	ocrEventSatisfy(B, B_m->getDataBlock());
	ocrEventSatisfy(x, x_m->getDataBlock());
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Guids: A, B, x" << endl;
#endif

	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 4, 2);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 2, 2);

#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt templates: product, subtract" << endl;
#endif

	u64 edtParams[4];

	ocrGuid_t Ax;
	ocrGuid_t edtA;
	edtParams[0] = A_m->getRows();
	edtParams[1] = A_m->getColumns();
	edtParams[2] = x_m->getRows();
	edtParams[3] = x_m->getColumns();
	ocrEdtCreate(&edtA, productEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &Ax);
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt: Ax" << endl;
#endif

	ocrGuid_t r_old;
	ocrGuid_t edtB;
	edtParams[0] = B_m->getRows();
	edtParams[1] = B_m->getColumns();
	edtParams[2] = x_m->getRows();
	edtParams[3] = x_m->getColumns();
	ocrEdtCreate(&edtB, subtractEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &r_old);
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

	ocrAddDependence(A, myEdt, 0, DB_MODE_RO);
	ocrAddDependence(B, myEdt, 1, DB_MODE_RO);
	ocrAddDependence(x, myEdt, 2, DB_MODE_RO);
	ocrAddDependence(r_old, myEdt, 3, DB_MODE_RO); //p_old
	ocrAddDependence(r_old, myEdt, 4, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: CgEdt" << endl;
#endif

	ocrAddDependence(B, edtB, 0, DB_MODE_RO);
	ocrAddDependence(Ax, edtB, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: edtB" << endl;
#endif

	ocrAddDependence(A, edtA, 0, DB_MODE_RO);
	ocrAddDependence(x, edtA, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: edtA" << endl;
#endif
}

void readMatrixFromFile(const char* fileName, Matrix &A)
{
	ifstream f;
	f.open(fileName);
	double n;
	for (int row = 0; row < MATRIX_N; row++) {
		for (int column = 0; column < MATRIX_N; column++) {
			f >> n;
			A.setValue(row, column, n);
		}
	}
	f.close();
}

extern "C" ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	Matrix A(MATRIX_N, MATRIX_N);
	readMatrixFromFile("cg.mat", A);
	Matrix x(MATRIX_N, 1);
	Matrix B(MATRIX_N, 1);
	for (int row = 0; row < MATRIX_N; row++) {
		x.setValue(row, 0, 0.0f); //Initial guess of all zeroes
		B.setValue(row, 0, fRand(-10.0f, 10.0f));
	}

	//timeval cgTime;
	//tick(cgTime);
	//cout << "Testing Conjugate Gradient..." << endl;
	//Matrix *fishy = conjugateGradient(&A, &x, &B);

	//cout << "Solution:" << endl;
	//for (int i = 0; i < MATRIX_N; i++)
	//	cout << fishy->getValue(i, 0) << " ";
	//cout << endl;

	//double cgTimeElapsed = tock(cgTime);

	//cout << "Time to test Conjugate Gradient: " << cgTimeElapsed << " ms" << endl << endl;

	timeval ocrCgTime;
	tick(ocrCgTime);
	cout << "Testing Conjugate Gradient OCR..." << endl;
	ocrGuid_t result;
	ocrEventCreate(&result, OCR_EVENT_STICKY_T, true);
	conjugateGradient_OCR(&A, &x, &B, result);

	u64 nparamv[3];
	nparamv[0] = x.getRows();
	nparamv[1] = ocrCgTime.tv_sec;
	nparamv[2] = ocrCgTime.tv_usec;

	ocrGuid_t printEdtTemplate;
	ocrEdtTemplateCreate(&printEdtTemplate, printEdt, 3, 1);

	ocrGuid_t printEdt;
	ocrEdtCreate(&printEdt, printEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_FINISH, NULL_GUID, NULL);
	ocrAddDependence(result, printEdt, 0, DB_MODE_RO);

	return NULL_GUID;
}

