/*
 * File:   main.cpp
 * Author: commlaptop
 *
 * Created on February 7, 2014, 1:38 PM
 */

//todo: create events for each matrix datablock. satisfy created datablocks in called Edts

#include <iostream>
#include <cmath>
#include <time.h>

extern "C" {
	#include <pthread.h>
}

#ifndef __OCR__
#define __OCR__
#endif

#include "ocr.h"

#include "matrix.h"

#define K_ITERATIONS 10

using namespace std;

void sleep( time_t delay )
{
	time_t timer0, timer1;
	time(&timer0);
	do {
		time(&timer1);
	} while ((timer1 - timer0) < delay);
}

//Parameters: A_rows, A_columns
//Dependencies: result_event, A_event, scalar_event (double)
//Output: Guid data block with result

extern "C" ocrGuid_t scaleEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t result = depv[0].guid;
	double *a = (double*) depv[1].ptr;
	double *b = (double*) depv[2].ptr;
	ocrGuid_t dataBlock;
	matrixScale(&dataBlock, a, (int) paramv[0], (int) paramv[1], b);
	ocrEventSatisfy(result, dataBlock);
	return NULL_GUID;
}

//Parameters: A_rows, A_columns, B_rows, B_columns
//Dependencies: result_event, A_event, B_event
//Output: Guid data block with result

extern "C" ocrGuid_t productEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t result = depv[0].guid;
	ocrGuid_t dataBlock;
	matrixProduct(&dataBlock, (double*) depv[1].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[2].ptr, (int) paramv[2], (int) paramv[3]);
	ocrEventSatisfy(result, dataBlock);
	return NULL_GUID;
}

//Parameters: A_rows, A_columns
//Dependencies: result_event, A_event
//Output: Guid data block with result

extern "C" ocrGuid_t transposeEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t result = depv[0].guid;
	ocrGuid_t dataBlock;
	matrixTranspose(&dataBlock, (double*) depv[1].ptr, (int) paramv[0], (int) paramv[1]);
	ocrEventSatisfy(result, dataBlock);
	return NULL_GUID;
}

//Parameters: A_rows, A_columns
//Dependencies: result_event, A_event, B_event (A and B need same rows/columns)
//Output: Guid data block with result

extern "C" ocrGuid_t addEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t result = depv[0].guid;
	ocrGuid_t dataBlock;
	matrixAdd(&dataBlock, (double*) depv[1].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[2].ptr);
	ocrEventSatisfy(result, dataBlock);
	return NULL_GUID;
}

//Parameters: A_rows, A_columns
//Dependencies: result_event, A_event, B_event (A and B need same rows/columns)
//Output: Guid data block with result

extern "C" ocrGuid_t subtractEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t result = depv[0].guid;
	ocrGuid_t dataBlock;
	matrixSubtract(&dataBlock, (double*) depv[1].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[2].ptr);
	ocrEventSatisfy(result, dataBlock);
	return NULL_GUID;
}

//Parameters: None
//Dependencies: result_event, A_event (double), B_event (double)
//Output: Guid data block with result (double)

extern "C" ocrGuid_t divideEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t result = depv[0].guid;
	double* r;
	ocrGuid_t dataBlock;
	DBCREATE(&dataBlock, (void**) &r, 1 * sizeof (double), DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	double* a = (double*) depv[1].ptr;
	double* b = (double*) depv[2].ptr;
	r[0] = a[0] / b[0];
	ocrEventSatisfy(result, dataBlock);
	return NULL_GUID;
}

//Parameters: x_rows
//Dependencies: A_event
//Output: NULL_GUID
//Not thread safe
extern "C" ocrGuid_t printEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	int rows = (int)paramv[0];
	double* r = (double*)depv[0].ptr;

	cout << "Solution:" << endl;
	for (int i = 0; i < rows; i++)
		cout << i << ",1: " << r[i] << endl;

	cout << "OCR Shutdown..." << endl;
	ocrShutdown();
	return NULL_GUID;
}

ocrGuid_t createEdt(ocrGuid_t templateGuid, ocrGuid_t &outputGuid,
	int A_rows = 0, int A_columns = 0, int B_rows = 0, int B_columns = 0)
{
	u64 nparamv[4];
	nparamv[0] = (u64) A_rows;
	nparamv[1] = (u64) A_columns;
	nparamv[2] = (u64) B_rows;
	nparamv[3] = (u64) B_columns;
	ocrGuid_t myEdt;
	ocrEdtCreate(&myEdt, templateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL_GUID, EDT_PROP_NONE, NULL_GUID, &outputGuid);
	ocrAddDependence(outputGuid, myEdt, 0, DB_MODE_ITW);
	return myEdt;
}

//Parameters: A_rows, A_columns, X_old_rows, X_old_columns, k, doneEvent
//Dependencies: A_db, B_db, X_old_db, P_old_db, R_old_db
//Output: Guid data block with result of x_new after K_ITERATIONS

extern "C" ocrGuid_t CgEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	int A_rows = (int) paramv[0];
	int A_columns = (int) paramv[1];
	int X_old_rows = (int) paramv[2];
	int X_old_columns = (int) paramv[3];
	int k = (int) paramv[4];

	ocrGuid_t result = depv[0].guid;
	ocrGuid_t A = depv[1].guid;
	ocrGuid_t B = depv[2].guid;
	ocrGuid_t x_old = depv[3].guid;
	ocrGuid_t p_old = depv[4].guid;
	ocrGuid_t r_old = depv[5].guid;
	
	cout << "k = " << k << endl;
	
	if (k == K_ITERATIONS) {
		ocrEventSatisfy(result, x_old);
		return NULL_GUID;
	}

	ocrGuid_t scaleEdtTemplate;
	ocrEdtTemplateCreate(&scaleEdtTemplate, scaleEdt, 2, 3);

	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 4, 3);

	ocrGuid_t transposeEdtTemplate;
	ocrEdtTemplateCreate(&transposeEdtTemplate, transposeEdt, 2, 2);

	ocrGuid_t addEdtTemplate;
	ocrEdtTemplateCreate(&addEdtTemplate, addEdt, 2, 3);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 2, 3);

	ocrGuid_t divideEdtTemplate;
	ocrEdtTemplateCreate(&divideEdtTemplate, divideEdt, 0, 3);

	ocrGuid_t rT;
	ocrEventCreate(&rT, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtA = createEdt(transposeEdtTemplate, rT, X_old_rows, X_old_columns);

	ocrGuid_t rTr;
	ocrEventCreate(&rTr, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtB = createEdt(productEdtTemplate, rTr, X_old_columns,
		X_old_rows, X_old_rows, X_old_columns);

	ocrGuid_t pT;
	ocrEventCreate(&pT, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtC = createEdt(transposeEdtTemplate, pT, X_old_rows, X_old_columns);

	ocrGuid_t pTA;
	ocrEventCreate(&pTA, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtD = createEdt(productEdtTemplate, pTA, X_old_columns,
		X_old_rows, A_rows, A_columns);

	ocrGuid_t pTAp;
	ocrEventCreate(&pTAp, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtE = createEdt(productEdtTemplate, pTAp, X_old_columns,
		X_old_rows, X_old_rows, X_old_columns);

	ocrGuid_t alpha;
	ocrEventCreate(&alpha, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtF = createEdt(divideEdtTemplate, alpha, 1, 1);

	//ocrDbDestroy(rT);
	//ocrDbDestroy(pT);
	//ocrDbDestroy(pTA);
	//ocrDbDestroy(pTAp);

	ocrGuid_t ap;
	ocrEventCreate(&ap, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtG = createEdt(scaleEdtTemplate, ap, X_old_rows, X_old_columns);

	ocrGuid_t x_new;
	ocrEventCreate(&x_new, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtH = createEdt(addEdtTemplate, x_new, X_old_rows, X_old_columns);

	//ocrDbDestroy(ap);

	ocrGuid_t aA;
	ocrEventCreate(&aA, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtI = createEdt(scaleEdtTemplate, aA, A_rows, A_columns);

	ocrGuid_t aAp;
	ocrEventCreate(&aAp, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtJ = createEdt(productEdtTemplate, aAp, A_rows,
		A_columns, X_old_rows, X_old_columns);

	ocrGuid_t r_new;
	ocrEventCreate(&r_new, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtK = createEdt(subtractEdtTemplate, r_new, X_old_rows, X_old_columns);

	//ocrDbDestroy(aA);
	//ocrDbDestroy(aAp);
	
	ocrGuid_t rT_new;
	ocrEventCreate(&rT_new, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtL = createEdt(transposeEdtTemplate, rT_new, X_old_rows, X_old_columns);

	ocrGuid_t rT_newr;
	ocrEventCreate(&rT_newr, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtM = createEdt(productEdtTemplate, rT_newr, X_old_columns,
		X_old_rows, X_old_rows, X_old_columns);

	ocrGuid_t beta;
	ocrEventCreate(&beta, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtN = createEdt(divideEdtTemplate, beta, rT_newr, 1, 1, rTr);

	//ocrDbDestroy(rT_new);
	//ocrDbDestroy(rT_newr);
	//ocrDbDestroy(rTr);
	
	ocrGuid_t bp;
	ocrEventCreate(&bp, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtO = createEdt(scaleEdtTemplate, bp, X_old_rows, X_old_columns);

	ocrGuid_t p_new;
	ocrEventCreate(&p_new, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtP = createEdt(addEdtTemplate, p_new, X_old_rows, X_old_columns);

	//ocrDbDestroy(bp);
	k++;

	//ocrDbDestroy(r_old);
	//r_old = r_new;
	//		if (k != 0) {
	//			ocrDbDestroy(x_old);
	//			ocrDbDestroy(p_old);
	//		}
	//x_old = x_new;
	//p_old = p_new;

	ocrGuid_t CgEdtTemplate;
	ocrEdtTemplateCreate(&CgEdtTemplate, CgEdt, 5, 6);

	u64 nparamv[5];
	nparamv[0] = (u64) A_rows;
	nparamv[1] = (u64) A_columns;
	nparamv[2] = (u64) X_old_rows;
	nparamv[3] = (u64) X_old_columns;
	nparamv[4] = (u64) k;

	ocrGuid_t myEdt;
	//Satisfy with A, B, x_new, p_new, r_new
	ocrEdtCreate(&myEdt, CgEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
	
	ocrAddDependence(r_old, edtA, 1, DB_MODE_RO);

	ocrAddDependence(rT, edtB, 1, DB_MODE_RO);
	ocrAddDependence(r_old, edtB, 2, DB_MODE_RO);

	ocrAddDependence(p_old, edtC, 1, DB_MODE_RO);

	ocrAddDependence(pT, edtD, 1, DB_MODE_RO);
	ocrAddDependence(A, edtD, 2, DB_MODE_RO);

	ocrAddDependence(pTA, edtE, 1, DB_MODE_RO);
	ocrAddDependence(p_old, edtE, 2, DB_MODE_RO);

	ocrAddDependence(rTr, edtF, 1, DB_MODE_RO);
	ocrAddDependence(pTAp, edtF, 2, DB_MODE_RO);

	ocrAddDependence(p_old, edtG, 1, DB_MODE_RO);
	ocrAddDependence(alpha, edtG, 2, DB_MODE_RO);

	ocrAddDependence(x_old, edtH, 1, DB_MODE_RO);
	ocrAddDependence(ap, edtH, 2, DB_MODE_RO);

	ocrAddDependence(A, edtI, 1, DB_MODE_RO);
	ocrAddDependence(alpha, edtI, 2, DB_MODE_RO);

	ocrAddDependence(aA, edtJ, 1, DB_MODE_RO);
	ocrAddDependence(p_old, edtJ, 2, DB_MODE_RO);

	ocrAddDependence(r_old, edtK, 1, DB_MODE_RO);
	ocrAddDependence(aAp, edtK, 2, DB_MODE_RO);
	
	ocrAddDependence(r_new, edtL, 1, DB_MODE_RO);

	ocrAddDependence(rT_new, edtM, 1, DB_MODE_RO);
	ocrAddDependence(r_new, edtM, 2, DB_MODE_RO);

	ocrAddDependence(rT_newr, edtN, 1, DB_MODE_RO);
	ocrAddDependence(rTr, edtN, 2, DB_MODE_RO);

	ocrAddDependence(p_old, edtO, 1, DB_MODE_RO);
	ocrAddDependence(beta, edtO, 2, DB_MODE_RO);

	ocrAddDependence(r_new, edtP, 1, DB_MODE_RO);
	ocrAddDependence(bp, edtP, 2, DB_MODE_RO);
	
	ocrAddDependence(result, myEdt, 0, DB_MODE_ITW);
	ocrAddDependence(A, myEdt, 1, DB_MODE_RO);
	ocrAddDependence(B, myEdt, 2, DB_MODE_RO);
	ocrAddDependence(x_new, myEdt, 3, DB_MODE_RO);
	ocrAddDependence(p_new, myEdt, 4, DB_MODE_RO);
	ocrAddDependence(r_new, myEdt, 5, DB_MODE_RO);

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

void conjugateGradient_OCR(Matrix *A, Matrix *x, Matrix *B, ocrGuid_t result)
{
	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 4, 3);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 2, 3);

	ocrGuid_t Ax;
	ocrEventCreate(&Ax, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtA = createEdt(productEdtTemplate, Ax, A->getRows(), A->getColumns(), 
		x->getRows(), x->getColumns());

	ocrGuid_t r_old;
	ocrEventCreate(&r_old, OCR_EVENT_IDEM_T, true);
	ocrGuid_t edtB = createEdt(subtractEdtTemplate, r_old, B->getRows(), B->getColumns(), 
		x->getRows(), x->getColumns());

	//ocrDbDestroy(Ax);

	int k = 0;

	//ocrGuid_t p_old = r_old;

	ocrGuid_t CgEdtTemplate;
	ocrEdtTemplateCreate(&CgEdtTemplate, CgEdt, 5, 6);

	u64 nparamv[5];
	nparamv[0] = (u64) A->getRows();
	nparamv[1] = (u64) A->getColumns();
	nparamv[2] = (u64) x->getRows();
	nparamv[3] = (u64) x->getColumns();
	nparamv[4] = (u64) k;

	ocrGuid_t myEdt;
	//Satisfy with A, B, x, p_old, r_old
	ocrEdtCreate(&myEdt, CgEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
	
	ocrAddDependence(A->getDataBlock(), edtA, 1, DB_MODE_RO);
	ocrAddDependence(x->getDataBlock(), edtA, 2, DB_MODE_RO);

	ocrAddDependence(B->getDataBlock(), edtB, 1, DB_MODE_RO);
	ocrAddDependence(Ax, edtB, 1, DB_MODE_RO);
	
	ocrAddDependence(result, myEdt, 0, DB_MODE_ITW);
	ocrAddDependence(A->getDataBlock(), myEdt, 1, DB_MODE_RO);
	ocrAddDependence(B->getDataBlock(), myEdt, 2, DB_MODE_RO);
	ocrAddDependence(x->getDataBlock(), myEdt, 3, DB_MODE_RO);
	ocrAddDependence(r_old, myEdt, 4, DB_MODE_RO); //p_old
	ocrAddDependence(r_old, myEdt, 5, DB_MODE_RO);
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
	ocrEventCreate(&result, OCR_EVENT_IDEM_T, true);
	conjugateGradient_OCR(&A, &x, &B, result);
	
	u64 nparamv[1];
	nparamv[0] = (u64)x.getRows();

	ocrGuid_t printEdtTemplate;
	ocrEdtTemplateCreate(&printEdtTemplate, printEdt, 1, 1);

	ocrGuid_t printEdt;
	ocrEdtCreate(&printEdt, printEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_FINISH, NULL_GUID, NULL);
	ocrAddDependence(result, printEdt, 0, DB_MODE_RO);

	return NULL_GUID;
}

