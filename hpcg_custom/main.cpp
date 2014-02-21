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

#define RESIDUAL_LIMIT 1.0e-9

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
//Dependencies: A_event, scalar_event (double), output_event
//Output: Guid event satisfied with result

extern "C" ocrGuid_t scaleEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t e = depv[2].guid;
	ocrGuid_t result;
	double *a = (double*) depv[0].ptr;
	double *b = (double*) depv[1].ptr;
	matrixScale(&result, a, (int) paramv[0], (int) paramv[1], b);
	ocrEventSatisfy(e, result);
	return e;
}

//Parameters: A_rows, A_columns, B_rows, B_columns
//Dependencies: A_event, B_event, output_event
//Output: Guid event satisfied with result

extern "C" ocrGuid_t productEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t e = depv[2].guid;
	ocrGuid_t result;
	matrixProduct(&result, (double*) depv[0].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[1].ptr, (int) paramv[2], (int) paramv[3]);
	ocrEventSatisfy(e, result);
	return e;
}

//Parameters: A_rows, A_columns
//Dependencies: A_event, output_event
//Output: Guid event satisfied with result

extern "C" ocrGuid_t transposeEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t e = depv[1].guid;
	ocrGuid_t result;
	matrixTranspose(&result, (double*) depv[0].ptr, (int) paramv[0], (int) paramv[1]);
	ocrEventSatisfy(e, result);
	return e;
}

//Parameters: A_rows, A_columns
//Dependencies: A_event, B_event (A and B need same rows/columns), output_event
//Output: Guid event satisfied with result

extern "C" ocrGuid_t addEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t e = depv[2].guid;
	ocrGuid_t result;
	matrixAdd(&result, (double*) depv[0].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[1].ptr);
	ocrEventSatisfy(e, result);
	return e;
}

//Parameters: A_rows, A_columns
//Dependencies: A_event, B_event (A and B need same rows/columns), output_event
//Output: Guid event satisfied with result

extern "C" ocrGuid_t subtractEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t e = depv[2].guid;
	ocrGuid_t result;
	matrixSubtract(&result, (double*) depv[0].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[1].ptr);
	ocrEventSatisfy(e, result);
	return e;
}

//Parameters: None
//Dependencies: A_event (double), B_event (double), output_event
//Output: Guid event satisfied with result (double)

extern "C" ocrGuid_t divideEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t e = depv[2].guid;
	ocrGuid_t result;
	double* r;
	DBCREATE(&result, (void**) &r, 1 * sizeof (double), DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	double* a = (double*) depv[0].ptr;
	double* b = (double*) depv[1].ptr;
	if (b[0] != 0.000000f) {
		r[0] = a[0] / b[0];
		ocrEventSatisfy(e, result);
		return e;
	} else {
		r[0] = 0.0f;
		ocrEventSatisfy(e, result);
		return -1;
	}
}

//Parameters: x_rows
//Dependencies: A_event, output_event
//Output: Guid event satisfied
//Not thread safe
extern "C" ocrGuid_t printEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t e = depv[2].guid;
	int rows = (int)paramv[0];
	double* r = (double*)depv[0].ptr;

	for (int i = 0; i < rows; i++)
		cout << i << ",1: " << r[i] << endl;
	
	ocrEventSatisfy(e, NULL_GUID);
	return e;
}

void createEdt(ocrGuid_t templateGuid, ocrGuid_t &outputGuid, ocrGuid_t A,
	int A_rows, int A_columns, ocrGuid_t B = -1, int B_rows = -1,
	int B_columns = -1)
{
	u64 nparamv[4];
	nparamv[0] = (u64) A_rows;
	nparamv[1] = (u64) A_columns;
	nparamv[2] = (u64) B_rows;
	nparamv[3] = (u64) B_columns;
	ocrGuid_t myEdt;
	ocrEdtCreate(&myEdt, templateGuid, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &outputGuid);
	ocrAddDependence(A, myEdt, 0, DB_MODE_RO);
	if (B != -1)
		ocrAddDependence(B, myEdt, 1, DB_MODE_RO);
}

extern "C" ocrGuid_t CgEdt_StepB(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[]);

//Parameters: A_rows, A_columns, X_old_rows, X_old_columns, k
//Dependencies: A_db, B_db, X_old_db, P_old_db, R_old_db, output_event
//Output: Guid event satisfied with result of conjugate gradient

extern "C" ocrGuid_t CgEdt_StepA(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	int A_rows = (int) paramv[0];
	int A_columns = (int) paramv[1];
	int X_old_rows = (int) paramv[2];
	int X_old_columns = (int) paramv[3];
	int k = (int) paramv[4];

	ocrGuid_t A = depv[0].guid;
	ocrGuid_t B = depv[1].guid;
	ocrGuid_t x_old = depv[2].guid;
	ocrGuid_t p_old = depv[3].guid;
	ocrGuid_t r_old = depv[4].guid;
	ocrGuid_t e = depv[5].guid;

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
	createEdt(transposeEdtTemplate, rT, r_old, X_old_rows, X_old_columns);

	ocrGuid_t rTr;
	ocrEventCreate(&rTr, OCR_EVENT_IDEM_T, true);
	createEdt(productEdtTemplate, rTr, rT, X_old_columns,
		X_old_rows, r_old, X_old_rows, X_old_columns);

	ocrGuid_t pT;
	ocrEventCreate(&pT, OCR_EVENT_IDEM_T, true);
	createEdt(transposeEdtTemplate, pT, p_old, X_old_rows, X_old_columns);

	ocrGuid_t pTA;
	ocrEventCreate(&pTA, OCR_EVENT_IDEM_T, true);
	createEdt(productEdtTemplate, pTA, pT, X_old_columns,
		X_old_rows, A, A_rows, A_columns);

	ocrGuid_t pTAp;
	ocrEventCreate(&pTAp, OCR_EVENT_IDEM_T, true);
	createEdt(productEdtTemplate, pTAp, pTA, X_old_columns,
		X_old_rows, p_old, X_old_rows, X_old_columns);

	ocrGuid_t alpha;
	ocrEventCreate(&alpha, OCR_EVENT_IDEM_T, true);
	createEdt(divideEdtTemplate, alpha, rTr, 1, 1, pTAp);

	//ocrDbDestroy(rT);
	//ocrDbDestroy(pT);
	//ocrDbDestroy(pTA);
	//ocrDbDestroy(pTAp);

	ocrGuid_t ap;
	ocrEventCreate(&ap, OCR_EVENT_IDEM_T, true);
	createEdt(scaleEdtTemplate, ap, p_old, X_old_rows, X_old_columns, alpha);

	ocrGuid_t x_new;
	ocrEventCreate(&x_new, OCR_EVENT_IDEM_T, true);
	createEdt(addEdtTemplate, x_new, x_old, X_old_rows, X_old_columns, ap);

	//ocrDbDestroy(ap);

	ocrGuid_t aA;
	ocrEventCreate(&aA, OCR_EVENT_IDEM_T, true);
	createEdt(scaleEdtTemplate, aA, A, A_rows, A_columns, alpha);

	ocrGuid_t aAp;
	ocrEventCreate(&aAp, OCR_EVENT_IDEM_T, true);
	createEdt(productEdtTemplate, aAp, aA, A_rows,
		A_columns, p_old, X_old_rows, X_old_columns);

	ocrGuid_t r_new;
	ocrEventCreate(&r_new, OCR_EVENT_IDEM_T, true);
	createEdt(subtractEdtTemplate, r_new, r_old, X_old_rows, X_old_columns, aAp);

	//ocrDbDestroy(aA);
	//ocrDbDestroy(aAp);

	ocrGuid_t stepBEdtTemplate;
	ocrEdtTemplateCreate(&stepBEdtTemplate, CgEdt_StepB, 5, 7);

	u64 nparamv[5];
	nparamv[0] = (u64) A_rows;
	nparamv[1] = (u64) A_columns;
	nparamv[2] = (u64) X_old_rows;
	nparamv[3] = (u64) X_old_columns;
	nparamv[4] = (u64) k;

	ocrGuid_t stepBEdt;
	ocrEdtCreate(&stepBEdt, stepBEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &e);
	ocrAddDependence(A, stepBEdt, 0, DB_MODE_RO);
	ocrAddDependence(B, stepBEdt, 1, DB_MODE_RO);
	ocrAddDependence(x_new, stepBEdt, 2, DB_MODE_RO);
	ocrAddDependence(p_old, stepBEdt, 3, DB_MODE_RO);
	ocrAddDependence(r_new, stepBEdt, 4, DB_MODE_RO);
	ocrAddDependence(rTr, stepBEdt, 5, DB_MODE_RO);
	ocrAddDependence(e, stepBEdt, 6, DB_MODE_RO);

	return e;
}

//Parameters: A_rows, A_columns, X_new_rows, X_new_columns, k
//Dependencies: A_db, B_db, x_new_db, p_old_db, r_new_db, rTr_db, output_event
//Output: Guid event satisfied with result of conjugate gradient

extern "C" ocrGuid_t CgEdt_StepB(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	int A_rows = (int) paramv[0];
	int A_columns = (int) paramv[1];
	int X_new_rows = (int) paramv[2];
	int X_new_columns = (int) paramv[3];
	int k = (int) paramv[4];

	ocrGuid_t A = depv[0].guid;
	ocrGuid_t B = depv[1].guid;
	ocrGuid_t x_new = depv[2].guid;
	ocrGuid_t p_old = depv[3].guid;
	ocrGuid_t r_new = depv[4].guid;
	ocrGuid_t rTr = depv[5].guid;
	ocrGuid_t e = depv[6].guid;

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

	double val = 0;
	double *p = (double*) depv[4].ptr;
	for (int i = 0; i < X_new_columns; i++)
		val += abs(p[i]);
	if (val <= RESIDUAL_LIMIT) {
		//ocrDbDestroy(r_old)
		//ocrDbDestroy(r_new)
		//if (k != 0) {
		//      ocrDbDestroy(x_old)
		//      ocrDbDestroy(p_old)
		//}
		//ocrDbDestroy(rTr)
		ocrEventSatisfy(e, x_new);
		return e;
	}

	ocrGuid_t rT_new;
	ocrEventCreate(&rT_new, OCR_EVENT_IDEM_T, true);
	createEdt(transposeEdtTemplate, rT_new, r_new, X_new_rows, X_new_columns);

	ocrGuid_t rT_newr;
	ocrEventCreate(&rT_newr, OCR_EVENT_IDEM_T, true);
	createEdt(productEdtTemplate, rT_newr, rT_new, X_new_columns,
		X_new_rows, r_new, X_new_rows, X_new_columns);

	ocrGuid_t beta;
	ocrEventCreate(&beta, OCR_EVENT_IDEM_T, true);
	createEdt(divideEdtTemplate, beta, rT_newr, 1, 1, rTr);

	//ocrDbDestroy(rT_new);
	//ocrDbDestroy(rT_newr);
	//ocrDbDestroy(rTr);
	ocrGuid_t bp;
	ocrEventCreate(&bp, OCR_EVENT_IDEM_T, true);
	createEdt(scaleEdtTemplate, bp, p_old, X_new_rows, X_new_columns, beta);

	ocrGuid_t p_new;
	ocrEventCreate(&p_new, OCR_EVENT_IDEM_T, true);
	createEdt(addEdtTemplate, p_new, r_new, X_new_rows, X_new_columns, bp);

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

	ocrGuid_t stepAEdtTemplate;
	ocrEdtTemplateCreate(&stepAEdtTemplate, CgEdt_StepA, 5, 6);

	u64 nparamv[5];
	nparamv[0] = (u64) A_rows;
	nparamv[1] = (u64) A_columns;
	nparamv[2] = (u64) X_new_rows;
	nparamv[3] = (u64) X_new_columns;
	nparamv[4] = (u64) k;

	ocrGuid_t stepAEdt;
	ocrEdtCreate(&stepAEdt, stepAEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &e);
	ocrAddDependence(A, stepAEdt, 0, DB_MODE_RO);
	ocrAddDependence(B, stepAEdt, 1, DB_MODE_RO);
	ocrAddDependence(x_new, stepAEdt, 2, DB_MODE_RO);
	ocrAddDependence(p_new, stepAEdt, 3, DB_MODE_RO);
	ocrAddDependence(r_new, stepAEdt, 4, DB_MODE_RO);
	ocrAddDependence(e, stepAEdt, 5, DB_MODE_RO);

	return e;
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
	while (1) {
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

		double val = 0;
		for (int i = 0; i < r_new->getColumns(); i++)
			val += abs(r_new->getValue(0, i));
		if (val <= RESIDUAL_LIMIT) {
			delete r_old;
			delete r_new;
			if (k != 0) {
				delete x_old;
				delete p_old;
			}
			delete rTr;
			break;
		}

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

ocrGuid_t conjugateGradient_OCR(Matrix *A, Matrix *x, Matrix *B)
{
	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 4, 2);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 2, 2);

	ocrGuid_t Ax;
	ocrEventCreate(&Ax, OCR_EVENT_IDEM_T, true);
	createEdt(productEdtTemplate, Ax, A->getDataBlock(), A->getRows(),
		A->getColumns(), x->getDataBlock(), x->getRows(), x->getColumns());

	ocrGuid_t r_old;
	ocrEventCreate(&r_old, OCR_EVENT_IDEM_T, true);
	createEdt(subtractEdtTemplate, r_old, B->getDataBlock(), B->getRows(),
		B->getColumns(), Ax, x->getRows(), x->getColumns());

	//ocrDbDestroy(Ax);

	int k = 0;

	//ocrGuid_t p_old = r_old;

	ocrGuid_t stepAEdtTemplate;
	ocrEdtTemplateCreate(&stepAEdtTemplate, CgEdt_StepA, 5, 5);

	u64 nparamv[5];
	nparamv[0] = (u64) A->getRows();
	nparamv[1] = (u64) A->getColumns();
	nparamv[2] = (u64) x->getRows();
	nparamv[3] = (u64) x->getColumns();
	nparamv[4] = (u64) k;

	ocrGuid_t e;
	ocrEventCreate(&e, OCR_EVENT_IDEM_T, true);
	ocrGuid_t stepAEdt;
	ocrEdtCreate(&stepAEdt, stepAEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &e);
	ocrAddDependence(A->getDataBlock(), stepAEdt, 0, DB_MODE_RO);
	ocrAddDependence(B->getDataBlock(), stepAEdt, 1, DB_MODE_RO);
	ocrAddDependence(x->getDataBlock(), stepAEdt, 2, DB_MODE_RO);
	ocrAddDependence(r_old, stepAEdt, 3, DB_MODE_RO); //p_old
	ocrAddDependence(r_old, stepAEdt, 4, DB_MODE_RO);
	ocrAddDependence(e, stepAEdt, 5, DB_MODE_RO);

	return e;
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

	Matrix *fishy = conjugateGradient(&A, &x, &B);
	cout << "Solution:" << endl;
	cout << fishy->getValue(0, 0) << endl;
	cout << fishy->getValue(1, 0) << endl << endl;

	ocrGuid_t result = conjugateGradient_OCR(&A, &x, &B);
	
	//return x_new;
//	u64 nparamv[1];
//	nparamv[0] = (u64)x.getRows();
//
//	ocrGuid_t printEdtTemplate;
//	ocrEdtTemplateCreate(&printEdtTemplate, printEdt, 1, 2);
//
//	ocrGuid_t printEdt;
//	ocrGuid_t e;
//	ocrEdtCreate(&printEdt, printEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
//		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
//	ocrAddDependence(result, printEdt, 0, DB_MODE_RO);
//	ocrAddDependence(e, printEdt, 1, DB_MODE_RO);

	ocrShutdown();
	cout << "done" << endl;
	return NULL_GUID;
}

