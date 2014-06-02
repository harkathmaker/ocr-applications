/*
 * File:   main.cpp
 * Author: Grady W. Ellison
 *
 * Created on February 7, 2014, 1:38 PM
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <cstdlib>
#include <vector>

extern "C" {
#include <pthread.h>
}

//Define DEBUG_MESSAGES to view debugging messages on what the program is doing.
//The messages can (and will) print out of order, but it will give an idea on
//what is going on in the program
//#define DEBUG_MESSAGES

#ifndef __OCR__
#define __OCR__
#endif

#include "ocr.h"

#include "matrix.h"
#include "timer.h"


//Matrix size n by n
unsigned int MATRIX_N = 500;
//Amount of iterations of the algorithm on a given problem
unsigned int K_ITERATIONS = MATRIX_N * 2;
//Whether to treat A as a sparse matrix
bool isSparse = false;
//Amount of elements in sparse matrix A
int elementAmount = 0;
//Whether the residual is low enough to trust the current iterated solution
bool isResidualLow = false;
//If all the residuals in the residual matrix are below this amount, the current iterated solution is returned
double RESIDUAL_LIMIT = 1e-10;
int finalIterationAmount = 0;

using namespace std;

//Generates a pseudorandom double from fMin to fMax
//Parameters:
//  fMin: Lowest value that the random double can be
//  fMin: Highest value that the random double can be
//Output:
//  Random double from fMin to fMax
double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

//Multiplies datablock A by a constant value and returns a new datablock R as the result
//Parameters:
//  A_rows: Number of rows in matrix A
//  A_columns: Number of columns in matrix A
//Dependencies:
//  A_db: Matrix A
//  scalar: Scalar double constant to scale A by
//Output:
//  Matrix A that is scaled by the cosntant scalar.
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

//Multiplies datablock A by datablock B and returns a new datablock R as the result
//Parameters:
//  A_rows: Number of rows in matrix A
//  A_columns: Number of columns in matrix A
//  B_rows: Number of rows in matrix B
//  B_columns: Number of columns in matrix B
//Dependencies: A_db, B_db
//  A_db: Matrix A
//  B_db: Matrix B
//Output:
//  The resulting matrix of dense matrix A multiplied by dense matrix B.
extern "C" ocrGuid_t productEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixProduct((double*) depv[0].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[1].ptr, (int) paramv[2], (int) paramv[3]);
#ifdef DEBUG_MESSAGES
	cout << "productEdt()" << endl;
#endif
	return dataBlock;
}

//Multiplies sparse datablock matrix A by dense datablock matrix B and returns a new datablock R as the result
//Parameters:
//  A_rows: Number of rows in matrix A
//  A_columns: Number of columns in matrix A
//  B_rows: Number of rows in matrix B
//  B_columns: Number of columns in matrix B
//Dependencies:
//  A_db: Scalar matrix A
//  A_elementList: Element list for sparse matrix A. Contains positions on the values for matrix A.
//  B_db: Dense Matrix B
//Output:
// The resulting dense matrix of sparse matrix A multiplied by dense matrix B
extern "C" ocrGuid_t productEdt_sparse(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixProduct_sparse((double*) depv[0].ptr, (unsigned int*) depv[1].ptr,
                                               elementAmount, (int)paramv[0], (int)paramv[1],
                                               (double*) depv[2].ptr, (int)paramv[2], (int)paramv[3]);
#ifdef DEBUG_MESSAGES
	cout << "productEdt_sparse()" << endl;
#endif
	return dataBlock;
}

//Transposes datablock A and returns a new datablock R as the result
//Parameters:
//  A_rows: Number of rows in matrix A
//  A_columns: Number of columns in matrix A
//Dependencies:
//  A_db: Scalar matrix A
//Output:
//  The resulting matrix that is a transpose of matrix A
extern "C" ocrGuid_t transposeEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixTranspose((double*) depv[0].ptr, (int) paramv[0], (int) paramv[1]);
#ifdef DEBUG_MESSAGES
	cout << "transposeEdt()" << endl;
#endif
	return dataBlock;
}

//Adds datablock A by datablock B and returns a new datablock R as the result
//Parameters:
//  A_rows: Number of rows in matrix A
//  A_columns: Number of columns in matrix A
//Dependencies: (A and B need same rows/columns)
//  A_db: Matrix A
//  B_db: Matrix B
//Output:
//  The resulting matrix from matrix A + matrix B.
extern "C" ocrGuid_t addEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixAdd((double*) depv[0].ptr, (int) paramv[0], (int) paramv[1], (double*) depv[1].ptr);
#ifdef DEBUG_MESSAGES
	cout << "addEdt()" << endl;
#endif
	return dataBlock;
}

//Subtracts datablock A by datablock B and returns a new datablock R as the result
//Parameters:
//  A_rows: Number of rows in matrix A
//  A_columns: Number of columns in matrix A
//Dependencies: (A and B need same rows/columns)
//  A_db: Matrix A
//  B_db: Matrix B
//Output:
//  The resulting matrix from matrix A - matrix B.
extern "C" ocrGuid_t subtractEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	ocrGuid_t dataBlock = matrixSubtract((double*) depv[0].ptr, (int) paramv[0], (int) paramv[1],
		(double*) depv[1].ptr);
#ifdef DEBUG_MESSAGES
	cout << "subtractEdt()" << endl;
#endif
	return dataBlock;
}

//Divides datablock A by datablock B and returns a new datablock R as the result
//Parameters:
//  None
//Dependencies:
//  A_db: Constant double A
//  B_db: Constant double B
//Output:
//  The resulting double from double A divided by double B
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

//Checks if the residuals are low enough to trust the current iterated solution
//and sets isResidualLow to true if it is. Satisfies residualEvent on completion.
//Parameters:
//  R_rows: The amount of rows in the residual matrix. The columns of the residual matrix is always 1.
//  residualEvent: Edt event that is satisfied on completion of this Edt
//Dependencies: R_db (double)
//  R_db: The residual matrix
//Output:
//  Sets isResidualLow to true if residuals are low enough. It is set to false otherwise.
extern "C" ocrGuid_t residualEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	double* r = (double*) depv[0].ptr;
	int R_rows = (int) paramv[0];
	bool tempFlag = true;
#ifdef DEBUG_MESSAGES
	cout << "residualEdt()" << endl;
#endif
	for (int i = 0; i < R_rows; i++) {
		if (r[i] > RESIDUAL_LIMIT) {
			ocrEventSatisfy((ocrGuid_t)paramv[1], NULL_GUID);
			return NULL_GUID;
		}
	}
	isResidualLow = true;
	ocrEventSatisfy((ocrGuid_t)paramv[1], NULL_GUID);
	return NULL_GUID;
}

//Prints the contents of a matrix datablock and stops a given timer for benchmarking purposes.
//Parameters:
//  x_rows: The amount of rows in the result matrix
//  t1.sec: A time_t value in seconds
//  t1.usec: A time_t value in milliseconds
//Dependencies:
//  A_db: The result matrix
//Output:
//  Prints the values of the result matrix as well as some benchmarking data
extern "C" ocrGuid_t printEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
#ifdef DEBUG_MESSAGES
	cout << "printEdt()" << endl;
#endif
	int rows = (int) paramv[0];
	time_t sec = (time_t) paramv[1];
	time_t usec = (suseconds_t) paramv[2];
	double* r = (double*) depv[0].ptr;

	timeval t2;
	t2.tv_sec = sec;
	t2.tv_usec = usec;
	double cgOcrTimeElapsed = tock(t2);

	cout << "Solution:" << endl;
	for (int i = 0; i < rows; i++)
		cout << r[i] << " ";
	cout << endl;

	cout << "Time to test Conjugate Gradient OCR: " << cgOcrTimeElapsed << " ms" << endl;
	int flops;
	if (isSparse == false)
		flops = ((2*rows*rows + rows*rows) + (finalIterationAmount * (2*rows*rows + 6*rows*rows + 3*rows + 2)));
	else
		flops = ((elementAmount*2 + rows*rows) + (finalIterationAmount * (elementAmount*2 + 6*rows*rows + 3*rows + 2)));
	cout << "Iterations: " << finalIterationAmount << endl; //<< "Iterations per second: " << ((float)finalIterationAmount / (cgOcrTimeElapsed / 1000.)) << endl;
	//if (isSparse == true)
	//	cout << "elementAmount: " << elementAmount << endl;
	cout << "MFLOPS: " << (((float)flops / 1000000.) / (cgOcrTimeElapsed / 1000.)) << endl;

#ifdef DEBUG_MESSAGES
	cout << "OCR Shutdown..." << endl;
#endif
	ocrShutdown();
	return NULL_GUID;
}

//Performs the conjugate gradient algorithm with matrix A, x, and B. Will
//recursively call itself for a number of iterations until the residual matrix
//contains values below a set threshold, which is when the solution matrix is returned
//Parameters:
//  A_rows: Number of rows in matrix A
//  A_columns: Number of columns in matrix A
//  X_old_rows: Number of rows in matrix X (and B)
//  X_old_columns: Number of columns in matrix X (and B)
//  k: Number of algorithm iterations currently processed
//  doneEvent: Event that is satisfied when the conjugate gradient algorithm is finished
//Dependencies: All datablock dependencies will be destroyed in the end
//  A_db: Matrix A
//  B_db: Matrix B
//  X_old_db: Matrix X  (X_new_db from previous iteration)
//  P_old_db: Search direction matrix (P_new_db from previous iteration)
//  R_old_db: Residual matrix (R_new_db from previous iteration)
//  elementList: Element list for sparse matrix A (if A_db is sparse)
//  For iterations k > 0, pass in the following datablocks (generated by
//    this Edt) to cleanup at the start of each iteration (15 in total):
//      rT, rTr, pT, Ap, pTAp, alpha, alpha_p, aAp, rT_new, rT_newr, beta, bp
//Output:
//  The solution matrix
extern "C" ocrGuid_t CgEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	//Edt parameters
	ocrGuid_t result = (ocrGuid_t) paramv[0];
	int A_rows = (int) paramv[1];
	int A_columns = (int) paramv[2];
	int X_old_rows = (int) paramv[3];
	int X_old_columns = (int) paramv[4];
	int k = (int) paramv[5];
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Retrieved parameters: CgEdt" << endl;
#endif

	//Edt dependencies
	ocrGuid_t A = depv[0].guid;
	ocrGuid_t B = depv[1].guid;
	ocrGuid_t x_old = depv[2].guid;
	//p0 = r0
	ocrGuid_t p_old = depv[3].guid;
	ocrGuid_t r_old = depv[4].guid;
	ocrGuid_t elementList;
	int m = 5;
	if (isSparse == true) {
		elementList = depv[5].guid;
		m++;
	}

	//Cleans up datablocks from previous iteration
	if (k > 0) {	
		for (int i = 0; i < 12; i++) {
			ocrDbDestroy(depv[m].guid);
			m++;
		}
	}
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Retrieved dependencies: CgEdt" << endl;
#endif

	//If the algorithm is executed K times on the conjugate gradient
	//problem or if the residuals are low enough, return the current result x
	if ((k == K_ITERATIONS) || (isResidualLow == true)){
		finalIterationAmount = k;
		ocrEventSatisfy(result, x_old);
#ifdef DEBUG_MESSAGES
		cout << "k = " << k << ". Satisfied result of CgEdt" << endl;
#endif
		return NULL_GUID;
	}

	u64 edtParams[4];

	//Edt templates
	ocrGuid_t scaleEdtTemplate;
	ocrEdtTemplateCreate(&scaleEdtTemplate, scaleEdt, 2, 2);

	ocrGuid_t scaleEdtTemplate_N;
	ocrEdtTemplateCreate(&scaleEdtTemplate_N, scaleEdt, 2, 3);

	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 4, 2);

	ocrGuid_t productEdtSparseTemplate;
	ocrEdtTemplateCreate(&productEdtSparseTemplate, productEdt_sparse, 4, 3);

	ocrGuid_t transposeEdtTemplate;
	ocrEdtTemplateCreate(&transposeEdtTemplate, transposeEdt, 2, 1);

	ocrGuid_t transposeEdtTemplate_K;
	ocrEdtTemplateCreate(&transposeEdtTemplate_K, transposeEdt, 2, 2);

	ocrGuid_t addEdtTemplate;
	ocrEdtTemplateCreate(&addEdtTemplate, addEdt, 2, 2);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 2, 2);

	ocrGuid_t divideEdtTemplate;
	ocrEdtTemplateCreate(&divideEdtTemplate, divideEdt, 0, 2);

	ocrGuid_t residualEdtTemplate;
	ocrEdtTemplateCreate(&residualEdtTemplate, residualEdt, 2, 1);

#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt templates: scale, product, transpose, add, subtract, divide" << endl;
#endif

	//EdtA through EdtB calculates rTr
	//EdtC through EdtE calculates pTAp
	//EdtF calculates alpha = rTr / pTAp

	//rT
	ocrGuid_t rT;
	ocrGuid_t edtA;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtA, transposeEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &rT);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtA" << endl;
#endif

	//rTr
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

	//pT
	ocrGuid_t pT;
	ocrGuid_t edtC;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtC, transposeEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &pT);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtC" << endl;
#endif

	//A * p_old
	ocrGuid_t Ap;
	ocrGuid_t edtD;
	edtParams[0] = A_rows;
	edtParams[1] = A_columns;
	edtParams[2] = X_old_rows;
	edtParams[3] = X_old_columns;
	if (isSparse == false)
		ocrEdtCreate(&edtD, productEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
			NULL, EDT_PROP_NONE, NULL_GUID, &Ap);
	else
		ocrEdtCreate(&edtD, productEdtSparseTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
			NULL, EDT_PROP_NONE, NULL_GUID, &Ap);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtD" << endl;
#endif

	//pTAp
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

	//alpha = rTr / pTAp
	ocrGuid_t alpha;
	ocrGuid_t edtF;
	ocrEdtCreate(&edtF, divideEdtTemplate, EDT_PARAM_DEF, NULL, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &alpha);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtF" << endl;
#endif

	//EdtG through EdtH calculates x_new = x_old + alpha * p_old

	//alpha * p_old
	ocrGuid_t alpha_p;
	ocrGuid_t edtG;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtG, scaleEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &alpha_p);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtG" << endl;
#endif

	//x_new = x_old + alpha * p_old
	ocrGuid_t x_new;
	ocrGuid_t edtH;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtH, addEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &x_new);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtH" << endl;
#endif

	//EdtI through EdtJ calculates r_new = r_old - alpha * A * p_old
	//Note how A * p_old is reused from before

	//alpha * A * p_old
	ocrGuid_t aAp;
	ocrGuid_t edtI;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtI, scaleEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &aAp);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtI" << endl;
#endif

	//r_new = r_old - alpha * A * p_old
	ocrGuid_t r_new;
	ocrGuid_t edtJ;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtJ, subtractEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &r_new);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtJ" << endl;
#endif

	//Residual check
	ocrGuid_t residualEvent;
	ocrEventCreate(&residualEvent, OCR_EVENT_ONCE_T, false);

	ocrGuid_t edtResidual;
	edtParams[0] = X_old_rows;
	edtParams[1] = (ocrGuid_t)residualEvent;
	ocrEdtCreate(&edtResidual, residualEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtResidual" << endl;
#endif

	//EdtK through EdtM calculates beta = (rT_new * r_new) / (rT_old * r_old)
	//Note how rT_old * r_old is reused from before

	//rT_new
	ocrGuid_t rT_new;
	ocrGuid_t edtK;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtK, transposeEdtTemplate_K, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &rT_new);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtK" << endl;
#endif

	//rT_new * r_new
	ocrGuid_t rT_newr;
	ocrGuid_t edtL;
	edtParams[0] = X_old_columns;
	edtParams[1] = X_old_rows;
	edtParams[2] = X_old_rows;
	edtParams[3] = X_old_columns;
	ocrEdtCreate(&edtL, productEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &rT_newr);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtL" << endl;
#endif

	//beta = (rT_new * r_new) / (rT_old * r_old)
	ocrGuid_t beta;
	ocrGuid_t edtM;
	ocrEdtCreate(&edtM, divideEdtTemplate, EDT_PARAM_DEF, NULL, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &beta);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtM" << endl;
#endif

	//EdtN through EdtO calculates p_new = r_new + beta * p_old

	//beta * p_old
	ocrGuid_t bp;
	ocrGuid_t edtN;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtN, scaleEdtTemplate_N, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &bp);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtN" << endl;
#endif

	//p_new = r_new + beta * p_old
	ocrGuid_t p_new;
	ocrGuid_t edtO;
	edtParams[0] = X_old_rows;
	edtParams[1] = X_old_columns;
	ocrEdtCreate(&edtO, addEdtTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &p_new);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: edtO" << endl;
#endif

	//CgEdt template
	ocrGuid_t CgEdtTemplate;
	if (isSparse == false)
		ocrEdtTemplateCreate(&CgEdtTemplate, CgEdt, 6, 5+12);
	else
		ocrEdtTemplateCreate(&CgEdtTemplate, CgEdt, 6, 6+12);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt template: CgEdt" << endl;
#endif

	//Sets up parameters for the next CgEdt iteration
	//k = k + 1
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

	//Create another CgEdt to spawn (until K iterations are completed)
	ocrGuid_t myEdt;
	//Satisfy with A, B, x_new, p_new, r_new
	ocrEdtCreate(&myEdt, CgEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Created Edt: CgEdt" << endl;
#endif

	//Add dependencies of the Edts in reverse order to ensure they don't complete early
	ocrAddDependence(A, myEdt, 0, DB_MODE_RO);
	ocrAddDependence(B, myEdt, 1, DB_MODE_RO);
	ocrAddDependence(x_new, myEdt, 2, DB_MODE_RO);
	ocrAddDependence(p_new, myEdt, 3, DB_MODE_RO);
	ocrAddDependence(r_new, myEdt, 4, DB_MODE_RO);
	m = 5;	
	if (isSparse == true) {
		ocrAddDependence(elementList, myEdt, 5, DB_MODE_RO);
		m++;
	}
	ocrAddDependence(rT, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(rTr, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(pT, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(Ap, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(pTAp, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(alpha, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(alpha_p, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(aAp, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(rT_new, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(rT_newr, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(beta, myEdt, m, DB_MODE_RO); m++;
	ocrAddDependence(bp, myEdt, m, DB_MODE_RO); m++;
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: CgEdt" << endl;
#endif

	ocrAddDependence(r_new, edtO, 0, DB_MODE_RO);
	ocrAddDependence(bp, edtO, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtO" << endl;
#endif

	ocrAddDependence(p_old, edtN, 0, DB_MODE_RO);
	ocrAddDependence(beta, edtN, 1, DB_MODE_RO);
	ocrAddDependence(residualEvent, edtN, 2, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtN" << endl;
#endif

	ocrAddDependence(rT_newr, edtM, 0, DB_MODE_RO);
	ocrAddDependence(rTr, edtM, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtM" << endl;
#endif

	ocrAddDependence(rT_new, edtL, 0, DB_MODE_RO);
	ocrAddDependence(r_new, edtL, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtL" << endl;
#endif

	ocrAddDependence(r_new, edtK, 0, DB_MODE_RO);
	ocrAddDependence(residualEvent, edtK, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtK" << endl;
#endif

	ocrAddDependence(r_new, edtResidual, 0, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtResidual" << endl;
#endif

	ocrAddDependence(r_old, edtJ, 0, DB_MODE_RO);
	ocrAddDependence(aAp, edtJ, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtJ" << endl;
#endif

	ocrAddDependence(Ap, edtI, 0, DB_MODE_RO);
	ocrAddDependence(alpha, edtI, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtI" << endl;
#endif

	ocrAddDependence(x_old, edtH, 0, DB_MODE_RO);
	ocrAddDependence(alpha_p, edtH, 1, DB_MODE_RO);
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

	ocrAddDependence(pT, edtE, 0, DB_MODE_RO);
	ocrAddDependence(Ap, edtE, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "k = " << k << ". Added dependencies: edtE" << endl;
#endif

	ocrAddDependence(A, edtD, 0, DB_MODE_RO);
	if (isSparse == false) {
		ocrAddDependence(p_old, edtD, 1, DB_MODE_RO);
	} else {
		ocrAddDependence(elementList, edtD, 1, DB_MODE_RO);
		ocrAddDependence(p_old, edtD, 2, DB_MODE_RO);
	}

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

//DEPRECATED: Sets up the conjugate gradient algorithm for matrices A, x, and B
//without OCR. Recently untested and may not work. Currently unused.
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

//Sets up the conjugate gradient algorithm Edt for matrices A, x, and B, where A is dense
//Parameters:
//  A_m: Matrix A
//  x_m: Matrix x
//  B_m: Matrix B
//  result: Where the result matrix is stored
//Output:
//  None
void conjugateGradient_OCR(Matrix *A_m, Matrix *x_m, Matrix *B_m, ocrGuid_t result)
{
	//Copies data from matrices A, x, B to their respective datablocks
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

	//Edt templates
	ocrGuid_t productEdtTemplate;
	ocrEdtTemplateCreate(&productEdtTemplate, productEdt, 4, 2);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 2, 2);

#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt templates: product, subtract" << endl;
#endif

	u64 edtParams[4];

	//Ax
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

	//residual (r0) = b - Ax0
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

	//Parameters for CgEdt
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

	//Creates the main algorithm Edt for the conjugate gradient problem
	ocrGuid_t myEdt;
	//Satisfy with A, B, x, p_old, r_old
	ocrEdtCreate(&myEdt, CgEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt: CgEdt" << endl;
#endif

	//Add dependencies for CgEdt
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

//Sets up the conjugate gradient algorithm Edt for matrices A, x, and B, where A is sparse
//Parameters:
//  A: Matrix A
//  elementList: Element list for sparse matrix A. Contains positions on the values for matrix A.
//  x_m: Matrix x
//  B_m: Matrix B
//  result: Where the result matrix is stored
//Output:
//  None
void conjugateGradient_OCR_sparse(ocrGuid_t A, ocrGuid_t elementList, Matrix *x_m, Matrix *B_m, ocrGuid_t result)
{
	//Copies data from matrices x, B to their respective datablocks
	ocrGuid_t B;
	ocrGuid_t x;
	ocrEventCreate(&B, OCR_EVENT_STICKY_T, true);
	ocrEventCreate(&x, OCR_EVENT_STICKY_T, true);
	ocrEventSatisfy(B, B_m->getDataBlock());
	ocrEventSatisfy(x, x_m->getDataBlock());
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Guids: B, x" << endl;
#endif

	//Edt templates
	ocrGuid_t productEdtSparseTemplate;
	ocrEdtTemplateCreate(&productEdtSparseTemplate, productEdt_sparse, 4, 3);

	ocrGuid_t subtractEdtTemplate;
	ocrEdtTemplateCreate(&subtractEdtTemplate, subtractEdt, 2, 2);

#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt templates: product, subtract" << endl;
#endif

	u64 edtParams[4];

	//Ax
	ocrGuid_t Ax;
	ocrGuid_t edtA;
	edtParams[0] = MATRIX_N; //A_rows
	edtParams[1] = MATRIX_N; //A_columns
	edtParams[2] = x_m->getRows();
	edtParams[3] = x_m->getColumns();
	ocrEdtCreate(&edtA, productEdtSparseTemplate, EDT_PARAM_DEF, edtParams, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, &Ax);
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt: Ax" << endl;
#endif

	//residual (r0) = b - Ax0
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
	ocrEdtTemplateCreate(&CgEdtTemplate, CgEdt, 6, 6);
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt template: CgEdt" << endl;
#endif

	//Parameters for CgEdt
	u64 nparamv[6];
	nparamv[0] = result;
	nparamv[1] = MATRIX_N; //A_rows
	nparamv[2] = MATRIX_N; //A_columns
	nparamv[3] = x_m->getRows();
	nparamv[4] = x_m->getColumns();
	nparamv[5] = k;

#ifdef DEBUG_MESSAGES
	cout << "Init. Initialized parameters: CgEdt" << endl;
#endif

	//Creates the main algorithm Edt for the conjugate gradient problem
	ocrGuid_t myEdt;
	//Satisfy with A, B, x, p_old, r_old, elementList
	ocrEdtCreate(&myEdt, CgEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_NONE, NULL_GUID, NULL);
#ifdef DEBUG_MESSAGES
	cout << "Init. Created Edt: CgEdt" << endl;
#endif

	//Add dependencies for CgEdt
	ocrAddDependence(A, myEdt, 0, DB_MODE_RO);
	ocrAddDependence(B, myEdt, 1, DB_MODE_RO);
	ocrAddDependence(x, myEdt, 2, DB_MODE_RO);
	ocrAddDependence(r_old, myEdt, 3, DB_MODE_RO); //p_old
	ocrAddDependence(r_old, myEdt, 4, DB_MODE_RO);
	ocrAddDependence(elementList, myEdt, 5, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: CgEdt" << endl;
#endif

	ocrAddDependence(B, edtB, 0, DB_MODE_RO);
	ocrAddDependence(Ax, edtB, 1, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: edtB" << endl;
#endif

	ocrAddDependence(A, edtA, 0, DB_MODE_RO);
	ocrAddDependence(elementList, edtA, 1, DB_MODE_RO);
	ocrAddDependence(x, edtA, 2, DB_MODE_RO);
#ifdef DEBUG_MESSAGES
	cout << "Init. Added dependencies: edtA" << endl;
#endif
}

//Reads a file containing data values and creates a matrix datablock from the data values
//Parameters:
//  fileName: The name of the matrix datafile
//  A: Where the loaded matrix is stored
//Output:
//  0 on success, -1 if file cannot be opened
int readMatrixFromFile(const char* fileName, Matrix &A)
{
	ifstream f;
	f.open(fileName);
	if (f.is_open() == false)
		return -1;

	double n;
	for (int row = 0; row < MATRIX_N; row++) {
		for (int column = 0; column < MATRIX_N; column++) {
			f >> n;
			A.setValue(row, column, n);
		}
	}
	f.close();
	return 0;
}

// Parses command line arguments, sets up initial matrices, calls the conjugate
// gradient algorithm helper functions, and sets up the printEdt
extern "C" ocrGuid_t mainEdt(u32 paramc, u64* paramv, u32 depc, ocrEdtDep_t depv[])
{
	//Command line argument parsing
	u64 argc = getArgc(depv[0].ptr);
	int i;
	if(argc < 3) {
		cout << "Argument 1: Whether A is to be treated as a sparse matrix. 0 is False, Non-zero is True." << endl;
		cout << "Argument 2: Matrix filename (filename.m)" << endl;
		cout << "Argument 3: Matrix size (in rows)" << endl;
		ocrShutdown();
		return NULL_GUID;
	}
#ifdef DEBUG_MESSAGES
	for(i=0;i<argc;i++) {
		char *argv = getArgv(depv[0].ptr,i);
		PRINTF("argv[%d]: %s\n",i,argv);
	}
#endif

	int foo = strtol(getArgv(depv[0].ptr,1),NULL,10);
	if (foo == 0)
		isSparse = false;
	else
		isSparse = true;
	MATRIX_N = strtol(getArgv(depv[0].ptr,3),NULL,10);

	//Reads data from a data file to matrix A
	Matrix A(MATRIX_N, MATRIX_N);
	if (readMatrixFromFile(getArgv(depv[0].ptr,2), A) == -1) {
		cout << "Error opening matrix file" << endl;
		ocrShutdown();
		return NULL_GUID;
	}

	//If matrix A is sparse, find non-zero elements in A and store their positions in an array
	ocrGuid_t A_sparse;
	ocrGuid_t A_sparseElementList;
	if (isSparse == true) {
		vector<unsigned int> tempPositions;
		vector<double> tempMatrix;
		for (int row = 0; row < MATRIX_N; row++) {
			for (int column = 0; column < MATRIX_N; column++) {
				double value = A.getValue(row, column);
				if (value != 0.00000000f) {
					tempPositions.push_back(row * MATRIX_N + column);
					tempMatrix.push_back(value);
				}
			}
		}
		elementAmount = tempMatrix.size();
		double *mat;
		unsigned int *positions;
		DBCREATE(&A_sparse, (void**) &mat, elementAmount * sizeof (double),
		         DB_PROP_NONE, NULL_GUID, NO_ALLOC);
		DBCREATE(&A_sparseElementList, (void**) &positions, elementAmount * sizeof (unsigned int),
		         DB_PROP_NONE, NULL_GUID, NO_ALLOC);
		for (int i = 0; i < elementAmount; i++) {
			mat[i] = tempMatrix[i];
			positions[i] = tempPositions[i];
		}
	}

	//Creates matrix x and B with a guess of all zeroes for x and random values for b
	Matrix x(MATRIX_N, 1);
	Matrix B(MATRIX_N, 1);
	for (int row = 0; row < MATRIX_N; row++) {
		x.setValue(row, 0, 1.0f);
		B.setValue(row, 0, fRand(-10.0f, 10.0f));
	}

	//Starts benchmarking timer
	timeval ocrCgTime;
	tick(ocrCgTime);
	cout << "Testing Conjugate Gradient OCR..." << endl;
	//Calls the conjugate gradient algorithm for OCR and stores the solutions in the result datablock
	ocrGuid_t result;
	ocrEventCreate(&result, OCR_EVENT_STICKY_T, true);
	if (isSparse == false)
		conjugateGradient_OCR(&A, &x, &B, result);
	else
		conjugateGradient_OCR_sparse(A_sparse, A_sparseElementList, &x, &B, result);

	//Parameters for printEdt
	u64 nparamv[3];
	nparamv[0] = x.getRows();
	nparamv[1] = ocrCgTime.tv_sec;
	nparamv[2] = ocrCgTime.tv_usec;

	//Creates print Edt (and template) and satisfies dependencies
	ocrGuid_t printEdtTemplate;
	ocrEdtTemplateCreate(&printEdtTemplate, printEdt, 3, 1);

	ocrGuid_t printEdt;
	ocrEdtCreate(&printEdt, printEdtTemplate, EDT_PARAM_DEF, nparamv, EDT_PARAM_DEF,
		NULL, EDT_PROP_FINISH, NULL_GUID, NULL);
	ocrAddDependence(result, printEdt, 0, DB_MODE_RO);

	return NULL_GUID;
}

