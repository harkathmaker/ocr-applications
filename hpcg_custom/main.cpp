/*
 * File:   main.cpp
 * Author: commlaptop
 *
 * Created on February 7, 2014, 1:38 PM
 */

#include <iostream>
#include <cmath>

#include "matrix.h"

#define RESIDUAL_LIMIT 1.0e-9

using namespace std;

//do error checks later
Matrix* conjugateGradient(Matrix *A, Matrix *x, Matrix *B)
{
	Matrix *r_old;
	Matrix *r_new;
	Matrix *p_new;
	Matrix *x_new;
	Matrix bar = Matrix(x);
	Matrix *x_old = &bar;
	Matrix *Ax = matrixProduct(A, x_old);
	r_old = matrixSubtract(B, Ax);
	delete Ax;
	Matrix foo = Matrix(r_old);
	Matrix *p_old = &foo;

	int k = 0;
	while (1) {
		Matrix *rT = matrixTranspose(r_old);
		Matrix *rTr = matrixProduct(rT, r_old);
		Matrix *pT = matrixTranspose(p_old);
		Matrix *pTA = matrixProduct(pT, A);
		Matrix *pTAp = matrixProduct(pTA, p_old);
		double alpha = rTr->getValue(0,0) / pTAp->getValue(0,0);
		delete rT;
		delete pT;
		delete pTA;
		delete pTAp;

		Matrix *ap = matrixScale(p_old, alpha);
		x_new = matrixAdd(x_old, ap);
		delete ap;

		Matrix *aA = matrixScale(A, alpha);
		Matrix *aAp = matrixProduct(aA, p_old);
		r_new = matrixSubtract(r_old, aAp);
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
		double beta = rT_newr->getValue(0,0) / rTr->getValue(0,0);
		delete rT_new;
		delete rT_newr;
		delete rTr;

		Matrix *bp = matrixScale(p_old, beta);
		p_new = matrixAdd(r_new, bp);
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


int main(int argc, char** argv) {
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
	cout << fishy->getValue(0,0) << endl;
	cout << fishy->getValue(1,0) << endl;
	return 0;
}

