/*
 * File:   matrix.h
 * Author: commlaptop
 *
 * Created on February 7, 2014, 1:39 PM
 */

#ifndef MATRIX_H
#define	MATRIX_H

//#include <iostream>

extern "C" {
	#include <pthread.h>
}

#ifndef __OCR__
#define __OCR__
#endif

#include "ocr.h"

class Matrix {
public:
	Matrix();
	Matrix(int row, int column);
	Matrix(Matrix *mat);
	~Matrix();
	int getRows();
	int getColumns();
	double getValue(int row, int column);
	void setValue(int row, int column, double value);
	ocrGuid_t getDataBlock();
	//double* getMatrixData();

private:
	int rows;
	int columns;
	double* mat;
	ocrGuid_t dataBlock;
};

Matrix* matrixScale(Matrix *A, double scalar);
ocrGuid_t matrixScale(double *A_mat, int A_rows, int A_columns, double *scalar);
Matrix* matrixProduct(Matrix *A, Matrix *B);
ocrGuid_t matrixProduct(double *A_mat, int A_rows, int A_columns, double *B_mat, int B_rows, int B_columns);
Matrix* matrixTranspose(Matrix *A);
ocrGuid_t matrixTranspose(double *A_mat, int A_rows, int A_columns);
Matrix* matrixAdd(Matrix *A, Matrix *B);
ocrGuid_t matrixAdd(double *A_mat, int A_rows, int A_columns, double *B_mat);
Matrix* matrixSubtract(Matrix *A, Matrix *B);
ocrGuid_t matrixSubtract(double *A_mat, int A_rows, int A_columns, double *B_mat);

#endif	/* MATRIX_H */

