/*
 * File:   matrix.cpp
 * Author: commlaptop
 *
 * Created on February 7, 2014, 1:39 PM
 */

#include "matrix.h"

using namespace std;

Matrix::Matrix() {
	Matrix(3, 3);
}

Matrix::Matrix(int row, int column) {
	this->columns = column;
	this->rows = row;
	if (columns <= 0) this->columns = 1;
	if (rows <= 0) this->rows = 1;
	DBCREATE(&dataBlock, (void**)&mat, this->columns * this->rows * sizeof(double), 
             DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	for (int i = 0; i < (this->columns * this->rows); i++)
		mat[i] = 0.0f;
}

Matrix::Matrix(Matrix *mat) {
	columns = mat->getColumns();
	rows = mat->getRows();
	DBCREATE(&dataBlock, (void**)&this->mat, this->columns * this->rows * sizeof(double), 
             DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	for (int j = 0; j < rows; j++)
		for (int i = 0; i < columns; i++)
			setValue(j, i, mat->getValue(j, i));
}

Matrix::~Matrix()
{
	//delete mat? wont work for some reason
}

 int Matrix::getRows()
{
	return rows;
}

 int Matrix::getColumns()
{
	return columns;
}

double Matrix::getValue(int row, int column)
{
	if (column >= columns) return 0.0f;
	if (row >= rows) return 0.0f;
	return mat[row*columns + column];
}

void Matrix::setValue(int row, int column, double value)
{
	if ((column < columns) && (row < rows))
		mat[row*columns + column] = value;
}

ocrGuid_t Matrix::getDataBlock() {
	return dataBlock;
}



Matrix* matrixScale(Matrix *A, double scalar)
{
	Matrix* r = new Matrix(A->getRows(), A->getColumns());
	for (int row = 0; row < r->getRows(); row++) {
		for (int column = 0; column < r->getColumns(); column++) {
			double v = A->getValue(row, column);
			r->setValue(row, column, v * scalar);
		}
	}
	return r;
}

double* matrixScale(ocrGuid_t *dataBlock, double *A_mat, int A_rows, int A_columns, double scalar)
{
	double* r;
	DBCREATE(dataBlock, (void**)&r, A_columns * A_rows * sizeof(double), 
             DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	for (int row = 0; row < A_rows; row++) {
		for (int column = 0; column < A_columns; column++) {
			double v = A_mat[row*A_columns+column];
			r[row*A_columns+column] = v * scalar;
		}
	}
	return r;
}

Matrix* matrixProduct(Matrix *A, Matrix *B)
{
	if (A->getColumns() != B->getRows()) return NULL;

	Matrix* r = new Matrix(A->getRows(), B->getColumns());

	for (int i = 0; i < A->getRows(); i++) {
		for (int j = 0; j < A->getColumns(); j++) {
			double val = 0;
			for (int x = 0; x < A->getColumns(); x++) {
				val += A->getValue(i, x) * B->getValue(x, j);
			}
			r->setValue(i, j, val);
		}
	}
	return r;
}

double* matrixProduct(ocrGuid_t *dataBlock, double *A_mat, int A_rows,
                      int A_columns, double *B_mat, int B_rows, int B_columns)
{
	if (A_columns != B_rows) return NULL;

	double* r;
	DBCREATE(dataBlock, (void**)&r, A_rows * B_columns * sizeof(double), 
             DB_PROP_NONE, NULL_GUID, NO_ALLOC);

	for (int i = 0; i < A_rows; i++) {
		for (int j = 0; j < A_columns; j++) {
			r[i*B_columns+j] = 0;
			for (int x = 0; x < A_columns; x++) {
				r[i*B_columns+j] += A_mat[i*A_columns+x] * B_mat[x*B_columns+j];
			}
		}
	}
	return r;
}

Matrix* matrixTranspose(Matrix *A)
{
	Matrix* r = new Matrix(A->getColumns(), A->getRows());
	for (int row = 0; row < A->getRows(); row++) {
		for (int column = 0; column < A->getColumns(); column++) {
			double val = A->getValue(row, column);
			r->setValue(column, row, val);
		}
	}
	return r;
}

double* matrixTranspose(ocrGuid_t *dataBlock, double *A_mat, int A_rows, int A_columns)
{
	double* r;
	DBCREATE(dataBlock, (void**)&r, A_columns * A_rows * sizeof(double), 
             DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	for (int row = 0; row < A_rows; row++) {
		for (int column = 0; column < A_columns; column++) {
			double v = A_mat[row*A_columns+column];
			r[column*A_columns+row] = v;
		}
	}
	return r;
}

Matrix* matrixAdd(Matrix *A, Matrix *B)
{
	if (A->getRows() != B->getRows()) return NULL;
	if (A->getColumns() != B->getColumns()) return NULL;
	Matrix* r = new Matrix(A->getRows(), A->getColumns());
	for (int row = 0; row < A->getRows(); row++) {
		for (int column = 0; column < A->getColumns(); column++) {
			double valA = A->getValue(row, column);
			double valB = B->getValue(row, column);
			r->setValue(row, column, valA + valB);
		}
	}
	return r;
}

double* matrixAdd(ocrGuid_t *dataBlock, double *A_mat, int A_rows,
                  int A_columns, double *B_mat, int B_rows, int B_columns)
{
	if (A_rows != B_rows) return NULL;
	if (A_columns != B_columns) return NULL;
	double* r;
	DBCREATE(dataBlock, (void**)&r, A_columns * A_rows * sizeof(double), 
             DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	for (int row = 0; row < A_rows; row++) {
		for (int column = 0; column < A_columns; column++) {
			double valA = A_mat[row*A_columns+column];
			double valB = B_mat[row*A_columns+column];
			r[row*A_columns+column] = valA + valB;
		}
	}
	return r;
}

Matrix* matrixSubtract(Matrix *A, Matrix *B)
{
	if (A->getRows() != B->getRows()) return NULL;
	if (A->getColumns() != B->getColumns()) return NULL;
	Matrix* r = new Matrix(A->getRows(), A->getColumns());
	for (int row = 0; row < A->getRows(); row++) {
		for (int column = 0; column < A->getColumns(); column++) {
			double valA = A->getValue(row, column);
			double valB = B->getValue(row, column);
			r->setValue(row, column, valA - valB);
		}
	}
	return r;
}

double* matrixSubtract(ocrGuid_t *dataBlock, double *A_mat, int A_rows,
                       int A_columns, double *B_mat, int B_rows, int B_columns)
{
	if (A_rows != B_rows) return NULL;
	if (A_columns != B_columns) return NULL;
	double* r;
	DBCREATE(dataBlock, (void**)&r, A_columns * A_rows * sizeof(double), 
             DB_PROP_NONE, NULL_GUID, NO_ALLOC);
	for (int row = 0; row < A_rows; row++) {
		for (int column = 0; column < A_columns; column++) {
			double valA = A_mat[row*A_columns+column];
			double valB = B_mat[row*A_columns+column];
			r[row*A_columns+column] = valA - valB;
		}
	}
	return r;
}