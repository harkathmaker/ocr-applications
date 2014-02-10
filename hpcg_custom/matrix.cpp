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
	this->mat = new double[this->columns * this->rows];
	for (int i = 0; i < (this->columns * this->rows); i++)
		mat[i] = 0.0f;
}

Matrix::Matrix(Matrix *mat) {
	columns = mat->getColumns();
	rows = mat->getRows();
	this->mat = new double[columns * rows];
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

void Matrix::matrixScale(double value)
{
	for (int row = 0; row < getRows(); row++) {
		for (int column = 0; column < getColumns(); column++) {
			double v = getValue(row, column);
			setValue(row, column, v * value);
		}
	}
}

void Matrix::matrixAdd(Matrix *B)
{
	if (getRows() != B->getRows()) return;
	if (getColumns() != B->getColumns()) return;
	for (int row = 0; row < getRows(); row++) {
		for (int column = 0; column < getColumns(); column++) {
			double valA = getValue(row, column);
			double valB = B->getValue(row, column);
			setValue(row, column, valA + valB);
		}
	}
}

void Matrix::matrixSubtract(Matrix *B)
{
	if (getRows() != B->getRows()) return;
	if (getColumns() != B->getColumns()) return;
	for (int row = 0; row < getRows(); row++) {
		for (int column = 0; column < getColumns(); column++) {
			double valA = getValue(row, column);
			double valB = B->getValue(row, column);
			setValue(row, column, valA - valB);
		}
	}
}

Matrix* matrixScale(Matrix *A, double value)
{
	Matrix* r = new Matrix(A->getRows(), A->getColumns());
	for (int row = 0; row < r->getRows(); row++) {
		for (int column = 0; column < r->getColumns(); column++) {
			double v = A->getValue(row, column);
			r->setValue(row, column, v * value);
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
