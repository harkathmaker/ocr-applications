/*
 * File:   matrix.h
 * Author: commlaptop
 *
 * Created on February 7, 2014, 1:39 PM
 */

#ifndef MATRIX_H
#define	MATRIX_H

#include <iostream>

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
	void matrixScale(double value);
	void matrixAdd(Matrix* B);
	void matrixSubtract(Matrix* B);

private:
	int rows;
	int columns;
	double *mat;
};

Matrix* matrixScale(Matrix* A, double value);
Matrix* matrixProduct(Matrix* A, Matrix* B);
Matrix* matrixTranspose(Matrix* B);
Matrix* matrixAdd(Matrix* A, Matrix* B);
Matrix* matrixSubtract(Matrix* A, Matrix* B);

#endif	/* MATRIX_H */

