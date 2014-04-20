#include "OCRDenseMatrix.h"

double OCRDenseMatrix::getElement(unsigned int row, unsigned int column) const {
    return values[row*columns + column];
}

void OCRDenseMatrix::setElement(unsigned int row, unsigned int column, double val) const {
    values[row*columns + column] = val;
}
