#include "OCRDenseMatrix.h"

OCRDenseMatrix::OCRDenseMatrix(void *_data, ocrGuid_t _datablock) {
    ASSERT(_data != NULL);
    data = _data;
    datablock = _datablock;
    rows = *((unsigned int*)_data);
    columns = *(((unsigned int*)_data)+1);
    values = (double*)(_data+sizeof(unsigned int)*2);
}

OCRDenseMatrix::OCRDenseMatrix(unsigned int rs, unsigned int cols) {
    rows = rs;
    columns = cols;
    DBCREATE(&datablock, &data, sizeof(unsigned int)*2 + sizeof(double)*cols*rs,
             0, NULL_GUID, NO_ALLOC);
    *((unsigned int*)data) = rows;
    *(((unsigned int*)data)+1) = columns;
    values = (double*)(data+sizeof(unsigned int)*2);
}

double OCRDenseMatrix::getElement(unsigned int row, unsigned int column) const {
    return values[row*columns + column];
}

void OCRDenseMatrix::setElement(unsigned int row, unsigned int column, double val) {
    values[row*columns + column] = val;
}
