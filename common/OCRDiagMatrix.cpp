#include "OCRDiagMatrix.h"

OCRDiagMatrix::OCRDiagMatrix(void *_data, ocrGuid_t _datablock) {
    ASSERT(_data != NULL);
    data = _data;
    datablock = _datablock;
    rows = *((unsigned int*)_data);
    columns = rows;
    values = (double*)(_data+sizeof(unsigned int));
}

OCRDiagMatrix::OCRDiagMatrix(unsigned int sz) {
    rows = columns = sz;
    DBCREATE(&datablock, &data, sizeof(unsigned int) + sizeof(double)*sz,
             0, NULL_GUID, NO_ALLOC);
    *((unsigned int*)data) = sz;
    values = (double*)(data+sizeof(unsigned int));
}

double OCRDiagMatrix::getDeterminant() const {
    double prod = 1.0;
    for(unsigned int i=0;i<rows;i++) {
        prod *= values[i];
    }
    return prod;
}

OCRMatrix *OCRDiagMatrix::getInverse() const {
    OCRDiagMatrix *ret = new OCRDiagMatrix(columns);
    // The inverse of a diagonal matrix is just 1 / x for each entry x 
    for(unsigned int i=0;i<rows;i++) {
        ret->setElement(i,i,1.0/values[row]);
    }
    return ret;
}

double OCRDiagMatrix::getElement(unsigned int row, unsigned int column) const {
    if(row != column) {
        return 0;
    } else {
        return values[row];
    }
}

void OCRDiagMatrix::setElement(unsigned int row, unsigned int column, double val) {
    assert(row == column);
    values[row] = val;
}
