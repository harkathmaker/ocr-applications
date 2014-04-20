/*
 * OCRDenseMatrix.h
 * Dense (conventional) representation of a matrix, constructed
 * with an OCR datablock.
 */

#ifndef __OCR_APPLICATIONS_DENSE_MATRIX
#define __OCR_APPLICATIONS_DENSE_MATRIX

#include "OCRBaseMatrix.h"

class OCRDenseMatrix: public OCRMatrix {
public:
    OCRDenseMatrix(void *_data, ocrGuid_t _datablock=NULL_GUID);

private:
    double *values;
};

#endif
