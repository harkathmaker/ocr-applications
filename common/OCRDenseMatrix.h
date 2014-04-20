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
    // Constructs the dense matrix from a datablock.
    OCRDenseMatrix(void *_data, ocrGuid_t _datablock=NULL_GUID);
    // Constructs a dense matrix from a sparse matrix.
    //OCRDenseMatrix(OCRSparseMatrix *spm);
    // Constructs a dense matrix with the given number of rows and columns.
    // A new datablock is allocated to hold this data.
    OCRDenseMatrix(unsigned int rows, unsigned int column);

    double getElement(unsigned int row, unsigned int column) const;
    void setElement(unsigned int row, unsigned int column, double val);
private:
    double *values;
};

#endif
