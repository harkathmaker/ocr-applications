/*
 * OCRDiagMatrix.h
 * Diagonal representation of a matrix, constructed
 * with an OCR datablock.
 */

#ifndef __OCR_APPLICATIONS_DIAG_MATRIX
#define __OCR_APPLICATIONS_DIAG_MATRIX

#include "OCRBaseMatrix.h"

class OCRDiagMatrix: public OCRMatrix {
public:
    // Constructs the dense matrix from a datablock.
    OCRDiagMatrix(void *_data, ocrGuid_t _datablock=NULL_GUID);
    // Constructs a dense matrix from a sparse matrix.
    //OCRDiagMatrix(OCRSparseMatrix *spm);
    // Constructs a dense matrix with the given number of rows and columns.
    // A new datablock is allocated to hold this data.
    OCRDiagMatrix(unsigned int size);

    double getDeterminant() const;

    OCRMatrix *getInverse() const;

    double getElement(unsigned int row, unsigned int column) const;
    void setElement(unsigned int row, unsigned int column, double val);
private:
    double *values;
};

#endif
