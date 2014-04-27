/*
 * OCRSparseMatrix.h
 * Sparse representation of a matrix, constructed with an OCR
 * datablock.
 */

#ifndef __OCR_APPLICATIONS_SPARSE_MATRIX
#define __OCR_APPLICATIONS_SPARSE_MATRIX

#include "OCRBaseMatrix.h"
#include "OCRDenseMatrix.h"

#include <map>

class OCRSparseMatrix: public OCRMatrix {
public:
    friend class OCRMatrix;

    // Constructs a sparse matrix from an already existing datablock.
    OCRSparseMatrix(void *_data, double defaultVal, ocrGuid_t _datablock=NULL_GUID);

    // Constructs a new sparse matrix, allocating a datablock in the process.
    // The specified capacity should be large enough to contain all future modifications
    // to the matrix.
    OCRSparseMatrix(unsigned int rows, unsigned int columns, unsigned int capacity);

    // Converts a dense matrix to a sparse matrix.
    OCRSparseMatrix(OCRDenseMatrix *m, double defaultVal);

    // Set the value for non-specified entries (normally this is set to 0).
    void setDefaultValue(double val);

    // Overridden from OCRMatrix
    void setElement(unsigned int row, unsigned int col, double val);
    double getElement(unsigned int row, unsigned int col) const;
private:
    // Indicates whether entries has been populated or not.
    // It is lazily allocated when required by a multiply operation.
    mutable bool cached;
    mutable std::map<unsigned int, double> entries;

    // Caches entries from array.
    void cacheEntries() const;

    double defaultValue;

    unsigned int capacity;
    // Reference to the size of the matrix. If entries are added this is updated
    // in the datablock.
    unsigned int *size;

    // Array of "keys" (row * columns + col) for all non-default entries
    unsigned int *indices;
    // Array of the non-default values
    double *values;
};

#endif
