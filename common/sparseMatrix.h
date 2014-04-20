/*
 * OCRSparseMatrix.h
 * Sparse representation of a matrix, constructed with an OCR
 * datablock.
 */

#ifndef __OCR_APPLICATIONS_SPARSE_MATRIX
#define __OCR_APPLICATIONS_SPARSE_MATRIX

#include "OCRBaseMatrix.h"

class OCRSparseMatrix: public OCRMatrix {
public:
    OCRSparseMatrix(void *_data, double defaultVal, ocrGuid_t _datablock=NULL_GUID);

    // Set the value for non-specified entries (normally this is set to 0).
    void setDefaultValue(double val);
private:
    // Indicates whether entries has been populated or not.
    // It is lazily allocated when required by a multiply operation.
    bool cached;
    std::map<unsigned int, double> entries;

    // Caches entries from array.
    void mapEntries();

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
