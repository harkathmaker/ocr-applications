/*
 * OCRBaseMatrix.h
 * Abstract base class for OCR matrix representation.
 * This class and its subclasses wrap a datablock containing a matrix,
 * allowing for easier data manipulation. All data is accessed through
 * the datablock, so no copying is done.
 */

#ifndef __OCR_APPLICATIONS_BASE_MATRIX
#define __OCR_APPLICATIONS_BASE_MATRIX

#define __OCR__

#include <ocr.h>
#include <ocr-errors.h>

class OCRSparseMatrix;
class OCRDenseMatrix;

class OCRMatrix {
public:
    OCRMatrix();
    virtual double getElement(unsigned int row, unsigned int column) const=0;
    virtual void setElement(unsigned int row, unsigned int column, double val)=0;

    void print() const;

    // Returns the determinant of the matrix. Only valid for square matrices.
    // This function can be overridden in child classes that can exploit special
    // patterns in their respective structures.
    // NOTE: computation of the determinant is almost always inefficient, and requires
    // additional memory to be allocated as well.
    virtual double getDeterminant() const;

    // Returns a pointer to a new matrix that is the inverse of the current one.
    // NOTE: inverse of a sparse matrix is usually dense. To solve systems such as
    // Ax = b, use LU decomposition or other techniques instead.
    virtual OCRMatrix *getInverse() const;

    // Performs matrix multiplication, outputting a sparse matrix.
    // If data is null, the result will be stored in a newly allocated datablock.
    // If data is non-null, the result will be stored in that address and no
    // datablock is created.
    // The returned OCRMatrix is created with new, and must be cleaned up with
    // delete.
    OCRSparseMatrix *multiplySparse(OCRMatrix *other, void *data=NULL, double defVal=0.0, unsigned int capacity=0) const;
    OCRDenseMatrix *multiplyDense(OCRMatrix *other, void *data=NULL) const;

    ocrGuid_t getDatablock() const;
    void *getData() const;

    unsigned int getRows() const;
    unsigned int getColumns() const;

protected:
    unsigned int rows;
    unsigned int columns;

    // If the matrix created a new datablock, this holds its guid value.
    // By default, its value is NULL_GUID.
    ocrGuid_t datablock;
    // Pointer to the matrix data.
    void *data;
};

#endif
