#include "OCRBaseMatrix.h"
#include "OCRSparseMatrix.h"
#include "OCRDenseMatrix.h"

OCRMatrix::OCRMatrix(): datablock(NULL_GUID), data(NULL) {
}
unsigned int OCRMatrix::getRows() const {
    return rows;
}
unsigned int OCRMatrix::getColumns() const {
    return columns;
}
ocrGuid_t OCRMatrix::getDatablock() const {
    return datablock;
}
void *OCRMatrix::getData() const {
    return data;
}

void OCRMatrix::print() const {
    for(unsigned int r=0;r<rows;r++) {
        PRINTF("[ ");
        for(unsigned int c=0;c<columns;c++) {
            PRINTF(" %.4e ",getElement(r,c));
        }
        PRINTF(" ]\n");
    }
}

OCRSparseMatrix *OCRMatrix::multiplySparse(OCRMatrix *other, void *dest, double defVal, unsigned int resCap) const {
    // Check that the matrix is a valid size
    ASSERT(columns == other->rows);
    unsigned int resultColumns = other->columns;
    unsigned int resultRows = rows;
    ocrGuid_t db = NULL_GUID;
    unsigned int matrixSize = resCap;
    // Keep computed entries while figuring out size of matrix
    std::map<unsigned int, double> rEntries;

    // TODO: not operating with adjacent data members. If this is a priority
    // we might need to represent sparse matrices in a different way
    unsigned int sparseCount = 0;
    for(unsigned int r=0;r<rows;r++) {
        for(unsigned int c=0;c<other->columns;c++) {
            double entry = 0;
            for(unsigned int i=0;i<columns;i++) {
                entry += getElement(r,i)*other->getElement(i,c);
            }
            if(entry != defVal) {
                rEntries.insert(std::pair<unsigned int, double>(r*resultColumns+c,entry));
                sparseCount++;
            }
        }
    }

    if(resCap < sparseCount) {
        matrixSize = sparseCount;
        if(resCap != 0) {
            // Warn the user their request has been overridden
            PRINTF("Requested capacity %lu too small to fit sparse matrix; set size to %lu\n",resCap,matrixSize);
        }
    }

    if(dest == NULL) {
        DBCREATE(&db, (void**)&dest,
                 sizeof(unsigned int)*(4+matrixSize)+sizeof(double)*matrixSize,
                 0,NULL_GUID,NO_ALLOC);
    }

    // Assign rows, columns, capacity and size
    *((unsigned int*)dest) = resultRows;
    *(((unsigned int*)dest)+1) = resultColumns;
    *(((unsigned int*)dest)+2) = matrixSize;
    *(((unsigned int*)dest)+3) = sparseCount;
    unsigned int *resultIndices = ((unsigned int*)dest)+4;
    double *resultVals = (double*)(dest+sizeof(unsigned int)*4);

    // Keep track of where we are in the sparse array
    OCRSparseMatrix *ret = new OCRSparseMatrix(dest,defVal,db);
    unsigned int k = 0;
    for(std::map<unsigned int,double>::iterator it=rEntries.begin();it!=rEntries.end();it++) {
        resultIndices[k] = it->first;
        resultVals[k] = it->second;
        k++;
    }
    // We've already calculated the entries, so copy them in
    ret->entries.swap(rEntries);
    return ret;
}

OCRDenseMatrix *OCRMatrix::multiplyDense(OCRMatrix *other, void *dest) const {
    // Check that the matrix is a valid size
    ASSERT(columns == other->rows);
    unsigned int resultColumns = other->columns;
    unsigned int resultRows = rows;
    ocrGuid_t db = NULL_GUID;
    if(dest == NULL) {
        unsigned int matrixSize = resultRows * resultColumns;
        
        DBCREATE(&db, (void**)&dest,
                 sizeof(unsigned int)*2+sizeof(double)*matrixSize,
                 0,NULL_GUID,NO_ALLOC);
    }

    // Assign rows and columns
    *((unsigned int*)dest) = resultRows;
    *(((unsigned int*)dest)+1) = resultColumns;
    double *resultVals = (double*)(dest+sizeof(unsigned int)*2);

    // TODO: not operating with adjacent data members. If this is a priority
    // we might need to represent sparse matrices in a different way
    for(unsigned int r=0;r<rows;r++) {
        for(unsigned int c=0;c<other->columns;c++) {
            resultVals[r*other->columns+c] = 0;
            for(unsigned int i=0;i<columns;i++) {
                resultVals[r*other->columns+c] += getElement(r,i)*other->getElement(i,c);
            }
        }
    }
    
    return new OCRDenseMatrix(dest,db);
}
