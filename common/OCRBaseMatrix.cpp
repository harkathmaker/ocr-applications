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

double OCRMatrix::getDeterminant() const {
    ASSERT(columns == rows);
    double determinant = 1.0;
    // Store a dense copy of the matrix for calculating the determinant. NOT EFFICIENT,
    // use only at the uttermost end of need.
    double *copy;
    ocrGuid_t copyDb;
    DBCREATE(&copyDb, (void**)&copy,sizeof(double)*rows*columns,0,NULL_GUID,NO_ALLOC);
    for(unsigned int r=0;r<rows;r++) {
        for(unsigned int c=0;c<columns;c++) {
            copy[r*columns+c] = getElement(r,c);
        }
    }
    
    // Access rows through pointer array to simplify swapping.
    double *copyRows[rows];
    for(unsigned int i=0;i<rows;i++) {
        copyRows[i] = copy+columns*i;
    }

    // For every column j:
    for(unsigned int j=0;j<columns;j++) {
        // Find a row m with a nonzero value in the column. If none exist, det = 0
        unsigned int m = j;
        while(m < rows && copyRows[m][j] == 0.0) {
            m++;
        }
        if(m == rows) {
            return 0.0;
        }
        // Swap the rows if not the current one, j. Multiply det by -1
        if(m != j) {
            std::swap(copyRows[m],copyRows[j]);
            determinant *= -1;
        }
        // For every row k below j:
        for(unsigned int k=j+1;k<rows;k++) {
            // Subtract multiples of j such that M_k,j == 0
            // Only values >= j need be subtracted, since those to the left are zero
            if(copyRows[k][j] == 0.0) {
                continue;
            }
            double mult = copyRows[k][j] / copyRows[j][j];
            for(unsigned int e=j;e<columns;e++) {
                copyRows[k][e] -= copyRows[j][e] * mult;
            }
        }
    }

    // Determinant = product of diagonals
    for(unsigned int i=0;i<columns;i++) {
        determinant *= copyRows[i][i];
    }

    ocrDbDestroy(copyDb);

    return determinant;
}

OCRMatrix *OCRMatrix::getInverse() const {
    ASSERT(columns == rows);

    // Store a dense copy of the matrix for calculating the determinant. NOT EFFICIENT,
    // use only at the uttermost end of need.
    double *copy;
    ocrGuid_t copyDb;
    DBCREATE(&copyDb, (void**)&copy,sizeof(double)*rows*columns*2,0,NULL_GUID,NO_ALLOC);
    for(unsigned int r=0;r<rows;r++) {
        for(unsigned int c=0;c<columns;c++) {
            copy[r*columns*2+c] = getElement(r,c);
            copy[r*columns*2+c+columns] = (r == c ? 1 : 0);
        }
    }
    
    // Access rows through pointer array to simplify swapping.
    double *invRows[rows];
    for(unsigned int i=0;i<rows;i++) {
        invRows[i] = copy+columns*2*i;
    }
    // For every column j:
    for(unsigned int j=0;j<columns;j++) {
        unsigned int m = j;
        // Find a row m w/ nonzero value. If none exist, no inverse exists.
        while(m < rows && invRows[m][j] == 0) {
            m++;
        }
        ASSERT(m != rows);
        // Swap row into position if necessary
        if(m != j) {
            std::swap(invRows[m],invRows[j]);
        }
        
        // Normalize row
        double normalizeFactor = 1.0 / invRows[j][j];
        for(unsigned int e=j;e<columns*2;e++) {
            invRows[j][e] *= normalizeFactor;
        }

        // for every row k != j:
        for(unsigned int k=0;k<rows;k++) {
            if(k == j) {
                continue;
            }
            // subtract multiples of row j such that M_k,j == 0
            // Only values >= j need be subtracted.
            normalizeFactor = -invRows[k][j];
            for(unsigned int e=j;e<columns*2;e++) {
                invRows[k][e] += invRows[j][e] * normalizeFactor;
            }
        }
    }

    // Build matrix from augmented right side
    OCRDenseMatrix *inverse = new OCRDenseMatrix(rows,columns);
    for(unsigned int r=0;r<rows;r++) {
        for(unsigned int c=columns;c<columns*2;c++) {
            inverse->setElement(r,c-columns,invRows[r][c]);
        }
    }

    ocrDbDestroy(copyDb);

    return inverse;
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
