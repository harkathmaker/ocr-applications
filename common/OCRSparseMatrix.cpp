#include "OCRSparseMatrix.h"
#include "OCRDenseMatrix.h"

OCRSparseMatrix::OCRSparseMatrix(void *_data, double defVal, ocrGuid_t _datablock): defaultValue(defVal) {
    data = _data;
    datablock = _datablock;
    cached = false;
    rows = *((unsigned int*)_data);
    columns = *(((unsigned int*)_data)+1);
    capacity = *(((unsigned int*)_data)+2);
    size = ((unsigned int*)_data)+3;
    indices = ((unsigned int*)_data)+4;
    values = (double*)(_data+sizeof(unsigned int)*(4+capacity));
}

OCRSparseMatrix::OCRSparseMatrix(unsigned int rs, unsigned int cols, unsigned int cap): capacity(cap), defaultValue(0.0) {
    rows = rs;
    columns = cols;
    cached = true;
    DBCREATE(&datablock, &data, sizeof(unsigned int)*(4+capacity) + sizeof(double)*capacity, 0,
             NULL_GUID, NO_ALLOC);
    *((unsigned int*)data) = rows;
    *(((unsigned int*)data)+1) = columns;
    *(((unsigned int*)data)+2) = capacity;
    size = ((unsigned int*)data)+3;
    *size = 0;
    indices = ((unsigned int*)data)+4;
    values = (double*)(data+sizeof(unsigned int)*(4+capacity));
}

OCRSparseMatrix::OCRSparseMatrix(OCRDenseMatrix *m, double defVal): defaultValue(defVal) {
    rows = m->getRows();
    columns = m->getColumns();
    for(unsigned int c=0;c<columns;c++) {
        for(unsigned int r=0;r<rows;r++) {
            if(m->getElement(r,c) != defVal) {
                entries.insert(std::pair<unsigned int,double>(r*columns+c,m->getElement(r,c)));
            }
        }
    }
    capacity = entries.size();
    DBCREATE(&datablock, &data, sizeof(unsigned int)*(4+capacity) + sizeof(double)*capacity, 0,
             NULL_GUID, NO_ALLOC);
    *((unsigned int*)data) = rows;
    *(((unsigned int*)data)+1) = columns;
    *(((unsigned int*)data)+2) = capacity;
    size = ((unsigned int*)data)+3;
    *size = capacity;
    indices = ((unsigned int*)data)+4;
    values = (double*)(data+sizeof(unsigned int)*(4+capacity));

    unsigned int k = 0;
    for(std::map<unsigned int,double>::iterator it=entries.begin();it!=entries.end();it++) {
        indices[k] = it->first;
        values[k] = it->second;
        k++;
    }
    cached = true;
}

void OCRSparseMatrix::cacheEntries() const {
    for(unsigned int i=0;i<*size;i++) {
        entries.insert(std::pair<unsigned int,double>(indices[i],values[i]));
    }
    cached = true;
}

double OCRSparseMatrix::getElement(unsigned int row, unsigned int col) const {
    ASSERT(row >= 0 && row < getRows() && col >= 0 && col <= getColumns());
    if(!cached) {
        cacheEntries();
    }
    unsigned int key = row * columns + col;
    std::map<unsigned int,double>::const_iterator it=entries.find(key);
    if(it != entries.end()) {
        return it->second;
    } else {
        return defaultValue;
    }
}

void OCRSparseMatrix::setElement(unsigned int row, unsigned int col, double val) {
    ASSERT(row >= 0 && row < getRows() && col >= 0 && col <= getColumns());
    if(!cached) {
        cacheEntries();
    }
    unsigned int key = row * columns + col;
    std::map<unsigned int,double>::iterator it=entries.find(key);
    if(val != defaultValue) {
        if(it != entries.end()) {
            it->second = val;
        } else {
            ASSERT(*size < capacity);
            entries.insert(std::pair<unsigned int,double>(key,val));
            indices[*size] = key;
            values[*size] = val;
            *size += 1;
        }
    } else {
        if(it != entries.end()) {
            // Remove entry since it is now default value
            *size -= 1;
            indices[it->first] = indices[*size];
            values[it->first] = values[*size];
        }
    }
}
