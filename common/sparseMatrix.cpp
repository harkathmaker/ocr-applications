#include "OCRSparseMatrix.h"
#include "OCRDenseMatrix.h"

OCRSparseMatrix::OCRSparseMatrix(void *_data, ocrGuid_t _datablock): data(_data), datablock(_datablock) {
    cached = false;
    rows = *((unsigned int*)_data);
    columns = *(((unsigned int*)_data)+1);
    capacity = *(((unsigned int*)_data)+2);
    size = ((unsigned int*)_data)+3;
    indices = ((unsigned int*)_data)+4;
    values = (double*)(_data+sizeof(unsigned int)*(4+capacity));
}

void OCRSparseMatrix::cacheEntries() {
    for(unsigned int i=0;i<size;i++) {
        entries.insert(std::pair<unsigned int,double>(indices[i],values[i]));
    }
    cached = true;
}

float OCRSparseMatrix::getElement(unsigned int row, unsigned int col) const {
    assert(row >= 0 && row < getRows() && col >= 0 && col <= getColumns());
    if(!cached) {
        cacheEntries();
    }
    unsigned int key = row * columns + column;
    std::map<unsigned int,double> const_iterator it=entries.find(key);
    if(it != entries.end()) {
        return it->second;
    } else {
        return defaultValue;
    }
}

void OCRSparseMatrix::setElement(unsigned int row, unsigned int col, double val) {
    assert(row >= 0 && row < getRows() && col >= 0 && col <= getColumns());
    if(!cached) {
        cacheEntries();
    }
    unsigned int key = row * columns + column;
    std::map<unsigned int,double> const_iterator it=entries.find(key);
    if(val != defaultValue) {
        if(it != entries.end()) {
            it->second = val;
        } else {
            assert(*size < capacity);
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
            values[it->second] = values[*size];
        }
    }
}
