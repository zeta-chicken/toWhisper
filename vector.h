#ifndef _VECTOR_H_20170809_
#define _VECTOR_H_20170809_

#include <stddef.h>

typedef struct _vector {
    double* elem;
    size_t sz;
} VECTOR;

VECTOR* initVector(size_t);
int freeVector(VECTOR*);
int resizeVector(VECTOR*, size_t);

#endif
