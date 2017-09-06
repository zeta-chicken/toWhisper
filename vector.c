#include "vector.h"
#include <stdlib.h>
#include <stdio.h>

VECTOR* initVector(size_t sz)
{
    VECTOR* v = (VECTOR*)malloc(sizeof(VECTOR));
    if (v == NULL) {
        fprintf(stderr, "In function, initVector: \
                can not allocate vecotr memory. return null pointer.\n");
        return v;
    }
    v->sz = sz;
    v->elem = calloc(v->sz, sizeof(double));
    return v;
}

int freeVector(VECTOR* v)
{
    if (v == NULL) return 0;
    free(v->elem);
    free(v);
    v = NULL;
    return 0;
}

int resizeVector(VECTOR* v, size_t sz)
{
    double* newElem = realloc(v->elem, sz);
    if (newElem == NULL) {
        freeVector(v);
        fprintf(stderr, "In function, resizeVector: \
                can not change size. free vector memory, done.\n");
        return -1;
    }
    v->sz = sz;
    v->elem = newElem;
    return 0;
}
