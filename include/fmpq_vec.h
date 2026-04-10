#ifndef FMPQ_VEC_H
#define FMPQ_VEC_H

#include <flint/fmpq.h>
#include <stdlib.h>

typedef struct {
  fmpq_t *data;
  ulong size;
  ulong capacity;
} fmpq_vec_t;

// init / clear
void fmpq_vec_init(fmpq_vec_t *vec);
void fmpq_vec_clear(fmpq_vec_t *vec);

// push operations
void fmpq_vec_push(fmpq_vec_t *vec, const fmpq_t val);
void fmpq_vec_push_interval(fmpq_vec_t *vec, const fmpq_t a, const fmpq_t b);

#endif