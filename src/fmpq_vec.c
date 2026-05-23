#include "../include/fmpq_vec.h"

static void fmpq_vec_reserve(fmpq_vec_t *vec, ulong new_cap) {
  if (new_cap <= vec->capacity)
    return;

  vec->data = realloc(vec->data, new_cap * sizeof(fmpq_t));

  for (ulong i = vec->capacity; i < new_cap; i++) {
    fmpq_init(vec->data[i]);
  }

  vec->capacity = new_cap;
}

void fmpq_vec_init(fmpq_vec_t *vec) {
  vec->data = NULL;
  vec->size = 0;
  vec->capacity = 0;
}

void fmpq_vec_clear(fmpq_vec_t *vec) {
  for (ulong i = 0; i < vec->capacity; i++) {
    fmpq_clear(vec->data[i]);
  }
  free(vec->data);
}

void fmpq_vec_push(fmpq_vec_t *vec, const fmpq_t val) {
  if (vec->size + 1 > vec->capacity) {
    ulong new_cap = (vec->capacity == 0) ? 4 : vec->capacity * 2;
    fmpq_vec_reserve(vec, new_cap);
  }

  fmpq_set(vec->data[vec->size], val);
  vec->size++;
}

void fmpq_vec_push_interval(fmpq_vec_t *vec, const fmpq_t a, const fmpq_t b) {
  if (vec->size + 2 > vec->capacity) {
    ulong new_cap = (vec->capacity == 0) ? 4 : vec->capacity * 2;
    while (new_cap < vec->size + 2)
      new_cap *= 2;

    fmpq_vec_reserve(vec, new_cap);
  }

  fmpq_set(vec->data[vec->size], a);
  fmpq_set(vec->data[vec->size + 1], b);
  vec->size += 2;
}

void fmpq_vec_append(fmpq_vec_t *dest, const fmpq_vec_t *src) {
  for (ulong i = 0; i < src->size; i++) {
    fmpq_vec_push(dest, src->data[i]);
  }
}