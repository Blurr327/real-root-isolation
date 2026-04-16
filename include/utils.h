#ifndef UTILS_H
#define UTILS_H
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>

void p(fmpz_poly_t in, char *s);

void p_irvl(fmpq_t start, fmpq_t end);

#endif