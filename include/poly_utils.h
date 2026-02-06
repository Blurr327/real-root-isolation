//
// Created by tom on 2/5/26.
//

#ifndef POLY_UTILS_H
#define POLY_UTILS_H
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <gmp.h>
int count_sign_variations(fmpz_poly_t poly);

void reverse_coeffs(fmpz_poly_t outPoly, fmpz_poly_t poly);

void shift_in_proportions_by_k(fmpz_poly_t poly, int k);

void cauchy_bound(fmpz_t bound, fmpz_poly_t poly);

void random_dense_fmpz_poly(fmpz_poly_t poly, flint_rand_t state, slong degree, flint_bitcnt_t bits);
#endif //POLY_UTILS_H

