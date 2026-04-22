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

void shift_in_proportions_by_k(fmpz_poly_t outPoly, fmpz_poly_t poly, int k);

void cauchy_bound(fmpq_t bound, fmpz_poly_t poly);

void neg_varchange(fmpz_poly_t out_poly, fmpz_poly_t in_poly);

void random_dense_fmpz_poly(fmpz_poly_t poly, flint_rand_t state, slong degree,
                            flint_bitcnt_t bits);

slong fmpq_clog(fmpq_t n, int b);

void fmpz_poly_split_coeffs(fmpz_poly_t P1, fmpz_poly_t P2,
                            const fmpz_poly_t f);

void fmpz_poly_taylor_shift_split_parallel(fmpz_poly_t res, const fmpz_poly_t f,
                                           const fmpz_t a);
#endif // POLY_UTILS_H
