//
// Created by tom on 2/5/26.
//

#ifndef RRI_ALGO_H
#define RRI_ALGO_H
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <gmp.h>

// is expected to be of size at least 2 deg(P)
int subdiv_algo_ext(fmpz_poly_t in_poly, fmpq_t sol[], fmpq_t start, fmpq_t end,
                    ulong *next_index_p);

void subdiv_algo(fmpz_poly_t in_poly, fmpq_t sol[], ulong *next_index_p);

#endif // RRI_ALGO_H
