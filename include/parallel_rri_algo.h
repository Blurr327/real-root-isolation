#ifndef PARALLEL_RRI_ALGO_H
#define PARALLEL_RRI_ALGO_H
#include "fmpq_vec.h"
#include <flint/fmpz_poly.h>

#define DEBUG 0
#define MAX_PARALLELIZATION_DEPTH 16

void par_subdiv_algo_ext(fmpz_poly_t in_poly, fmpq_vec_t *sol, fmpq_t start,
                         fmpq_t end, int depth);

void par_subdiv_algo(fmpz_poly_t in_poly, fmpq_vec_t *sol);

#endif