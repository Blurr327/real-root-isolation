#include "../include/parallel_rri_algo.h"

#include "../include/fmpq_vec.h"
#include "../include/poly_utils.h"

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <omp.h>

fmpz_poly_t og_poly;

void par_subdiv_algo_ext(fmpz_poly_t in_poly, fmpq_vec_t *sol, fmpq_t start,
                         fmpq_t end, int depth) {
  fmpz_t tmp;
  fmpq_t mid;
  fmpq_t tmpq;
  fmpz_poly_t tmp_poly;

  fmpq_init(mid);
  fmpz_init_set_ui(tmp, 1);
  fmpq_init(tmpq);
  fmpz_poly_init(tmp_poly);

  // x -> 1/x
  fmpz_poly_reverse(tmp_poly, in_poly, fmpz_poly_degree(in_poly) + 1);

  // x -> x + 1
  fmpz_poly_taylor_shift(tmp_poly, tmp_poly, tmp);

  int c = count_sign_variations(tmp_poly);

  if (c == 1) {
#pragma omp critical
    fmpq_vec_push_interval(sol, start, end);
  }

  if (c == 0 || c == 1)
    goto cleanup;

  // mid = (start + end) / 2
  fmpq_add(mid, start, end);
  fmpq_div_2exp(mid, mid, 1);

  // check exact root
  fmpz_poly_evaluate_fmpq(tmpq, og_poly, mid);
  if (fmpq_is_zero(tmpq)) {
#pragma omp critical
    fmpq_vec_push_interval(sol, mid, mid);
  }

  // separate polynomials to ensure memory safety
  fmpz_poly_t left_poly, right_poly;
  fmpz_poly_init(left_poly);
  fmpz_poly_init(right_poly);

  // LEFT: x = y/2
  shift_in_proportions_by_k(left_poly, in_poly, -1);

  // RIGHT: x = (y+1)/2
  shift_in_proportions_by_k(right_poly, in_poly, -1);
  fmpz_poly_taylor_shift(right_poly, right_poly, tmp);

  // task creation
#pragma omp task shared(sol) if (depth < MAX_PARALLELIZATION_DEPTH)
  par_subdiv_algo_ext(left_poly, sol, start, mid, depth + 1);

#pragma omp task shared(sol) if (depth < MAX_PARALLELIZATION_DEPTH)
  par_subdiv_algo_ext(right_poly, sol, mid, end, depth + 1);

#pragma omp taskwait

  fmpz_poly_clear(left_poly);
  fmpz_poly_clear(right_poly);

cleanup:
  fmpz_clear(tmp);
  fmpq_clear(mid);
  fmpq_clear(tmpq);
  fmpz_poly_clear(tmp_poly);
}

void par_subdiv_algo(fmpz_poly_t in_poly, fmpq_vec_t *sol) {
  fmpq_t bound;
  fmpq_t start, end;
  fmpz_poly_t tmp_poly;

  fmpz_poly_init(og_poly);

  fmpz_poly_init(tmp_poly);
  fmpq_init(bound);
  fmpq_init(start);
  fmpq_init(end);

  fmpz_poly_set(og_poly, in_poly);

  cauchy_bound(bound, in_poly);

  int b = fmpq_clog(bound, 2);

  fmpq_set_ui(bound, 1, 1);
  if (b > 0)
    fmpq_mul_2exp(bound, bound, b);
  else
    fmpq_div_2exp(bound, bound, b);

  // scale to [0,1]
  shift_in_proportions_by_k(tmp_poly, in_poly, b);

  // entry point for parallelization
#pragma omp parallel
  {
#pragma omp single nowait
    {
      // positive roots
      fmpq_set_ui(start, 0, 1);
      fmpq_set(end, bound);

      par_subdiv_algo_ext(tmp_poly, sol, start, end, 0);

      // negative roots: x -> -x
      neg_varchange(tmp_poly, tmp_poly);

      fmpq_set_ui(start, 0, 1);
      fmpq_mul_si(bound, bound, -1);
      fmpq_set(end, bound);

      par_subdiv_algo_ext(tmp_poly, sol, start, end, 0);
    }
  }

  fmpq_clear(start);
  fmpq_clear(end);
  fmpq_clear(bound);
  fmpz_poly_clear(tmp_poly);
  fmpz_poly_clear(og_poly);
}