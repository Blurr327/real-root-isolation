#include "../include/parallel_rri_algo.h"

#include "../include/fmpq_vec.h"
#include "../include/poly_utils.h"

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <omp.h>

void par_subdiv_algo_ext(fmpz_poly_t in_poly, fmpz_poly_t original_poly,
                         fmpq_vec_t *sol, fmpq_t start, fmpq_t end,
                         int depth) {
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
  fmpz_poly_taylor_shift_multi_mod(tmp_poly, tmp_poly, tmp);
  // fmpz_poly_taylor_shift_multi_mod(tmp_poly, tmp_poly, tmp);

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
  fmpz_poly_evaluate_fmpq(tmpq, original_poly, mid);
  if (fmpq_is_zero(tmpq)) {
#pragma omp critical
    fmpq_vec_push_interval(sol, mid, mid);
  }

  if (depth >= MAX_PARALLELIZATION_DEPTH) {
    fmpz_poly_t left_poly;
    fmpz_poly_init(left_poly);

    // LEFT: x = y/2
    shift_in_proportions_by_k(left_poly, in_poly, -1);
    par_subdiv_algo_ext(left_poly, original_poly, sol, start, mid, depth + 1);

    fmpz_poly_clear(left_poly);

    fmpz_poly_t right_poly;
    fmpz_poly_init(right_poly);

    // RIGHT: x = (y+1)/2
    shift_in_proportions_by_k(right_poly, in_poly, -1);
    fmpz_poly_taylor_shift_multi_mod(right_poly, right_poly, tmp);
    // fmpz_poly_taylor_shift_multi_mod(right_poly, right_poly, tmp);
    par_subdiv_algo_ext(right_poly, original_poly, sol, mid, end, depth + 1);

    fmpz_poly_clear(right_poly);

    goto cleanup;
  }

  // task creation
#pragma omp task shared(in_poly, original_poly, sol, start, mid)
  {
    fmpz_poly_t left_poly;
    fmpz_poly_init(left_poly);

    // LEFT: x = y/2
    shift_in_proportions_by_k(left_poly, in_poly, -1);
    par_subdiv_algo_ext(left_poly, original_poly, sol, start, mid, depth + 1);

    fmpz_poly_clear(left_poly);
  }

#pragma omp task shared(in_poly, original_poly, sol, mid, end, tmp)
  {
    fmpz_poly_t right_poly;
    fmpz_poly_init(right_poly);

    // RIGHT: x = (y+1)/2
    shift_in_proportions_by_k(right_poly, in_poly, -1);
    fmpz_poly_taylor_shift_multi_mod(right_poly, right_poly, tmp);
    // fmpz_poly_taylor_shift_multi_mod(right_poly, right_poly, tmp);
    par_subdiv_algo_ext(right_poly, original_poly, sol, mid, end, depth + 1);

    fmpz_poly_clear(right_poly);
  }

#pragma omp taskwait

cleanup:
  fmpz_clear(tmp);
  fmpq_clear(mid);
  fmpq_clear(tmpq);
  fmpz_poly_clear(tmp_poly);
}

void par_subdiv_algo(fmpz_poly_t in_poly, fmpq_vec_t *sol,
                     int omp_num_threads, int flint_num_threads) {
  fmpq_t bound;
  fmpq_t pos_start, pos_end, neg_start, neg_end;
  fmpz_poly_t pos_poly, neg_poly, original_poly;

  if (omp_num_threads < 1)
    omp_num_threads = omp_get_max_threads();
  if (flint_num_threads < 1)
    flint_num_threads = 1;

  omp_set_dynamic(0);
  omp_set_num_threads(omp_num_threads);
  flint_set_num_threads(flint_num_threads);

  fmpz_poly_init(pos_poly);
  fmpz_poly_init(neg_poly);
  fmpz_poly_init(original_poly);
  fmpq_init(bound);
  fmpq_init(pos_start);
  fmpq_init(pos_end);
  fmpq_init(neg_start);
  fmpq_init(neg_end);

  fmpz_poly_set(original_poly, in_poly);

  cauchy_bound(bound, in_poly);

  int b = fmpq_clog(bound, 2);

  fmpq_set_ui(bound, 1, 1);
  if (b > 0)
    fmpq_mul_2exp(bound, bound, b);
  else
    fmpq_div_2exp(bound, bound, b);

  // scale to [0,1]
  shift_in_proportions_by_k(pos_poly, in_poly, b);

  // negative roots: x -> -x
  neg_varchange(neg_poly, pos_poly);

  fmpq_set_ui(pos_start, 0, 1);
  fmpq_set(pos_end, bound);

  fmpq_set_ui(neg_start, 0, 1);
  fmpq_neg(neg_end, bound);

  // entry point for parallelization
#pragma omp parallel num_threads(omp_num_threads)
  {
#pragma omp single nowait
    {
      // positive roots
#pragma omp task shared(pos_poly, original_poly, sol, pos_start, pos_end)
      par_subdiv_algo_ext(pos_poly, original_poly, sol, pos_start, pos_end, 0);

      // negative roots
#pragma omp task shared(neg_poly, original_poly, sol, neg_start, neg_end)
      par_subdiv_algo_ext(neg_poly, original_poly, sol, neg_start, neg_end, 0);

#pragma omp taskwait
    }
  }

  fmpq_clear(pos_start);
  fmpq_clear(pos_end);
  fmpq_clear(neg_start);
  fmpq_clear(neg_end);
  fmpq_clear(bound);
  fmpz_poly_clear(pos_poly);
  fmpz_poly_clear(neg_poly);
  fmpz_poly_clear(original_poly);
}
