#include "../include/rri_algo.h"
#include "../include/fmpq_vec.h"
#include "../include/poly_utils.h"
#include "../include/utils.h"

#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

int subdiv_algo_ext(fmpz_poly_t in_poly, fmpq_vec_t *sol, fmpq_t start,
                    fmpq_t end) {
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

  if (DEBUG) {
    p_irvl(start, end);
    p(tmp_poly, "tmp for counting : ");
  }

  int c = count_sign_variations(tmp_poly);

  if (c == 1) {
    fmpq_vec_push_interval(sol, start, end);
  }

  if (DEBUG) {
    printf("Number of sign variations counted : %d\n", c);
  }

  if (c == 1 || c == 0)
    goto cleanup;

  // mid = (start + end)/2
  fmpq_add(mid, start, end);
  fmpq_div_2exp(mid, mid, 1);

  // LEFT: x = y/2
  shift_in_proportions_by_k(tmp_poly, in_poly, -1);
  if (DEBUG) {
    p(tmp_poly, "for left : ");
  }
  int c1 = subdiv_algo_ext(tmp_poly, sol, start, mid);

  // RIGHT: x = (y+1)/2
  shift_in_proportions_by_k(tmp_poly, in_poly, -1);
  fmpz_poly_taylor_shift(tmp_poly, tmp_poly, tmp);
  if (DEBUG) {
    p(tmp_poly, "for right : ");
  }
  int c2 = subdiv_algo_ext(tmp_poly, sol, mid, end);

  if ((c1 + c2) % 2 != c % 2) {
    fmpq_vec_push_interval(sol, mid, mid);
  }

cleanup:
  fmpz_clear(tmp);
  fmpq_clear(mid);
  fmpq_clear(tmpq);
  fmpz_poly_clear(tmp_poly);

  return (c % 2);
}

void subdiv_algo(fmpz_poly_t in_poly, fmpq_vec_t *sol) {
  fmpq_t bound;
  fmpq_t start, end;
  fmpz_poly_t tmp_poly;

  fmpz_poly_init(tmp_poly);
  fmpq_init(bound);
  fmpq_init(start);
  fmpq_init(end);

  cauchy_bound(bound, in_poly);

  int b = fmpq_clog(bound, 2);

  fmpq_set_ui(bound, 1, 1);
  if (b > 0)
    fmpq_mul_2exp(bound, bound, b);
  else
    fmpq_div_2exp(bound, bound, b);

  // scale to [0,1]
  shift_in_proportions_by_k(tmp_poly, in_poly, b);

  // positive roots
  fmpq_set_ui(start, 0, 1);
  fmpq_set(end, bound);

  subdiv_algo_ext(tmp_poly, sol, start, end);

  // negative roots: x -> -x
  neg_varchange(tmp_poly, tmp_poly);

  fmpq_set_si(start, 0, 1);
  fmpq_mul_si(bound, bound, -1);
  fmpq_set(end, bound);

  subdiv_algo_ext(tmp_poly, sol, start, end);

  fmpq_clear(start);
  fmpq_clear(end);
  fmpq_clear(bound);
  fmpz_poly_clear(tmp_poly);
}