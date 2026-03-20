#include "../include/rri_algo.h"
#include "../include/poly_utils.h"
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <gmp.h>

// sol is expected to be of size at least 2*deg(P)
int subdiv_algo_ext(fmpz_poly_t in_poly, double sol[], double start, double end,
                    int *next_index_p, fmpz_t tmp, fmpz_poly_t tmp_poly) {
  int next_index = *next_index_p;

  // x = 1/(y + 1) ; roots in ]0, 1[ -> roots in ]0, +inf[
  // x -> 1/x
  reverse_coeffs(tmp_poly, in_poly);
  // x -> x + 1
  fmpz_poly_taylor_shift(tmp_poly, tmp_poly, tmp);

  // base case
  int c = count_sign_variations(tmp_poly);

  if (c == 1) {
    sol[*next_index_p] = start;
    sol[*next_index_p + 1] = end;
    *next_index_p += 2;
  }

  if (c == 1 || c == 0)
    return c;

  float mid = (start + end) / 2.0;
  // x = (1/2) * y variable change ; roots in ]0, 1/2[ -> roots in ]0, 1[
  shift_in_proportions_by_k(tmp_poly, in_poly, -1);

  int n1 =
      subdiv_algo_ext(tmp_poly, sol, start, mid, next_index_p, tmp, tmp_poly);

  // x = (y + 1)/2 variable change ; roots in ]1/2, 1[ -> roots in ]0, 1[
  // x -> x / 2
  shift_in_proportions_by_k(tmp_poly, in_poly, -1);
  // x -> x + 1
  fmpz_poly_taylor_shift(tmp_poly, tmp_poly, tmp);

  int n2 =
      subdiv_algo_ext(tmp_poly, sol, mid, end, next_index_p, tmp, tmp_poly);

  // checking if 1/2 is a root
  if (((n1 + n2) % 2) != (c % 2)) {
    sol[*next_index_p] = mid;
    sol[*next_index_p + 1] = mid;
    *next_index_p += 2;
  }

  return (c % 2);
}

void normalize_polyposroots(fmpz_poly_t out_poly, fmpz_poly_t in_poly) {
  fmpq_t bound;
  fmpq_init(bound);
  cauchy_bound(bound, in_poly);
  // closest 2^b to the bound
  int b = mpq_clog(bound,
                   2); // FIXME : implement function or find an implementation
  // normalize polynomial via variable change ]0;2^b[ => ]0;1[
  shift_in_proportions_by_k(out_poly, in_poly, b);
  fmpq_clear(bound);
}

void subdiv_algo(fmpz_poly_t in_poly, double sol[]) {
  fmpz_t tmp;
  fmpz_poly_t tmp_poly;
  int b;

  fmpz_init_set_ui(tmp, 1);
  fmpz_poly_init(tmp_poly);

  normalize_polyposroots(in_poly, in_poly);

  // execute subdiv_algo_ext
  int next_index = 0;
  subdiv_algo_ext(in_poly, sol, 0, 1, &next_index, tmp, tmp_poly);

  // FINDING NEGATIVE ROOTS
  // x -> -x variable change
  neg_varchange(in_poly,
                in_poly); // FIXME : implement or find implementation in docs
  subdiv_algo_ext(in_poly, sol, 0, 1, &next_index, tmp, tmp_poly);

  fmpz_clear(tmp);
  fmpz_poly_clear(tmp_poly);
}