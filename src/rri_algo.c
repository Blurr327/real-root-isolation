#include "../include/rri_algo.h"
#include "../include/poly_utils.h"
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <gmp.h>

// sol is expected to be of size at least 2*deg(P)
int subdiv_algo_ext(fmpz_poly_t in_poly, fmpq_t sol[], fmpq_t start, fmpq_t end,
                    ulong *next_index_p) {
  fmpz_t tmp;
  fmpq_t mid;
  fmpz_poly_t tmp_poly;

  // initializations
  fmpq_init(mid);
  fmpz_init_set_ui(tmp, 1);
  fmpz_poly_init(tmp_poly);

  // x = 1/(y + 1) ; roots in ]0, 1[ -> roots in ]0, +inf[
  // x -> 1/x
  fmpz_poly_reverse(tmp_poly, in_poly, fmpz_poly_degree(in_poly) + 1);
  // x -> x + 1
  fmpz_poly_taylor_shift(tmp_poly, tmp_poly, tmp);

  // char *x = "x";
  // printf("tmp_poly : ");
  // fmpz_poly_print_pretty(tmp_poly, x);
  // printf("\n");

  // base case
  int c = count_sign_variations(tmp_poly);

  // printf("counting sign_variations of c : %d\n", c);

  if (c == 1) {
    fmpq_set(sol[*next_index_p], start);
    fmpq_set(sol[*next_index_p + 1], end);
    *next_index_p += 2;
  }

  if (c == 1 || c == 0)
    goto cleanup;

  fmpq_add(mid, start, end);  // mid = start + end
  fmpq_div_2exp(mid, mid, 1); // mid = mid /2
  // x = (1/2) * y variable change ; roots in ]0, 1/2[ -> roots in ]0, 1[
  shift_in_proportions_by_k(tmp_poly, in_poly, -1);

  int n1 = subdiv_algo_ext(tmp_poly, sol, start, mid, next_index_p);

  // x = (y + 1)/2 variable change ; roots in ]1/2, 1[ -> roots in ]0, 1[
  // x -> x / 2
  shift_in_proportions_by_k(tmp_poly, in_poly, -1);
  // x -> x + 1
  fmpz_poly_taylor_shift(tmp_poly, tmp_poly, tmp);

  int n2 = subdiv_algo_ext(tmp_poly, sol, mid, end, next_index_p);

  // checking if 1/2 is a root
  if (((n1 + n2) % 2) != (c % 2)) {
    fmpq_set(sol[*next_index_p], mid);
    fmpq_set(sol[*next_index_p + 1], mid);
    *next_index_p += 2;
  }

cleanup:
  fmpz_clear(tmp);
  fmpq_clear(mid);
  fmpz_poly_clear(tmp_poly);

  return (c % 2);
}

void subdiv_algo(fmpz_poly_t in_poly, fmpq_t sol[], ulong *next_index_p) {
  fmpq_t bound;
  fmpq_t start, end;
  fmpz_poly_t tmp_poly;
  fmpz_poly_init(tmp_poly);
  fmpq_init(bound);
  fmpq_init(start);
  fmpq_init(end);

  cauchy_bound(bound, in_poly);

  // printf("INIT BOUND :");
  // fmpq_print(bound);
  // printf("\n");

  int b = fmpq_clog(bound, 2);
  fmpq_set_ui(bound, 1, 1);
  if (b > 0)
    fmpq_mul_2exp(bound, bound, b);
  else
    fmpq_div_2exp(bound, bound, b);

  // printf("BOUND ROUNDED TO 2^b :");
  // fmpq_print(bound);
  // printf("\n");

  // roots in [0, bound] -> [0, 1]
  shift_in_proportions_by_k(tmp_poly, in_poly, b);
  // printf("tmp_poly after scaling var change ");
  // fmpz_poly_print_pretty(tmp_poly, "x");
  // printf("\n");

  // search in [0, bound]
  fmpq_set_ui(start, 0, 1);
  fmpq_set(end, bound);
  // FINDING POSITIVE ROOTS
  // execute subdiv_algo_ext
  subdiv_algo_ext(tmp_poly, sol, start, end, next_index_p);
  // FINDING NEGATIVE ROOTS
  // x -> -x variable change
  // printf("tmp_poly before neg var change ");
  // fmpz_poly_print_pretty(tmp_poly, "x");
  // printf("\n");

  neg_varchange(tmp_poly, tmp_poly);

  // printf("tmp_poly after neg var change ");
  // fmpz_poly_print_pretty(tmp_poly, "x");
  // printf("\n");
  // search in [-bound, 0]
  fmpq_set_si(start, 0, 1);
  fmpq_mul_si(bound, bound, -1);
  fmpq_set(end, bound);

  subdiv_algo_ext(tmp_poly, sol, start, end, next_index_p);

  // for (int i = 0; i < *next_index_p; i++) {
  //   printf("sol[%d] = ", i);
  //   fmpq_print(sol[i]);
  //   printf("\n");
  // }

  fmpq_clear(start);
  fmpq_clear(end);
  fmpq_clear(bound);
  fmpz_poly_clear(tmp_poly);
}