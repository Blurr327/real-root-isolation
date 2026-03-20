#include "../include/rri_algo.h"
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <fmpz_extras.h>
#include <time.h>

void rand_poly_knownroots(fmpz_poly_t out_poly, int d, fmpq_t roots[d], int min_root_bound, int max_root_bound) {
  fmpq_t r;
  flint_rand_t s;
  fmpz_poly_t p;
  fmpz_poly_t term;

  fmpz_poly_init(term);
  flint_randinit(s);
  flint_randseed(s, time(NULL), time(NULL) + 2);
  fmpz_init_set_ui(tmp, 1);
  fmpq_init(r);
  fmpz_poly_init(p);
  fmpz_poly_set_coeff_ui(p, 0, 1);
}



int main() {
  int nb_roots = 1;
  int bits = 32;
  fmpz_t tmp;
  fmpq_t r;
  flint_rand_t s;
  fmpz_poly_t p;
  fmpz_poly_t term;

  fmpz_poly_init(term);
  flint_randinit(s);
  flint_randseed(s, time(NULL), time(NULL) + 2);
  fmpz_init_set_ui(tmp, 1);
  fmpq_init(r);
  fmpz_poly_init(p);
  fmpz_poly_set_coeff_ui(p, 0, 1);

  for (int i = 0; i < nb_roots; i++) {
    fmpq_randbits(r, s, bits);

    fmpq_print(r);

    fmpz_t min, max;
    fmpz_max(max, fmpq_numref(r), fmpq_denref(r));
    fmpz_min(min, fmpq_numref(r), fmpq_denref(r));
    fmpz_mul_si(min, min, -1);

    fmpz_poly_set_coeff_fmpz(term, 1, max);
    fmpz_poly_set_coeff_fmpz(term, 0, min);

    fmpz_poly_mul(p, p, term);
  }

  fmpz_poly_print(p);

  fmpz_clear(tmp);
  fmpq_clear(r);
  flint_randclear(s);
}