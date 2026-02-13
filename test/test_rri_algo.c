#include "../include/rri_algo.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_poly(const char *name, const fmpz_poly_t p) {
  printf("%s (deg=%ld): ", name, fmpz_poly_degree(p));
  fmpz_poly_print_pretty(p, "x");
  printf("\n");
}

/* ---------- interval helpers ---------- */

typedef struct {
  double l, u;
} interval_t;

static int cmp_interval(const void *p, const void *q) {
  const interval_t *a = p;
  const interval_t *b = q;
  return (a->l > b->l) - (a->l < b->l);
}

/* ---------- generate distinct roots in (0,1) ---------- */

static void gen_distinct_root(flint_rand_t state, fmpq_t r, fmpq_t *prev,
                              int prev_n, slong bits) {
  fmpz_t a, b;
  fmpz_init(a);
  fmpz_init(b);

  for (;;) {
    /* choose b >= 2 */
    fmpz_randbits(b, state, bits);
    if (fmpz_cmp_ui(b, 2) < 0)
      fmpz_set_ui(b, 2);

    /* choose 0 < a < b */
    fmpz_randm(a, state, b);
    if (fmpz_is_zero(a))
      continue;

    fmpq_set_fmpz_frac(r, a, b); /* canonicalizes automatically */

    /* ensure uniqueness */
    int dup = 0;
    for (int i = 0; i < prev_n; ++i) {
      if (fmpq_equal(r, prev[i])) {
        dup = 1;
        break;
      }
    }
    if (!dup)
      break;
  }

  fmpz_clear(a);
  fmpz_clear(b);
}

/* ---------- main randomized test ---------- */

int main(int argc, char **argv) {
  int deg = 20;    /* default degree */
  slong bits = 16; /* denominator size */
  ulong seed = 0;
  int verbose = 1;

  /* CLI parsing */
  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "-n") && i + 1 < argc)
      deg = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-bits") && i + 1 < argc)
      bits = atol(argv[++i]);
    else if (!strcmp(argv[i], "-seed") && i + 1 < argc)
      seed = strtoul(argv[++i], NULL, 10);
    else if (!strcmp(argv[i], "-q"))
      verbose = 0;
    else {
      fprintf(stderr, "Usage: %s [-n DEG] [-bits BITS] [-seed SEED] [-q]\n",
              argv[0]);
      return 2;
    }
  }

  if (deg <= 0) {
    fprintf(stderr, "Degree must be > 0\n");
    return 2;
  }

  /* RNG */
  flint_rand_t state;
  flint_randinit(state);
  if (seed)
    flint_randseed(state, seed, seed ^ 0x9e3779b97f4a7c15ULL);

  /* ---------- generate distinct rational roots ---------- */

  fmpq_t *roots = malloc(sizeof(fmpq_t) * deg);
  for (int i = 0; i < deg; ++i)
    fmpq_init(roots[i]);

  for (int i = 0; i < deg; ++i)
    gen_distinct_root(state, roots[i], roots, i, bits);

  /* ---------- build polynomial Π (b x − a) ---------- */

  fmpz_poly_t poly, factor;
  fmpz_poly_init(poly);
  fmpz_poly_init(factor);
  fmpz_poly_one(poly);

  for (int i = 0; i < deg; ++i) {
    const fmpz *a = fmpq_numref(roots[i]);
    const fmpz *b = fmpq_denref(roots[i]);

    fmpz_poly_zero(factor);
    fmpz_poly_set_coeff_fmpz(factor, 1, b);

    fmpz_t neg_a;
    fmpz_init(neg_a);
    fmpz_neg(neg_a, a);
    fmpz_poly_set_coeff_fmpz(factor, 0, neg_a);
    fmpz_clear(neg_a);

    fmpz_poly_mul(poly, poly, factor);
  }

  if (verbose) {
    print_poly("test poly", poly);
    printf("\nExact roots:\n");
    for (int i = 0; i < deg; ++i) {
      printf("  r[%d] = ", i);
      fmpq_print(roots[i]);
      printf("\n");
    }
  }

  /* ---------- run subdivision algorithm ---------- */

  double *sol = malloc(sizeof(double) * 2 * deg);
  for (int i = 0; i < 2 * deg; ++i)
    sol[i] = -12345.0f;

  int next_index = 0;
  int ret = subdiv_algo_ext(poly, sol, 0.0f, 1.0f, &next_index);

  printf("\nsubdiv_algo_ext returned %d, next_index=%d\n", ret, next_index);

  int interval_count = next_index / 2;

  if (interval_count <= 0) {
    printf("FAIL: no intervals returned\n");
    return 1;
  }

  if (verbose) {
    printf("\nReturned intervals (%d):\n", interval_count);
    for (int i = 0; i < interval_count; i++) {
      printf("  I[%d] = [%.17lf, %.17lf]\n", i, sol[2 * i], sol[2 * i + 1]);
    }
  }

  /* ---------- verify containment ---------- */

  int ok = 1;

  for (int i = 0; i < deg; ++i) {
    double r = fmpq_get_d(roots[i]);

    int found = 0;
    for (int j = 0; j < interval_count; j++) {
      if (sol[2 * j] <= r && r <= sol[2 * j + 1]) {
        found = 1;
        break;
      }
    }

    if (!found) {
      printf("FAIL: root %.17g not contained in any interval\n", r);
      ok = 0;
    }
  }

  if (ok)
    printf("\nPASS: all roots contained in returned intervals.\n");
  else
    printf("\nFAIL: some roots missing from intervals.\n");

  /* ---------- cleanup ---------- */

  free(sol);

  for (int i = 0; i < deg; ++i)
    fmpq_clear(roots[i]);
  free(roots);

  fmpz_poly_clear(poly);
  fmpz_poly_clear(factor);

  flint_randclear(state);

  return ok ? 0 : 1;
}
