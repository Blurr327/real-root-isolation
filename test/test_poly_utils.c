#include "../include/poly_utils.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int test_count_sign_variations_random(int N) {
  flint_rand_t state;
  flint_randinit(state);

  flint_randseed(state, time(NULL), clock());

  fmpz_poly_t poly;
  fmpz_poly_init(poly);

  int ok = 1;

  printf("==== count_sign_variations random test (%d runs) ====\n", N);

  for (int t = 0; t < N; t++) {


    do {
      fmpz_poly_randtest(poly, state, 10, 8);
    } while (fmpz_poly_is_zero(poly));

    int res = count_sign_variations(poly);
    slong deg = fmpz_poly_degree(poly);

    printf("[%d] ", t + 1);
    fmpz_poly_print_pretty(poly, "x");
    printf("  --> variations = %d\n", res);

    if (res < 0 || res > deg) {
      printf("\e[31mFAIL (invalid range)\e[0m\n");
      ok = 0;
      break;
    }
  }

  if (ok)
    printf("\e[32mcount_sign_variations completed\e[0m\n");

  fmpz_poly_clear(poly);
  flint_randclear(state);
  return ok;
}

int test_neg_varchange_random(int N) {
  flint_rand_t state;
  flint_randinit(state);

  fmpz_poly_t poly, transformed;
  fmpz_poly_init(poly);
  fmpz_poly_init(transformed);

  fmpq_t x, lhs, rhs;
  fmpq_init(x);
  fmpq_init(lhs);
  fmpq_init(rhs);

  int ok = 1;

  for (int t = 0; t < N; t++) {
    random_dense_fmpz_poly(poly, state, 8, 8);
    neg_varchange(transformed, poly);

    for (int k = 0; k < 5; k++) {
      // random rational x
      fmpq_set_si(x, (rand() % 5) + 1, (rand() % 5) + 1);

      // lhs = P(-x)
      fmpq_neg(lhs, x);
      fmpz_poly_evaluate_fmpq(lhs, poly, lhs);

      // rhs = transformed(x)
      fmpz_poly_evaluate_fmpq(rhs, transformed, x);

      if (fmpq_cmp(lhs, rhs) != 0) {
        printf("\e[31mFAIL neg_varchange\e[0m\n");
        ok = 0;
        goto cleanup;
      }
    }
  }

  printf("\e[32mneg_varchange passed\e[0m\n");

cleanup:
  fmpz_poly_clear(poly);
  fmpz_poly_clear(transformed);
  fmpq_clear(x);
  fmpq_clear(lhs);
  fmpq_clear(rhs);
  flint_randclear(state);
  return ok;
}

int test_cauchy_bound_random(int N) {
  flint_rand_t state;
  flint_randinit(state);

  fmpz_poly_t poly;
  fmpz_poly_init(poly);

  fmpq_t bound, val;
  fmpq_init(bound);
  fmpq_init(val);

  int ok = 1;

  for (int t = 0; t < N; t++) {
    random_dense_fmpz_poly(poly, state, 6, 6);
    cauchy_bound(bound, poly);

    // brute force integer roots in small range
    for (int r = -10; r <= 10; r++) {
      fmpq_set_si(val, r, 1);

      fmpq_t eval;
      fmpq_init(eval);
      fmpz_poly_evaluate_fmpq(eval, poly, val);

      if (fmpq_is_zero(eval)) {
        // check |r| <= bound
        fmpq_abs(val, val);
        if (fmpq_cmp(val, bound) > 0) {
          printf("\e[31mFAIL cauchy_bound\e[0m\n");
          ok = 0;
          fmpq_clear(eval);
          goto cleanup;
        }
      }
      fmpq_clear(eval);
    }
  }

  printf("\e[32mcauchy_bound passed\e[0m\n");

cleanup:
  fmpz_poly_clear(poly);
  fmpq_clear(bound);
  fmpq_clear(val);
  flint_randclear(state);
  return ok;
}

int test_shift_in_proportions_random(int N) {
  flint_rand_t state;
  flint_randinit(state);

  fmpz_poly_t poly, out;
  fmpz_poly_init(poly);
  fmpz_poly_init(out);

  fmpq_t x, lhs, rhs, tmp;
  fmpq_init(x);
  fmpq_init(lhs);
  fmpq_init(rhs);
  fmpq_init(tmp);

  int ok = 1;

  for (int t = 0; t < N; t++) {
    random_dense_fmpz_poly(poly, state, 6, 6);

    int k = (rand() % 7) - 3; // range [-3, 3], includes negative
    if (k == 0)
      k = 1;

    shift_in_proportions_by_k(out, poly, k);

    slong d = fmpz_poly_degree(poly);

    for (int i = 0; i < 5; i++) {
      // random rational x
      fmpq_set_si(x, (rand() % 5) + 1, (rand() % 5) + 1);

      // compute 2^k * x
      fmpq_set(tmp, x);
      if (k > 0)
        fmpq_mul_2exp(tmp, tmp, k);
      else
        fmpq_div_2exp(tmp, tmp, -k);

      // lhs = P(2^k x)
      fmpz_poly_evaluate_fmpq(lhs, poly, tmp);

      if (k < 0) {
        // multiply lhs by 2^{-k d}
        fmpq_mul_2exp(lhs, lhs, (-k) * d);
      }

      // rhs = out(x)
      fmpz_poly_evaluate_fmpq(rhs, out, x);

      if (fmpq_cmp(lhs, rhs) != 0) {
        printf("\e[31mFAIL shift_in_proportions_by_k\e[0m\n");
        printf("k = %d\n", k);
        printf("poly = ");
        fmpz_poly_print(poly);
        printf("\n");
        printf("out  = ");
        fmpz_poly_print(out);
        printf("\n");

        ok = 0;
        goto cleanup;
      }
    }
  }

  printf("\e[32mshift_in_proportions_by_k passed\e[0m\n");

cleanup:
  fmpz_poly_clear(poly);
  fmpz_poly_clear(out);
  fmpq_clear(x);
  fmpq_clear(lhs);
  fmpq_clear(rhs);
  fmpq_clear(tmp);
  flint_randclear(state);
  return ok;
}

int test_fmpq_clog_random(int N) {
  int ok = 1;

  for (int t = 0; t < N; t++) {
    fmpq_t x;
    fmpq_init(x);

    int num = (rand() % 50) + 1;
    int den = (rand() % 50) + 1;
    int b = 2 + (rand() % 5);

    fmpq_set_si(x, num, den);

    slong r = fmpq_clog(x, b);

    if (r < -100 || r > 100) {
      printf("\e[31mFAIL fmpq_clog (unstable)\e[0m\n");
      ok = 0;
      fmpq_clear(x);
      break;
    }

    fmpq_clear(x);
  }

  if (ok)
    printf("\e[32mfmpq_clog passed\e[0m\n");
  return ok;
}

int main(int argc, char *argv[]) {
  int N = 5;
  if (argc > 1)
    N = atoi(argv[1]);

  srand(time(NULL));

  printf("======= POLY UTILS TESTS (N=%d) =======\n", N);

  int ok = 1;

  ok &= test_count_sign_variations_random(N);
  ok &= test_neg_varchange_random(N);
  ok &= test_cauchy_bound_random(N);
  ok &= test_shift_in_proportions_random(N);
  ok &= test_fmpq_clog_random(N);

  if (ok)
    printf("\n\e[32mALL TESTS PASSED\e[0m\n");
  else
    printf("\n\e[31mSOME TESTS FAILED\e[0m\n");

  return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}