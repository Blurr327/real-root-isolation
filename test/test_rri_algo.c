#include "../include/fmpq_vec.h"
#include "../include/rri_algo.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// --- interval correctness test ---

int test_root_isolation_intervals(fmpz_poly_t test_poly, ulong degree,
                                  fmpq_t *sol, ulong length) {
  int verbose = 1;
  int ret = 1;

  fmpq_t tmp, eval_left, eval_right;
  fmpq_init(tmp);
  fmpq_init(eval_left);
  fmpq_init(eval_right);

  if (((fmpz_poly_degree(test_poly) - length / 2) % 2) == 1) {
    printf("\e[31mNumber of intervals is incoherent.\e[0m\n");
    ret = 0;
    goto cleanup;
  }

  for (ulong i = 0; i < length; i += 2) {
    // evaluate endpoints
    fmpz_poly_evaluate_fmpq(eval_left, test_poly, sol[i]);
    fmpz_poly_evaluate_fmpq(eval_right, test_poly, sol[i + 1]);

    if (verbose) {
      printf("[");
      fmpq_print(sol[i]);
      printf(", ");
      fmpq_print(sol[i + 1]);
      printf("]\n");
    }

    // check sign change
    fmpq_mul(tmp, eval_left, eval_right);

    if (fmpq_sgn(tmp) > 0) {
      printf("\e[31mInvalid interval.\e[0m\n");
      ret = 0;
      goto cleanup;
    }
  }

  printf("\e[32mtest_root_isolation_intervals passed.\e[0m\n");

cleanup:
  fmpq_clear(tmp);
  fmpq_clear(eval_left);
  fmpq_clear(eval_right);

  return ret;
}

// --- test subdiv_algo_ext ---

int test_subdiv_algo_ext(fmpz_poly_t test_poly, ulong degree) {
  int verbose = 1;

  fmpq_vec_t sol;
  fmpq_vec_init(&sol);

  fmpq_t start, end;
  fmpq_init(start);
  fmpq_init(end);

  fmpq_set_ui(start, 0, 1);
  fmpq_set_ui(end, 1, 1);

  subdiv_algo_ext(test_poly, &sol, start, end);

  if (verbose)
    printf("number of intervals : %lu\n", sol.size / 2);

  int r = test_root_isolation_intervals(test_poly, degree, sol.data, sol.size);

  fmpq_vec_clear(&sol);
  fmpq_clear(start);
  fmpq_clear(end);

  return r;
}

// --- test subdiv_algo ---

int test_subdiv_algo(fmpz_poly_t test_poly, ulong degree) {
  int verbose = 1;

  fmpq_vec_t sol;
  fmpq_vec_init(&sol);

  clock_t begin = clock();

  subdiv_algo(test_poly, &sol);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

  printf("\e[36mTime : %lfs\n\e[0m", time_spent);

  if (verbose)
    printf("number of intervals : %lu\n", sol.size / 2);

  int r = test_root_isolation_intervals(test_poly, degree, sol.data, sol.size);

  fmpq_vec_clear(&sol);

  return r;
}

// --- main ---

int main() {
  ulong bits = 8;
  ulong degree = 3;
  int number_of_tests = 10;
  int t = 1;
  int random = 1;
  char test_poly_str[] = "3  -1 -1 1";

  fmpz_poly_t test_poly;
  flint_rand_s randomio;

  flint_randinit(&randomio);
  flint_randseed(&randomio, time(NULL), time(NULL) + 2);

  fmpz_poly_init2(test_poly, degree + 1);

  while (t && (number_of_tests--)) {
    printf("================== TEST NUMBER %d ===============\n",
           number_of_tests + 1);

    if (random)
      fmpz_poly_randtest(test_poly, &randomio, degree + 1, bits);
    else
      fmpz_poly_set_str(test_poly, test_poly_str);

    if (degree < 15) {
      printf("==== TEST POLY ====> ");
      fmpz_poly_print_pretty(test_poly, "x");
      printf("\n");
    }

    if (!fmpz_poly_is_squarefree(test_poly) ||
        !fmpz_poly_get_coeff_ui(test_poly, 0)) {
      printf("\e[31mInvalid polynomial.\n\e[0m");
      continue;
    }

    printf("= test subdiv =\n");
    t = test_subdiv_algo(test_poly, degree);
  }

  fmpz_poly_clear(test_poly);
  flint_randclear(&randomio);

  return 0;
}