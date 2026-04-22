#include "../include/fmpq_vec.h"
#include "../include/poly_utils.h"
#include "../include/rri_algo.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define VERBOSE 1

int test_root_isolation_intervals(fmpz_poly_t test_poly, fmpq_t *sol,
                                  ulong length) {
  int ret = 1;

  fmpq_t tmp, eval_left, eval_right;
  fmpq_init(tmp);
  fmpq_init(eval_left);
  fmpq_init(eval_right);

  if (VERBOSE) {
    for (ulong i = 0; i < length; i += 2) {
      printf("[");
      fmpq_print(sol[i]);
      printf(", ");
      fmpq_print(sol[i + 1]);
      printf("]\n");
    }
  }
  slong degree = fmpz_poly_degree(test_poly);

  if (degree <= 300) {
    ulong num_real_roots = fmpz_poly_num_real_roots(test_poly);

    if (num_real_roots != (length / 2)) {
      printf("\e[31mExpected %lu roots, got %lu roots\e[0m\n", num_real_roots,
             length / 2);
      ret = 0;
      goto cleanup;
    }
  } else {
    if ((degree - (length / 2)) % 2 != 0) {
      printf("\e[31mIncoherent number of roots.\e[0m\n");
      ret = 0;
      goto cleanup;
    }
  }

  for (ulong i = 0; i < length; i += 2) {
    // evaluate endpoints
    fmpz_poly_evaluate_fmpq(eval_left, test_poly, sol[i]);
    fmpz_poly_evaluate_fmpq(eval_right, test_poly, sol[i + 1]);

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

int test_subdiv_algo_ext(fmpz_poly_t test_poly, ulong degree) {

  fmpq_vec_t sol;
  fmpq_vec_init(&sol);

  fmpq_t start, end;
  fmpq_init(start);
  fmpq_init(end);

  fmpq_set_ui(start, 0, 1);
  fmpq_set_ui(end, 1, 1);

  subdiv_algo_ext(test_poly, &sol, start, end);

  if (VERBOSE)
    printf("number of intervals : %lu\n", sol.size / 2);

  int r = test_root_isolation_intervals(test_poly, sol.data, sol.size);

  fmpq_vec_clear(&sol);
  fmpq_clear(start);
  fmpq_clear(end);

  return r;
}

int test_subdiv_algo(fmpz_poly_t test_poly, ulong degree, double *time) {

  fmpq_vec_t sol;
  fmpq_vec_init(&sol);

  clock_t begin = clock();

  subdiv_algo(test_poly, &sol);

  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  *time = time_spent;

  printf("\e[36mTime : %lfs\n\e[0m", time_spent);

  if (VERBOSE)
    printf("number of intervals : %lu\n", sol.size / 2);

  int r = test_root_isolation_intervals(test_poly, sol.data, sol.size);

  fmpq_vec_clear(&sol);

  return r;
}

int main(int argc, char *argv[]) {
  ulong bits = 8;
  ulong degree = 10;
  int number_of_tests = 100;
  int ok = 1;
  int random = 1;
  char test_poly_str[] = "4  -8 0 106 -52";

  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--bits") == 0 && i + 1 < argc) {
      bits = strtoul(argv[++i], NULL, 10);
    } else if (strcmp(argv[i], "--degree") == 0 && i + 1 < argc) {
      degree = strtoul(argv[++i], NULL, 10);
    } else if (strcmp(argv[i], "--tests") == 0 && i + 1 < argc) {
      number_of_tests = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--random") == 0 && i + 1 < argc) {
      random = atoi(argv[++i]);
    } else if (strcmp(argv[i], "--poly") == 0 && i + 1 < argc) {
      strncpy(test_poly_str, argv[++i], sizeof(test_poly_str) - 1);
      test_poly_str[sizeof(test_poly_str) - 1] = '\0';
      random = 0; // force non-random mode
    } else {
      printf("Unknown or incomplete argument: %s\n", argv[i]);
      return 1;
    }
  }

  printf("Configuration:\n");
  printf("  bits   = %lu\n", bits);
  printf("  degree = %lu\n", degree);
  printf("  tests  = %d\n", number_of_tests);
  printf("  random = %d\n", random);
  if (!random)
    printf("  poly   = %s\n", test_poly_str);

  fmpz_poly_t test_poly;
  flint_rand_s randomio;

  flint_randinit(&randomio);
  flint_randseed(&randomio, time(NULL), time(NULL) + 2);

  fmpz_poly_init2(test_poly, degree + 1);

  if (!random)
    number_of_tests = 1;

  double time_spent = 0;
  double avg_time = 0;
  int count = number_of_tests;
  while (ok && (count--)) {
    printf("================== TEST NUMBER %d ===============\n", count + 1);

    // generating/setting the polynomial
    if (!random)
      fmpz_poly_set_str(test_poly, test_poly_str);

    if (random) {
      do
        random_dense_fmpz_poly(test_poly, &randomio, degree, bits);
      while (!fmpz_poly_is_squarefree(test_poly) ||
             !fmpz_poly_get_coeff_ui(test_poly, 0));
    }

    // printing the polynomial
    if (degree < 15) {
      printf("==== TEST POLY ====> ");
      fmpz_poly_print_pretty(test_poly, "x");
      printf("\n");
      fmpz_poly_print(test_poly);
      printf("\n");
    }

    printf("= test subdiv =\n");
    ok = test_subdiv_algo(test_poly, degree, &time_spent);
    avg_time += time_spent;
  }

  avg_time /= number_of_tests;
  printf("\e[36mAverage Time : %lf\e[0m\n", avg_time);
  if (count == -1)
    printf("\e[32mAll tests passed.\e[0m\n");

  // clean up
  fmpz_poly_clear(test_poly);
  flint_randclear(&randomio);

  return 0;
}