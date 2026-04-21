#include "../include/fmpq_vec.h"
#include "../include/parallel_rri_algo.h"
#include "../include/poly_utils.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define VERBOSE 1

int test_root_isolation_intervals(fmpz_poly_t test_poly, fmpq_t *sol,
                                  ulong length) {
  int ret = 1;

  fmpq_t tmp, eval_left, eval_right;
  fmpq_init(tmp);
  fmpq_init(eval_left);
  fmpq_init(eval_right);

  ulong num_real_roots = fmpz_poly_num_real_roots(test_poly);

  if (num_real_roots != (length / 2)) {
    printf("\e[31mExpected %lu roots, got %lu roots\e[0m\n", num_real_roots,
           length / 2);
    ret = 0;
    goto cleanup;
  }

  for (ulong i = 0; i < length; i += 2) {
    // evaluate endpoints
    fmpz_poly_evaluate_fmpq(eval_left, test_poly, sol[i]);
    fmpz_poly_evaluate_fmpq(eval_right, test_poly, sol[i + 1]);

    if (VERBOSE) {
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

// int test_par_subdiv_algo_ext(fmpz_poly_t test_poly, ulong degree) {

//   fmpq_vec_t sol;
//   fmpq_vec_init(&sol);

//   fmpq_t start, end;
//   fmpq_init(start);
//   fmpq_init(end);

//   fmpq_set_ui(start, 0, 1);
//   fmpq_set_ui(end, 1, 1);

//   par_subdiv_algo_ext(test_poly, &sol, start, end, 0);

//   if (VERBOSE)
//     printf("number of intervals : %lu\n", sol.size / 2);

//   int r = test_root_isolation_intervals(test_poly, sol.data, sol.size);

//   fmpq_vec_clear(&sol);
//   fmpq_clear(start);
//   fmpq_clear(end);

//   return r;
// }

int test_par_subdiv_algo(fmpz_poly_t test_poly, ulong degree) {

  fmpq_vec_t sol;
  fmpq_vec_init(&sol);

  double begin = omp_get_wtime();

  par_subdiv_algo(test_poly, &sol);

  double end = omp_get_wtime();

  printf("\e[36mTime : %lfs\n\e[0m", end - begin);

  if (VERBOSE)
    printf("number of intervals : %lu\n", sol.size / 2);

  int r = test_root_isolation_intervals(test_poly, sol.data, sol.size);

  fmpq_vec_clear(&sol);

  return r;
}

int main() {
  ulong bits = 8;
  ulong degree = 5000;
  int number_of_tests = 1;
  int ok = 1;
  int random = 1;
  char test_poly_str[] = "4  -8 0 106 -52";

  fmpz_poly_t test_poly;
  flint_rand_s randomio;

  flint_randinit(&randomio);
  flint_randseed(&randomio, time(NULL), time(NULL) + 2);

  fmpz_poly_init2(test_poly, degree + 1);

  if (!random)
    number_of_tests = 1;

  while (ok && (number_of_tests--)) {
    printf("================== TEST NUMBER %d ===============\n",
           number_of_tests + 1);

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
    ok = test_par_subdiv_algo(test_poly, degree);
  }

  if (number_of_tests == -1)
    printf("\e[32mAll tests passed.\e[0m\n");

  // clean up
  fmpz_poly_clear(test_poly);
  flint_randclear(&randomio);

  return 0;
}