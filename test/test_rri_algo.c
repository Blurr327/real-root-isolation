#include "../include/rri_algo.h"
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz_poly.h>
#include <time.h>

int test_root_isolation_intervals(fmpz_poly_t test_poly, ulong degree,
                                  fmpq_t *sol, ulong length) {
  int verbose = 1;
  fmpq_t tmp;
  fmpq_t *evals = malloc((2 * (degree)) * sizeof(fmpq_t));
  fmpq_init(tmp);
  for (int i = 0; i < (2 * degree); i++)
    fmpq_init(evals[i]);

  for (int i = 0; i < length; i++) {
    fmpz_poly_evaluate_fmpq(evals[i], test_poly, sol[i]);
    if (i % 2 == 1) {
      if (verbose) {
        printf("[");
        fmpq_print(sol[i - 1]);
        printf(", ");
        fmpq_print(sol[i]);
        printf("],\n");
      }
      fmpq_mul(tmp, evals[i - 1], evals[i]);
      if (fmpq_sgn(tmp) > 0) {
        printf("\e[31mtest_root_isolation_intervals failed. Detected too many "
               "or no "
               "roots in one "
               "interval.\e[0m\n");
        goto cleanup;
      }
    }
  }
  printf("\e[32mtest_root_isolation_intervals passed.\e[0m\n");

cleanup:
  for (int i = 0; i < (2 * degree); i++)
    fmpq_clear(evals[i]);
  free(evals);
  fmpq_clear(tmp);

  return 1;
}

// test whether it detects all roots in ]0, 1[
int test_subdiv_algo_ext(fmpz_poly_t test_poly, ulong degree) {
  int verbose = 1;
  fmpq_t *sol = malloc((2 * (degree)) * sizeof(fmpq_t));
  fmpq_t start, end;
  for (int i = 0; i < (2 * degree); i++)
    fmpq_init(sol[i]);
  fmpq_init(start);
  fmpq_init(end);

  fmpq_set_ui(start, 0, 1);
  fmpq_set_ui(end, 1, 1);

  ulong next_index = 0;
  subdiv_algo_ext(test_poly, sol, start, end, &next_index);
  if (verbose)
    printf("number of intervals : %lu\n", next_index / 2);

  int r = test_root_isolation_intervals(test_poly, degree, sol, next_index);

  fmpq_clear(start);
  fmpq_clear(end);
  for (int i = 0; i < (2 * degree); i++)
    fmpq_clear(sol[i]);
  free(sol);
  return r;
}

int test_subdiv_algo(fmpz_poly_t test_poly, ulong degree) {
  int verbose = 1;
  fmpq_t *sol = malloc((2 * (degree)) * sizeof(fmpq_t));
  for (int i = 0; i < (2 * degree); i++)
    fmpq_init(sol[i]);

  ulong next_index = 0;
  subdiv_algo(test_poly, sol, &next_index);
  if (verbose)
    printf("number of intervals : %lu\n", next_index / 2);

  int r = test_root_isolation_intervals(test_poly, degree, sol, next_index);

  for (int i = 0; i < (2 * degree); i++)
    fmpq_clear(sol[i]);
  free(sol);
  return r;
}

int main() {
  ulong bits = 8;
  ulong degree = 20;
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
    printf("==== TEST POLY ====> ");
    fmpz_poly_print_pretty(test_poly, "x");
    printf("\n");

    if (!fmpz_poly_is_squarefree(test_poly)) {
      printf("Given polynomial is not square free.\n");
      continue;
    }

    // printf("test subdiv ext \n");
    // t = test_subdiv_algo_ext(test_poly, degree);
    printf("= test subdiv =\n");
    t = test_subdiv_algo(test_poly, degree);
  }

  fmpz_poly_clear(test_poly);
  flint_randclear(&randomio);
}