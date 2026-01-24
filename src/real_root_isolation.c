#include "flint/fmpz_poly.h"
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/random.h>
#include <sys/time.h>
#include <unistd.h>

double wtime() {
  struct timeval timecheck;
  gettimeofday(&timecheck, NULL);
  return (double)timecheck.tv_sec + (double)timecheck.tv_usec / 1000000;
}

void test_flint_taylor_shift() {
  int N = 5, nb_points = 1000;
  FILE *test = fopen("test.csv", "a");
  fmpz_poly_t p, o;
  flint_rand_t r;
  fmpz_t c;

  fmpz_poly_init(p);
  fmpz_poly_init(o);
  fmpz_init(c);
  fmpz_set_ui(c, 1);
  flint_randinit(r);

  int sweep_degree = 0;
  int tau = 512;
  int degree = 5000;

  for (int j = 0; j < nb_points; j++) {
    if (sweep_degree)
      degree = 5 + j;
    else
      tau = 32 + j;
    // generate N random polynomials and time the flint taylor shift function
    double sum_durations = 0;

    for (int i = 0; i < N; i++) {
      fmpz_poly_randtest(p, r, degree + 1, tau);
      double start = wtime();
      fmpz_poly_taylor_shift(o, p, c);
      double end = wtime();
      sum_durations += (end - start);
    }
    double avg = sum_durations / N;

    // write to test.csv : d*tau time (s)
    int xi = degree * tau;
    double yi = avg;
    char x[12], y[65];
    char space = ',', newline = '\n';
    sprintf(x, "%d", xi);
    sprintf(y, "%f", yi);
    fwrite(&x, sizeof(char), strlen(x), test);
    fwrite(&space, sizeof(char), 1, test);
    fwrite(&y, sizeof(char), strlen(y), test);
    fwrite(&newline, sizeof(char), 1, test);
  }
  fclose(test);
}

int main() { test_flint_taylor_shift(); }