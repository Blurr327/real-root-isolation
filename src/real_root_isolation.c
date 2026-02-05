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
  int N = 10, nb_points = 150;
  FILE *test = fopen("testxyz.csv", "w");
  fmpz_poly_t p, o;
  flint_rand_t r;
  fmpz_t c;

  fmpz_poly_init(p);
  fmpz_poly_init(o);
  fmpz_init(c);
  fmpz_set_ui(c, 1);
  flint_randinit(r);

  int sweep_degree = 1;
  int tau = 512;
  int degree = 5000;

  for (int j = 0; j < nb_points; j++) {
    degree = 15 + j;
    for (int k = 0; k < nb_points; k++) {
      tau = 1000 + k*100;
      // generate N random polynomials and time the flint taylor shift function
      double sum_durations = 0;

      for (int i = 0; i < N; i++) {
        fmpz_poly_randtest(p, r, degree + 1, tau);
        double start = wtime();
        fmpz_poly_taylor_shift_divconquer(o, p, c);
        double end = wtime();
        sum_durations += (end - start);
      }
      double avg = sum_durations / N;

      // write to test.csv : d*tau time (s)
      int xi = degree * tau;
      double yi = avg;
      char x[12], y[65], xd[12], zd[12];
      char space = ',', newline = '\n';
      sprintf(x, "%d", xi);
      sprintf(y, "%f", yi);
      sprintf(xd, "%d", degree);
      sprintf(zd, "%d", tau);

      fwrite(&xd, sizeof(char), strlen(xd), test);
      fwrite(&space, sizeof(char), 1, test);
      fwrite(&zd, sizeof(char), strlen(zd), test);
      fwrite(&space, sizeof(char), 1, test);
      fwrite(&y, sizeof(char), strlen(y), test);
      fwrite(&newline, sizeof(char), 1, test);
      /*
      fwrite(&x, sizeof(char), strlen(x), test);
      fwrite(&space, sizeof(char), 1, test);
      fwrite(&y, sizeof(char), strlen(y), test);
      fwrite(&newline, sizeof(char), 1, test);
      */
    }
    }


  fclose(test);
}



int main() { test_flint_taylor_shift(); }