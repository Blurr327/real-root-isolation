#include "../include/utils.h"

void p(fmpz_poly_t in, char *s) {
  printf("%s\n", s);
  fmpz_poly_print_pretty(in, "x");
  printf("\n");
}

void p_irvl(fmpq_t start, fmpq_t end) {
  printf("[");
  fmpq_print(start);
  printf(",");
  fmpq_print(end);
  printf("]\n");
}