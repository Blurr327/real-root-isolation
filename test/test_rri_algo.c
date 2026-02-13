#include "../include/rri_algo.h"

#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpz.h>
#include <flint/fmpz_poly.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---------- printing ---------- */

static void print_poly(const char *name, const fmpz_poly_t p) {
  printf("%s (deg=%ld): ", name, fmpz_poly_degree(p));
  fmpz_poly_print_pretty(p, "x");
  printf("\n");
}

/* ---------- generate distinct roots in (0,1) ---------- */

static void gen_distinct_root(flint_rand_t state, fmpq_t r, fmpq_t *prev,
                              int prev_n, slong bits) {
  fmpz_t a, b;
  fmpz_init(a);
  fmpz_init(b);

  for (;;) {
    fmpz_randbits(b, state, bits);
    if (fmpz_cmp_ui(b, 2) < 0)
      fmpz_set_ui(b, 2);

    fmpz_randm(a, state, b); /* 0 <= a < b */
    if (fmpz_is_zero(a))
      continue;

    fmpq_set_fmpz_frac(r, a, b); /* canonicalizes automatically */

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

/* ---------- read whole file ---------- */

static char *read_whole_file(const char *path) {
  FILE *f = fopen(path, "rb");
  if (!f)
    return NULL;

  if (fseek(f, 0, SEEK_END) != 0) {
    fclose(f);
    return NULL;
  }
  long n = ftell(f);
  if (n < 0) {
    fclose(f);
    return NULL;
  }
  rewind(f);

  char *buf = (char *)malloc((size_t)n + 1);
  if (!buf) {
    fclose(f);
    return NULL;
  }

  size_t got = fread(buf, 1, (size_t)n, f);
  fclose(f);
  if (got != (size_t)n) {
    free(buf);
    return NULL;
  }

  buf[n] = '\0';
  return buf;
}

/* ---------- main ---------- */

int main(int argc, char **argv) {
  /* random-mode defaults */
  int deg = 20;
  slong bits = 16;

  /* common */
  ulong seed = 0;
  int verbose = 1;
  double start = 0.0, end = 1.0;

  /* input poly mode */
  const char *poly_str = NULL;
  const char *poly_file = NULL;

  for (int i = 1; i < argc; ++i) {
    if (!strcmp(argv[i], "-n") && i + 1 < argc)
      deg = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-bits") && i + 1 < argc)
      bits = atol(argv[++i]);
    else if (!strcmp(argv[i], "-seed") && i + 1 < argc)
      seed = strtoul(argv[++i], NULL, 10);
    else if (!strcmp(argv[i], "-start") && i + 1 < argc)
      start = strtod(argv[++i], NULL);
    else if (!strcmp(argv[i], "-end") && i + 1 < argc)
      end = strtod(argv[++i], NULL);
    else if (!strcmp(argv[i], "-poly") && i + 1 < argc)
      poly_str = argv[++i];
    else if (!strcmp(argv[i], "-polyfile") && i + 1 < argc)
      poly_file = argv[++i];
    else if (!strcmp(argv[i], "-q"))
      verbose = 0;
    else {
      fprintf(
          stderr,
          "Usage:\n"
          "  %s [-n DEG] [-bits BITS] [-seed SEED] [-start A] [-end B] [-q]\n"
          "  %s -poly \"<flint fmpz_poly str>\" [-start A] [-end B] [-q]\n"
          "  %s -polyfile path.txt [-start A] [-end B] [-q]\n\n"
          "Poly string format (fmpz_poly_set_str): \"len c0 c1 ...\"\n"
          "Example: 36*x^2 - 54*x + 20  => \"3 20 -54 36\"\n",
          argv[0], argv[0], argv[0]);
      return 2;
    }
  }

  if (end <= start) {
    fprintf(stderr, "Error: need end > start\n");
    return 2;
  }

  fmpz_poly_t poly, factor;
  fmpz_poly_init(poly);
  fmpz_poly_init(factor);

  int have_known_roots = 0;
  fmpq_t *roots = NULL;

  /* ---------- build/load polynomial ---------- */

  if (poly_str || poly_file) {
    char *owned = NULL;
    const char *s = poly_str;

    if (poly_file) {
      owned = read_whole_file(poly_file);
      if (!owned) {
        fprintf(stderr, "Error: could not read file '%s'\n", poly_file);
        return 2;
      }
      s = owned;
    }
    if (fmpz_poly_set_str(poly, s) != 0) {
      fprintf(stderr, "Error: failed to parse polynomial string. : %s\n", s);
      if (owned)
        free(owned);
      return 2;
    }

    if (owned)
      free(owned);

    if (verbose)
      print_poly("input poly", poly);
    have_known_roots = 0; /* we don't know its roots */
  } else {
    if (deg <= 0) {
      fprintf(stderr, "Degree must be > 0\n");
      return 2;
    }

    flint_rand_t state;
    flint_randinit(state);
    if (seed)
      flint_randseed(state, seed, seed ^ 0x9e3779b97f4a7c15ULL);

    roots = (fmpq_t *)malloc(sizeof(fmpq_t) * (size_t)deg);
    for (int i = 0; i < deg; ++i)
      fmpq_init(roots[i]);

    for (int i = 0; i < deg; ++i)
      gen_distinct_root(state, roots[i], roots, i, bits);

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
      printf("\nExact roots (random mode):\n");
      for (int i = 0; i < deg; ++i) {
        printf("  r[%d] = ", i);
        fmpq_print(roots[i]);
        printf("\n");
      }
    }

    have_known_roots = 1;
    flint_randclear(state);
  }

  /* ---------- run subdivision algorithm ---------- */

  slong d = fmpz_poly_degree(poly);
  if (d < 0) {
    fprintf(stderr, "Error: zero polynomial?\n");
    return 2;
  }

  /* allocate plenty: in your original test you assumed <= deg intervals.
     Here we don't know, so allocate a conservative multiple. */
  int cap_pairs = (int)(2 * (d + 1) + 16);
  double *sol = (double *)malloc(sizeof(double) * (size_t)(2 * cap_pairs));
  for (int i = 0; i < 2 * cap_pairs; ++i)
    sol[i] = -12345.0;

  int next_index = 0;
  int ret = subdiv_algo_ext(poly, sol, (float)start, (float)end, &next_index);

  printf("\nsubdiv_algo_ext returned %d, next_index=%d\n", ret, next_index);

  int interval_count = next_index / 2;
  if (interval_count <= 0) {
    printf("FAIL: no intervals returned\n");
    free(sol);
    if (roots) {
      for (int i = 0; i < deg; ++i)
        fmpq_clear(roots[i]);
      free(roots);
    }
    fmpz_poly_clear(poly);
    fmpz_poly_clear(factor);
    return 1;
  }

  if (verbose) {
    printf("\nReturned intervals (%d):\n", interval_count);
    for (int i = 0; i < interval_count; i++) {
      printf("  I[%d] = [%.17lf, %.17lf]\n", i, sol[2 * i], sol[2 * i + 1]);
    }
  }

  /* ---------- verify containment (only in random mode) ---------- */

  int ok = 1;

  if (have_known_roots) {
    for (int i = 0; i < deg; ++i) {
      double r = fmpq_get_d(roots[i]);

      int found = 0;
      for (int j = 0; j < interval_count; j++) {
        double L = sol[2 * j], U = sol[2 * j + 1];
        if (L > U) {
          double t = L;
          L = U;
          U = t;
        }
        if (L <= r && r <= U) {
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
      printf("\nPASS: all generated roots contained in returned intervals.\n");
    else
      printf("\nFAIL: some generated roots missing from intervals.\n");
  } else {
    printf("\nNOTE: containment check skipped (input polynomial mode has no "
           "known roots).\n");
  }

  /* ---------- cleanup ---------- */

  free(sol);

  if (roots) {
    for (int i = 0; i < deg; ++i)
      fmpq_clear(roots[i]);
    free(roots);
  }

  fmpz_poly_clear(poly);
  fmpz_poly_clear(factor);

  return ok ? 0 : 1;
}
