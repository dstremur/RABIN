#ifndef BENCH_PLOT_H
#define BENCH_PLOT_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef struct {
  const char* name;   /* e.g., "BPSW", "QFT", "MR-Mont" */
  const char* color;  /* optional ANSI code, e.g., "\x1b[36m"; NULL for auto */
  double value;       /* latency (s) or throughput (tests/sec) */
} bench_series_t;

typedef struct {
  const char* group;          /* e.g., "100b", "256b", "1000b" */
  const bench_series_t* series;
  size_t num_series;
} bench_group_t;

typedef enum {
  BENCH_LOWER_IS_BETTER,
  BENCH_HIGHER_IS_BETTER
} bench_direction_t;

/* Callback signature: run `iters` reps on `ctx` and return an accumulator/result
 * so return values can be folded into a sink and avoid dead-code elimination. */
typedef uint64_t (*bench_op_fn)(void* ctx, size_t iters);

typedef struct {
  const char* name;   /* e.g., "BPSW", "QFT" */
  const char* color;  /* optional ANSI color; NULL for auto */
  bench_op_fn fn;
} bench_algo_t;

typedef struct {
  int bits;
  void* ctx;          /* per-case input passed to each bench_op_fn */
} bench_case_t;

static volatile uint64_t bench_sink = 0;

static const char* bench_auto_color(size_t idx)
{
  static const char* palette[] = {
      "\x1b[36m", /* cyan */
      "\x1b[33m", /* amber/yellow */
      "\x1b[32m", /* green */
      "\x1b[35m", /* magenta */
      "\x1b[34m", /* blue */
      "\x1b[96m", /* bright cyan */
  };
  return palette[idx % (sizeof(palette) / sizeof(palette[0]))];
}

static void bench_print_bar(double value, double max_val, int max_cols,
                            const char* color, bool use_color)
{
  static const char* eighths[9] = {
      "", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"
  };

  if (max_val <= 0.0 || value <= 0.0) {
    if (use_color && color) printf("%s▏\x1b[0m", color);
    else printf("▏");
    for (int c = 1; c < max_cols; c++) putchar(' ');
    return;
  }

  double scaled = (value / max_val) * (double)max_cols;
  if (scaled < 0.0) scaled = 0.0;
  if (scaled > (double)max_cols) scaled = (double)max_cols;

  long total_eighths = llround(scaled * 8.0);
  long max_eighths = (long)max_cols * 8L;
  if (total_eighths > max_eighths) total_eighths = max_eighths;

  int full = (int)(total_eighths / 8L);
  int rem = (int)(total_eighths % 8L);

  if (use_color && color) printf("%s", color);
  for (int c = 0; c < full; c++) printf("█");
  if (rem > 0) printf("%s", eighths[rem]);
  if (use_color && color) printf("\x1b[0m");

  int printed = full + (rem > 0 ? 1 : 0);
  for (int c = printed; c < max_cols; c++) putchar(' ');
}

static double bench_best_in_group(const bench_group_t* g,
                                  bench_direction_t dir)
{
  double best = (dir == BENCH_LOWER_IS_BETTER) ? INFINITY : -INFINITY;
  for (size_t s = 0; s < g->num_series; s++) {
    double v = g->series[s].value;
    if (v <= 0.0) continue;
    if (dir == BENCH_LOWER_IS_BETTER && v < best) best = v;
    if (dir == BENCH_HIGHER_IS_BETTER && v > best) best = v;
  }
  return isfinite(best) ? best : 0.0;
}

static void bench_plot_metric(const char* title, bench_direction_t dir,
                              const bench_group_t* groups, size_t num_groups,
                              const char* unit, int decimals)
{
  const int bar_width = 40;
  bool use_color = isatty(STDOUT_FILENO) != 0;

  double max_val = 0.0;
  size_t max_group_len = 5;
  size_t max_name_len = 4;

  for (size_t i = 0; i < num_groups; i++) {
    size_t glen = groups[i].group ? strlen(groups[i].group) : 0;
    if (glen > max_group_len) max_group_len = glen;
    for (size_t s = 0; s < groups[i].num_series; s++) {
      double v = groups[i].series[s].value;
      if (v > max_val) max_val = v;
      size_t nlen =
          groups[i].series[s].name ? strlen(groups[i].series[s].name) : 0;
      if (nlen > max_name_len) max_name_len = nlen;
    }
  }

  const char* dir_note =
      (dir == BENCH_LOWER_IS_BETTER) ? "lower is better" : "higher is better";

  printf("\n%s [%s]\n", title, dir_note);

  for (size_t i = 0; i < num_groups; i++) {
    if (groups[i].num_series > 0) {
      for (size_t s = 0; s < groups[i].num_series; s++) {
        const char* col = groups[i].series[s].color
                              ? groups[i].series[s].color
                              : bench_auto_color(s);
        if (use_color) printf("  %s■\x1b[0m %s", col, groups[i].series[s].name);
        else printf("  [#] %s", groups[i].series[s].name);
      }
      printf("\n");
      break;
    }
  }

  int prefix_width = (int)(max_group_len + max_name_len + 3);

  for (size_t i = 0; i < num_groups; i++) {
    const bench_group_t* g = &groups[i];
    double best = bench_best_in_group(g, dir);

    for (size_t s = 0; s < g->num_series; s++) {
      const bench_series_t* item = &g->series[s];
      const char* col = item->color ? item->color : bench_auto_color(s);

      if (s == 0) {
        printf("%-*s %-*s |", (int)max_group_len,
               g->group ? g->group : "", (int)max_name_len,
               item->name ? item->name : "");
      } else {
        printf("%*s %-*s |", (int)max_group_len, "", (int)max_name_len,
               item->name ? item->name : "");
      }

      bench_print_bar(item->value, max_val, bar_width, col, use_color);

      printf("| %.*f%s", decimals, item->value, unit ? unit : "");

      if (best > 0.0 && item->value > 0.0) {
        double ratio = (dir == BENCH_LOWER_IS_BETTER)
                           ? (item->value / best)
                           : (best / item->value);
        if (ratio < 1.05) printf("  (best)");
        else printf("  (%.1fx slower)", ratio);
      }
      printf("\n");
    }
  }

  printf("%*s+%0*d\n", prefix_width, "", bar_width, 0);
  printf("%*s0%*.*f%s\n", prefix_width, "", bar_width - 1, decimals,
         max_val, unit ? unit : "");
}

static double bench_elapsed_sec(struct timespec start, struct timespec end)
{
  return (double)(end.tv_sec - start.tv_sec) +
         (double)(end.tv_nsec - start.tv_nsec) / 1e9;
}

static void bench_free_groups(bench_group_t* time_groups,
                              bench_group_t* tps_groups,
                              char (*labels)[16], size_t count)
{
  if (time_groups) {
    for (size_t i = 0; i < count; i++) free((void*)time_groups[i].series);
  }
  if (tps_groups) {
    for (size_t i = 0; i < count; i++) free((void*)tps_groups[i].series);
  }
  free(time_groups);
  free(tps_groups);
  free(labels);
}

static void bench_run_and_plot(const char* suite_name,
                               const bench_algo_t* algos, size_t num_algos,
                               const bench_case_t* cases, size_t num_cases,
                               size_t iterations)
{
  if (!algos || !cases || num_algos == 0 || num_cases == 0 || iterations == 0)
    return;

  bench_group_t* time_groups = calloc(num_cases, sizeof(bench_group_t));
  bench_group_t* tps_groups = calloc(num_cases, sizeof(bench_group_t));
  char (*labels)[16] = calloc(num_cases, sizeof(char[16]));

  if (!time_groups || !tps_groups || !labels) {
    bench_free_groups(time_groups, tps_groups, labels, num_cases);
    return;
  }

  size_t built = 0;
  for (size_t c = 0; c < num_cases; c++) {
    snprintf(labels[c], sizeof(labels[c]), "%db", cases[c].bits);

    bench_series_t* t_series = calloc(num_algos, sizeof(bench_series_t));
    bench_series_t* p_series = calloc(num_algos, sizeof(bench_series_t));
    if (!t_series || !p_series) {
      free(t_series);
      free(p_series);
      break;
    }

    for (size_t a = 0; a < num_algos; a++) {
      struct timespec start, end;

      printf("Benchmarking %s (%s): %zu iters at %d bits...\n",
             suite_name ? suite_name : "Suite",
             algos[a].name ? algos[a].name : "algo",
             iterations, cases[c].bits);

      clock_gettime(CLOCK_MONOTONIC, &start);
      uint64_t acc = algos[a].fn(cases[c].ctx, iterations);
      clock_gettime(CLOCK_MONOTONIC, &end);

      bench_sink = bench_sink ^ acc;

      double secs = bench_elapsed_sec(start, end);
      double tps = (secs > 0.0) ? ((double)iterations / secs) : 0.0;

      printf("  Total Time: %.4fs | Throughput: %.2f tests/sec\n", secs, tps);

      t_series[a] = (bench_series_t){algos[a].name, algos[a].color, secs};
      p_series[a] = (bench_series_t){algos[a].name, algos[a].color, tps};
    }

    time_groups[built] = (bench_group_t){labels[built], t_series, num_algos};
    tps_groups[built] = (bench_group_t){labels[built], p_series, num_algos};
    built++;
  }

  if (built > 0) {
    char title_buf[128];
    const char* base = suite_name ? suite_name : "Benchmark";

    snprintf(title_buf, sizeof(title_buf), "%s — Total Time", base);
    bench_plot_metric(title_buf, BENCH_LOWER_IS_BETTER, time_groups, built,
                      "s", 4);

    snprintf(title_buf, sizeof(title_buf), "%s — Throughput", base);
    bench_plot_metric(title_buf, BENCH_HIGHER_IS_BETTER, tps_groups, built,
                      " tests/sec", 2);
  }

  bench_free_groups(time_groups, tps_groups, labels, built);
}

#endif /* BENCH_PLOT_H */
