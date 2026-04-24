/*
 * bench_step_vector.cpp: nanobench performance benchmarks for StepVector
 *
 * Usage: ./bench_step_vector
 *
 * Sections:
 *   1. add() — sequential non-overlapping intervals
 *   2. add() — random overlapping intervals (pre-generated)
 *   3. at()  — point query after N adds
 *   4. at_range() — full-range query after N adds
 */

#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

#include <random>
#include <string>
#include <utility>
#include <vector>

#include "StepVector.hpp"

namespace nb = ankerl::nanobench;

int main() {

  // -------------------------------------------------------------------------
  // 1. add() — sequential non-overlapping intervals
  //
  // Best-case pattern: each new interval starts after the previous one ends,
  // so there are no overlaps and the internal map grows monotonically.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("StepVector::add  sequential non-overlapping")
         .unit("add")
         .warmup(3);

    for (int n : {1000, 10000, 100000, 1000000}) {
      bench.run("N=" + std::to_string(n), [&] {
        StepVector<int> sv;
        for (int i = 0; i < n; ++i)
          sv.add(static_cast<size_t>(i) * 10,
                 static_cast<size_t>(i) * 10 + 5, 1);
        nb::doNotOptimizeAway(sv);
      });
    }
  }

  // -------------------------------------------------------------------------
  // 2. add() — random overlapping intervals
  //
  // Realistic pattern: random start positions and lengths over a fixed genome-
  // like range. Intervals overlap frequently, exercising the middle-of-range
  // update loop. Intervals are pre-generated to avoid measuring RNG overhead.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("StepVector::add  random overlapping")
         .unit("add")
         .warmup(3);

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 100000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n : {1000, 10000, 100000}) {
      std::vector<std::pair<size_t, size_t>> intervals(n);
      for (auto &iv : intervals) {
        size_t start = pos_dist(rng);
        iv = {start, start + len_dist(rng)};
      }

      bench.run("N=" + std::to_string(n), [&] {
        StepVector<int> sv;
        for (const auto &iv : intervals)
          sv.add(iv.first, iv.second, 1);
        nb::doNotOptimizeAway(sv);
      });
    }
  }

  // -------------------------------------------------------------------------
  // 3. at() — point query after N adds
  //
  // Measures upper_bound() lookup cost as the map grows. Query positions are
  // pre-generated and spread uniformly across the populated range.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("StepVector::at  point query")
         .unit("query")
         .warmup(3);

    std::mt19937 rng(42);

    for (int n : {1000, 10000, 100000}) {
      StepVector<int> sv;
      for (int i = 0; i < n; ++i)
        sv.add(static_cast<size_t>(i) * 10,
               static_cast<size_t>(i) * 10 + 5, 1);

      std::uniform_int_distribution<size_t> qpos(0, static_cast<size_t>(n) * 10);
      std::vector<size_t> queries(1000);
      for (auto &q : queries) q = qpos(rng);

      bench.run("N=" + std::to_string(n), [&] {
        for (size_t q : queries)
          nb::doNotOptimizeAway(sv.at(q));
      });
    }
  }

  // -------------------------------------------------------------------------
  // 4. at_range() — full-range query after N adds
  //
  // Measures the cost of collecting all step boundaries into an output vector.
  // The output vector is reused across iterations so allocation cost is not
  // included after the first run, reflecting real-world buffer-reuse patterns.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("StepVector::at_range  full range query")
         .unit("query")
         .warmup(3);

    for (int n : {1000, 10000, 100000}) {
      StepVector<int> sv;
      for (int i = 0; i < n; ++i)
        sv.add(static_cast<size_t>(i) * 10,
               static_cast<size_t>(i) * 10 + 5, 1);

      std::vector<std::pair<size_t, int>> out;

      bench.run("N=" + std::to_string(n), [&] {
        sv.at_range(0, static_cast<size_t>(n) * 10, out);
        nb::doNotOptimizeAway(out);
      });
    }
  }

  return 0;
}
