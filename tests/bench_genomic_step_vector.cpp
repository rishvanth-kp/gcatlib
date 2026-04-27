/*
 * bench_genomic_step_vector.cpp: nanobench benchmarks for GenomicStepVector
 *
 * Usage: ./bench_genomic_step_vector
 *
 * Sections:
 *   1. add() — single chromosome (unordered_map dispatch overhead baseline)
 *   2. add() — 1, 5, 20 chromosomes, fixed adds per chromosome
 *   3. add() — 1, 5, 20 chromosomes, random overlapping (extended N)
 *   4. at()  — query across 1, 5, 20 chromosomes
 *   5. STRESS: worst-case add() — wide interval over dense pre-built map
 *   6. STRESS: maximally dense add() — length-1 intervals
 */

#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

#include <random>
#include <string>
#include <utility>
#include <vector>

#include "GenomicStepVector.hpp"

namespace nb = ankerl::nanobench;

// Build a vector of chromosome name strings "chr1" .. "chrN".
static std::vector<std::string> make_chrom_names(int n) {
  std::vector<std::string> names;
  names.reserve(n);
  for (int i = 1; i <= n; ++i)
    names.push_back("chr" + std::to_string(i));
  return names;
}

int main() {

  // -------------------------------------------------------------------------
  // 1. add() — single chromosome, sequential non-overlapping
  //
  // Baseline: measures the cost of the unordered_map chromosome dispatch on
  // top of StepVector::add(). Compare directly with bench_step_vector results.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  single chromosome sequential")
         .unit("add")
         .warmup(3);

    for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
      bench.run("N=" + std::to_string(n), [&] {
        GenomicStepVector<int> gsv;
        for (int i = 0; i < n; ++i)
          gsv.add("chr1",
                  static_cast<size_t>(i) * 10,
                  static_cast<size_t>(i) * 10 + 5, 1);
        nb::doNotOptimizeAway(gsv);
      });
    }
  }

  // -------------------------------------------------------------------------
  // 2. add() — 1, 5, 20 chromosomes, fixed adds per chromosome
  //
  // Each chromosome receives the same number of sequential non-overlapping
  // adds. Measures how the unordered_map scales as chromosome count grows
  // while per-chromosome work stays constant.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  multi-chromosome sequential")
         .unit("add")
         .warmup(3);

    for (int n_chrom : {1, 5, 20}) {
      const std::vector<std::string> chroms = make_chrom_names(n_chrom);

      for (int n_per_chrom : {1000, 10000, 100000}) {
        const std::string label = "n_chrom=" + std::to_string(n_chrom) +
                                  " n_per_chrom=" + std::to_string(n_per_chrom);
        bench.run(label, [&] {
          GenomicStepVector<int> gsv;
          for (const auto &chr : chroms)
            for (int i = 0; i < n_per_chrom; ++i)
              gsv.add(chr,
                      static_cast<size_t>(i) * 10,
                      static_cast<size_t>(i) * 10 + 5, 1);
          nb::doNotOptimizeAway(gsv);
        });
      }
    }
  }

  // -------------------------------------------------------------------------
  // 3. add() — 1, 5, 20 chromosomes, random overlapping (pre-generated)
  //
  // Each chromosome receives the same pre-generated random intervals.
  // Measures multi-chromosome performance under realistic overlapping data.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  multi-chromosome random overlapping")
         .unit("add")
         .warmup(3);

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 100000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n_chrom : {1, 5, 20}) {
      const std::vector<std::string> chroms = make_chrom_names(n_chrom);

      for (int n_per_chrom : {1000, 10000, 100000, 1000000}) {
        // Pre-generate one set of intervals shared across chromosomes.
        std::vector<std::pair<size_t, size_t>> intervals(n_per_chrom);
        for (auto &iv : intervals) {
          size_t start = pos_dist(rng);
          iv = {start, start + len_dist(rng)};
        }

        const std::string label = "n_chrom=" + std::to_string(n_chrom) +
                                  " n_per_chrom=" + std::to_string(n_per_chrom);
        bench.run(label, [&] {
          GenomicStepVector<int> gsv;
          for (const auto &chr : chroms)
            for (const auto &iv : intervals)
              gsv.add(chr, iv.first, iv.second, 1);
          nb::doNotOptimizeAway(gsv);
        });
      }
    }
  }

  // -------------------------------------------------------------------------
  // 4. at() — query across 1, 5, 20 chromosomes
  //
  // Builds a populated GenomicStepVector, then benchmarks at() queries with
  // pre-generated random (chromosome, region) pairs. Measures chromosome
  // lookup cost and output merging as chromosome count grows.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::at  multi-chromosome query")
         .unit("query")
         .warmup(3);

    std::mt19937 rng(42);
    const int n_per_chrom = 10000;
    const int n_queries   = 1000;

    for (int n_chrom : {1, 5, 20}) {
      const std::vector<std::string> chroms = make_chrom_names(n_chrom);

      // Build the vector.
      GenomicStepVector<int> gsv;
      for (const auto &chr : chroms)
        for (int i = 0; i < n_per_chrom; ++i)
          gsv.add(chr,
                  static_cast<size_t>(i) * 10,
                  static_cast<size_t>(i) * 10 + 5, 1);

      // Pre-generate random queries: pick a random chromosome and a random
      // region within its populated range.
      std::uniform_int_distribution<int>    chrom_dist(0, n_chrom - 1);
      std::uniform_int_distribution<size_t> pos_dist(0,
                                  static_cast<size_t>(n_per_chrom) * 10);
      std::vector<GenomicRegion> queries;
      queries.reserve(n_queries);
      for (int q = 0; q < n_queries; ++q) {
        size_t start = pos_dist(rng);
        size_t end   = start + 10000;
        queries.push_back(GenomicRegion{chroms[chrom_dist(rng)], start, end});
      }

      std::vector<std::pair<GenomicRegion, int>> out;
      const std::string label = "n_chrom=" + std::to_string(n_chrom);
      bench.run(label, [&] {
        for (const auto &q : queries) {
          gsv.at(q, out);
          nb::doNotOptimizeAway(out);
        }
      });
    }
  }

  // -------------------------------------------------------------------------
  // 5. STRESS: worst-case add() — wide interval over dense pre-built map
  //
  // Pre-build a map with N non-overlapping length-1 intervals, creating N
  // step boundaries. Then benchmark a single add() call spanning the entire
  // range, which must traverse and update every existing boundary node.
  // Each benchmark iteration is O(N), clearly showing the update-loop cost.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  STRESS worst-case wide interval")
         .unit("add")
         .warmup(1);

    for (int n : {10000, 100000, 1000000}) {
      // Pre-build: N length-1 non-overlapping intervals create N map entries.
      GenomicStepVector<int> gsv;
      for (int i = 0; i < n; ++i)
        gsv.add("chr1",
                static_cast<size_t>(i) * 2,
                static_cast<size_t>(i) * 2 + 1, 1);

      const size_t wide_end = static_cast<size_t>(n) * 2;
      bench.run("N=" + std::to_string(n), [&] {
        // Each call traverses all N existing map nodes.
        gsv.add("chr1", 0, wide_end, 1);
        nb::doNotOptimizeAway(gsv);
      });
    }
  }

  // -------------------------------------------------------------------------
  // 6. STRESS: maximally dense add() — sequential length-1 intervals
  //
  // Each add creates a new map entry (no merging possible), so the map
  // grows to its maximum possible size. Measures performance under peak
  // memory and cache pressure from a large, dense std::map.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  STRESS maximally dense length-1")
         .unit("add")
         .warmup(1);

    for (int n : {100000, 1000000}) {
      bench.run("N=" + std::to_string(n), [&] {
        GenomicStepVector<int> gsv;
        for (int i = 0; i < n; ++i)
          gsv.add("chr1",
                  static_cast<size_t>(i) * 2,
                  static_cast<size_t>(i) * 2 + 1, 1);
        nb::doNotOptimizeAway(gsv);
      });
    }
  }

  return 0;
}
