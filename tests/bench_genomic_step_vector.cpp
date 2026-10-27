/*
 * bench_genomic_step_vector.cpp: nanobench benchmarks for GenomicStepVector
 *
 * Usage: ./bench_genomic_step_vector
 *
 * Results are printed to stdout and written as JSON to
 * genomic_step_vector_bench.json.
 *
 * Sections:
 *   1. add() — single chromosome, sequential non-overlapping
 *   2. add() — single chromosome, random overlapping
 *   3. add() — single chromosome, random overlapping, sorted
 *   4. add() — 1, 5, 10, 20 chromosomes, random overlapping, sorted,
 *              fixed total number of adds
 *   5. add() — 1, 5, 10, 20 chromosomes, random overlapping, sorted,
 *              fixed number of adds per chromosome
 *   6. at()  — single chromosome, unsorted queries, varying query count
 *   7. at()  — single chromosome, sorted queries, varying query count
 *   8. at()  — single chromosome, sorted queries, varying entry count
 *
 * Disabled (commented out): multi-chromosome sequential add(),
 * multi-chromosome at(), and the two STRESS add() benchmarks.
 */

#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

#include <random>
#include <string>
#include <utility>
#include <vector>

#include "GenomicStepVector.hpp"

using std::cout;
using std::cerr;
using std::endl;

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

  // log file
  std::ofstream log_file("genomic_step_vector_bench.json");

  // -------------------------------------------------------------------------
  // add() — single chromosome, sequential non-overlapping
  //
  // Baseline: measures the cost of the unordered_map chromosome dispatch on
  // top of StepVector::add(). Compare directly with bench_step_vector results.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  single chromosome sequential")
         .unit("add")
         .epochs(10)
         .warmup(3);

    for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
      bench.complexityN(n).run("N=" + std::to_string(n), [&] {
        GenomicStepVector<int> gsv;
        for (int i = 0; i < n; ++i)
          gsv.add("chr1",
                  static_cast<size_t>(i) * 10,
                  static_cast<size_t>(i) * 10 + 5, 1);
        nb::doNotOptimizeAway(gsv);
      });
    }

    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);
  }

/*
  // -------------------------------------------------------------------------
  // add() — 1, 5, 20 chromosomes, fixed adds per chromosome
  //
  // Each chromosome receives the same number of sequential non-overlapping
  // adds. Measures how the unordered_map scales as chromosome count grows
  // while per-chromosome work stays constant.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  multi-chromosome sequential")
         .unit("add")
         .epochs(10)
         .warmup(3);

    for (int n_chrom : {1, 5, 20}) {
      const std::vector<std::string> chroms = make_chrom_names(n_chrom);

      for (int n_per_chrom : {1000, 10000, 100000, 1000000, 10000000}) {
        const std::string label = "n_chrom=" + std::to_string(n_chrom) +
                                  " n_per_chrom=" + std::to_string(n_per_chrom);
        bench.complexityN(n_chrom * n_per_chrom).run(label, [&] {
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

    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);
  }
*/

  // -------------------------------------------------------------------------
  // add() — 1 chromosome, random overlapping (pre-generated)
  //
  // Each chromosome receives the same pre-generated random intervals.
  // Measures multi-chromosome performance under realistic overlapping data.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  single chromosome random overlapping")
         .unit("add")
         .epochs(10)
         .warmup(3);

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 10000000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
      // Pre-generate one set of intervals shared across chromosomes.
      std::vector<std::pair<size_t, size_t>> intervals(n);
      for (auto &iv : intervals) {
        size_t start = pos_dist(rng);
        iv = {start, start + len_dist(rng)};
      }

      const std::string label = "N=" + std::to_string(n);
      bench.complexityN(n).run(label, [&] {
        GenomicStepVector<int> gsv;
        for (const auto &iv : intervals)
          gsv.add("chr1", iv.first, iv.second, 1);
        nb::doNotOptimizeAway(gsv);
      });
    }
    
    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);
  }
  
  // -------------------------------------------------------------------------
  // add() — 1 chromosome, random, sorted, and overlapping (pre-generated)
  //
  // Each chromosome receives the same pre-generated random intervals.
  // Measures multi-chromosome performance under realistic overlapping data.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  single chromosome random sorted overlapping")
         .unit("add")
         .epochs(10)
         .warmup(3);

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 10000000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
      // Pre-generate one set of intervals shared across chromosomes.
      std::vector<std::pair<size_t, size_t>> intervals(n);
      for (auto &iv : intervals) {
        size_t start = pos_dist(rng);
        iv = {start, start + len_dist(rng)};
      }
        
      std::sort(intervals.begin(), intervals.end());

      const std::string label = "N=" + std::to_string(n);
      bench.complexityN(n).run(label, [&] {
        GenomicStepVector<int> gsv;
        for (const auto &iv : intervals)
          gsv.add("chr1", iv.first, iv.second, 1);
        nb::doNotOptimizeAway(gsv);
      });
    }

    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);  
  }
  
  // -------------------------------------------------------------------------
  // add() — 1, 5, 10, 20 chromosomes, random, sorted, and overlapping 
  // (pre-generated)
  //
  // Measures multi-chromosome performance under realistic overlapping data.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  multi-chromosome random sorted overlapping")
         .unit("add")
         .epochs(10)
         .warmup(3);

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 10000000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n_chrom : {1, 5, 10, 20}) {
      const std::vector<std::string> chroms = make_chrom_names(n_chrom);  
      
      for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
        // pre-generate intervels and assign to a random chrom
        std::uniform_int_distribution<size_t> chrom_dist(0, n_chrom - 1);
        std::vector<GenomicRegion> intervals;
        for (int i = 0; i < n; ++i) {
          size_t start = pos_dist(rng);
          size_t end = start + len_dist(rng);
          intervals.push_back(GenomicRegion{chroms[chrom_dist(rng)], start, end});
        }
 
        // sort the GenomicRanges
        std::sort(intervals.begin(), intervals.end(), 
                  [] (const GenomicRegion& a, const GenomicRegion &b) {
                    return std::tie(a.name, a.start) < std::tie(b.name, b.start);
                  });
        

        const std::string label = "n_chrom=" + std::to_string(n_chrom) +
                                  " n=" + std::to_string(n);
        bench.complexityN(n).run(label, [&] {
          GenomicStepVector<int> gsv;
          for (const auto &iv : intervals)
            gsv.add(iv, 1);
          nb::doNotOptimizeAway(gsv);
        });
      }
    }

    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);  
  } 
  
  // -------------------------------------------------------------------------
  // add() — 1, 5, 10, 20 chromosomes, random, sorted, and overlapping 
  // (pre-generated)
  //
  // Measures multi-chromosome performance under realistic overlapping data.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  multi-chromosome random sorted overlapping 2")
         .unit("add")
         .epochs(10)
         .warmup(3);

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 10000000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n_chrom : {1, 5, 10, 20}) {
      const std::vector<std::string> chroms = make_chrom_names(n_chrom);  
      
      for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
        // pre-generate intervels and assign to a random chrom
        std::uniform_int_distribution<size_t> chrom_dist(0, n_chrom - 1);
        std::vector<GenomicRegion> intervals;
        for (int i = 0; i < n * n_chrom; ++i) {
          size_t start = pos_dist(rng);
          size_t end = start + len_dist(rng);
          intervals.push_back(GenomicRegion{chroms[chrom_dist(rng)], start, end});
        }
 
        // sort the GenomicRanges
        std::sort(intervals.begin(), intervals.end(), 
                  [] (const GenomicRegion& a, const GenomicRegion &b) {
                    return std::tie(a.name, a.start) < std::tie(b.name, b.start);
                  });
        

        const std::string label = "n_chrom=" + std::to_string(n_chrom) +
                                  " n_per_chrom=" + std::to_string(n);
        bench.complexityN(n * n_chrom).run(label, [&] {
          GenomicStepVector<int> gsv;
          for (const auto &iv : intervals)
            gsv.add(iv, 1);
          nb::doNotOptimizeAway(gsv);
        });
      }
    }

    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);  
  } 

/* 
  // -------------------------------------------------------------------------
  // add() — 1, 5, 20 chromosomes, random, sorted, and overlapping (pre-generated)
  //
  // Each chromosome receives the same pre-generated random intervals.
  // Measures multi-chromosome performance under realistic overlapping data.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  multi-chromosome random sorted overlapping")
         .unit("add")
         .epochs(10)
         .warmup(3);

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 10000000);
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

        std::sort(intervals.begin(), intervals.end());

        const std::string label = "n_chrom=" + std::to_string(n_chrom) +
                                  " n_per_chrom=" + std::to_string(n_per_chrom);
        bench.complexityN(n_chrom * n_per_chrom).run(label, [&] {
          GenomicStepVector<int> gsv;
          for (const auto &chr : chroms)
            for (const auto &iv : intervals)
              gsv.add(chr, iv.first, iv.second, 1);
          nb::doNotOptimizeAway(gsv);
        });
      }
    }
    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);
  }
*/



  // -------------------------------------------------------------------------
  // at() — single chromosome unsorted query with fixed number of entires. 
  //
  // Builds a populated GenomicStepVector, then benchmarks at() queries with
  // pre-generated unsorted (chromosome, region) pairs. Measures query cost
  // as the number of queries grow.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::at single-chromosome unsorted query")
         .unit("query")
         .epochs(10)
         .warmup(3);

    const int n_entries = 1000000;

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 10000000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
      // build the genomic step vector
      std::vector<std::pair<size_t, size_t>> entries(n_entries);
      for (auto &e : entries) {
        size_t start = pos_dist(rng);
        e = {start, start + len_dist(rng)};
      }
      
      std::sort(entries.begin(), entries.end());

      GenomicStepVector<int> gsv;
      for(auto &e : entries) {
        gsv.add("chr1", e.first, e.second, 1);
      }

      // pre-generate random queries
      std::vector<GenomicRegion> queries(n); 
      for (auto &q : queries) {
        size_t start = pos_dist(rng);
        q = GenomicRegion{"chr1", start, start + len_dist(rng)};
      }

      // query the step vector
      std::vector<std::pair<GenomicRegion, int>> out;
      const std::string label = "N=" + std::to_string(n);
      bench.complexityN(n).run(label, [&] {
        for (const auto &q : queries) {
          gsv.at(q, out);
          nb::doNotOptimizeAway(out);
        }
      }); 
    }    
    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);
  }
  
  // -------------------------------------------------------------------------
  // at() — single chromosome sorted query with fixed number of entires. 
  //
  // Builds a populated GenomicStepVector, then benchmarks at() queries with
  // pre-generated unsorted (chromosome, region) pairs. Measures query cost
  // as the number of queries grow.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::at single-chromosome sorted query")
         .unit("query")
         .epochs(10)
         .warmup(3);

    const int n_entries = 1000000;

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 10000000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
      // build the genomic step vector
      std::vector<std::pair<size_t, size_t>> entries(n_entries);
      for (auto &e : entries) {
        size_t start = pos_dist(rng);
        e = {start, start + len_dist(rng)};
      }
      
      std::sort(entries.begin(), entries.end());

      GenomicStepVector<int> gsv;
      for(auto &e : entries) {
        gsv.add("chr1", e.first, e.second, 1);
      }

      // pre-generate random queries
      std::vector<GenomicRegion> queries(n); 
      for (auto &q : queries) {
        size_t start = pos_dist(rng);
        q = GenomicRegion{"chr1", start, start + len_dist(rng)};
      }

      std::sort(queries.begin(), queries.end(), 
                [] (const GenomicRegion& a, const GenomicRegion &b) {
                  return std::tie(a.name, a.start) < std::tie(b.name, b.start);
                });

      // query the step vector
      std::vector<std::pair<GenomicRegion, int>> out;
      const std::string label = "N=" + std::to_string(n);
      bench.complexityN(n).run(label, [&] {
        for (const auto &q : queries) {
          gsv.at(q, out);
          nb::doNotOptimizeAway(out);
        }
      }); 
    }    
    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);
  }

  // -------------------------------------------------------------------------
  // at() — single chromosome varying the number of entries with fixed number 
  // of queries. 
  //
  // Builds a populated GenomicStepVector, then benchmarks at() queries with
  // pre-generated unsorted (chromosome, region) pairs. Measures query cost
  // as the number of entries grow.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::at single-chromosome varying entries")
         .unit("query")
         .epochs(10)
         .warmup(3);

    const int n_queries = 1000000;

    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> pos_dist(0, 10000000);
    std::uniform_int_distribution<size_t> len_dist(1, 1000);

    for (int n : {1000, 10000, 100000, 1000000, 10000000}) {
      // build the genomic step vector
      std::vector<std::pair<size_t, size_t>> entries(n);
      for (auto &e : entries) {
        size_t start = pos_dist(rng);
        e = {start, start + len_dist(rng)};
      }
      
      std::sort(entries.begin(), entries.end());

      GenomicStepVector<int> gsv;
      for(auto &e : entries) {
        gsv.add("chr1", e.first, e.second, 1);
      }

      // pre-generate random queries
      std::vector<GenomicRegion> queries(n_queries); 
      for (auto &q : queries) {
        size_t start = pos_dist(rng);
        q = GenomicRegion{"chr1", start, start + len_dist(rng)};
      }

      std::sort(queries.begin(), queries.end(), 
                [] (const GenomicRegion& a, const GenomicRegion &b) {
                  return std::tie(a.name, a.start) < std::tie(b.name, b.start);
                });

      // query the step vector
      std::vector<std::pair<GenomicRegion, int>> out;
      const std::string label = "N=" + std::to_string(n);
      bench.complexityN(n).run(label, [&] {
        for (const auto &q : queries) {
          gsv.at(q, out);
          nb::doNotOptimizeAway(out);
        }
      }); 
    }    
    cout << bench.complexityBigO() << endl;
    bench.render(ankerl::nanobench::templates::json(), log_file);
  }

/*
  // -------------------------------------------------------------------------
  // at() — query across 1, 5, 20 chromosomes
  //
  // Builds a populated GenomicStepVector, then benchmarks at() queries with
  // pre-generated random (chromosome, region) pairs. Measures chromosome
  // lookup cost and output merging as chromosome count grows.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::at  multi-chromosome query old")
         .unit("query")
         .epochs(10)
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
*/


/*
  // -------------------------------------------------------------------------
  // STRESS: worst-case add() — wide interval over dense pre-built map
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
         .epochs(10)
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
  // STRESS: maximally dense add() — sequential length-1 intervals
  //
  // Each add creates a new map entry (no merging possible), so the map
  // grows to its maximum possible size. Measures performance under peak
  // memory and cache pressure from a large, dense std::map.
  // -------------------------------------------------------------------------
  {
    nb::Bench bench;
    bench.title("GenomicStepVector::add  STRESS maximally dense length-1")
         .unit("add")
         .epochs(10)
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
*/
  // close the log file
  log_file.close();

  return 0;
}
