/*
 * test_genomic_step_vector.cpp: Catch2 correctness tests for GenomicStepVector
 */

#include "catch_amalgamated.hpp"

#include <algorithm>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "GenomicStepVector.hpp"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// GenomicRegion has no operator==, so compare fields explicitly.
static bool region_eq(const GenomicRegion &a, const GenomicRegion &b) {
  return a.name == b.name && a.start == b.start && a.end == b.end;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: empty vector") {

  GenomicStepVector<int> gsv;

  SECTION("counters start at zero") {
    REQUIRE(gsv.chrom_count() == 0);
    REQUIRE(gsv.entry_count() == 0);
  }

  SECTION("query on non-existent chromosome returns empty") {
    std::vector<std::pair<GenomicRegion, int>> out;
    gsv.at(GenomicRegion{"chr1", 0, 100}, out);
    REQUIRE(out.empty());
  }
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: chrom_count and entry_count") {

  GenomicStepVector<int> gsv;

  gsv.add("chr1", 0, 100, 1);
  REQUIRE(gsv.chrom_count() == 1);
  REQUIRE(gsv.entry_count() == 1);

  // Second add on the same chromosome: new entry, no new chromosome.
  gsv.add("chr1", 50, 150, 1);
  REQUIRE(gsv.chrom_count() == 1);
  REQUIRE(gsv.entry_count() == 2);

  // Add on a new chromosome: both counters increment.
  gsv.add("chr2", 0, 100, 1);
  REQUIRE(gsv.chrom_count() == 2);
  REQUIRE(gsv.entry_count() == 3);
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: single chromosome, single interval") {

  GenomicStepVector<int> gsv;
  gsv.add("chr1", 10, 20, 5);

  std::vector<std::pair<GenomicRegion, int>> out;
  gsv.at(GenomicRegion{"chr1", 0, 30}, out);

  REQUIRE(out.size() == 1);
  REQUIRE(region_eq(out[0].first, GenomicRegion{"chr1", 10, 20}));
  REQUIRE(out[0].second == 5);
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: chromosome isolation") {
  // Data added to chr1 must not appear in chr2 queries and vice versa.

  GenomicStepVector<int> gsv;
  gsv.add("chr1", 10, 20, 5);
  gsv.add("chr2", 30, 40, 3);

  std::vector<std::pair<GenomicRegion, int>> out;

  SECTION("chr1 query sees only chr1 data") {
    gsv.at(GenomicRegion{"chr1", 0, 50}, out);
    REQUIRE(out.size() == 1);
    REQUIRE(region_eq(out[0].first, GenomicRegion{"chr1", 10, 20}));
    REQUIRE(out[0].second == 5);
  }

  SECTION("chr2 query sees only chr2 data") {
    gsv.at(GenomicRegion{"chr2", 0, 50}, out);
    REQUIRE(out.size() == 1);
    REQUIRE(region_eq(out[0].first, GenomicRegion{"chr2", 30, 40}));
    REQUIRE(out[0].second == 3);
  }

  SECTION("query on a third chromosome returns empty") {
    gsv.at(GenomicRegion{"chr3", 0, 50}, out);
    REQUIRE(out.empty());
  }
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: overlapping intervals on same chromosome") {
  //  [10 -------- 20)  val=5
  //          [15 ---------- 25)  val=3
  //
  //  Output: [10,15)=5, [15,20)=8, [20,25)=3

  GenomicStepVector<int> gsv;
  gsv.add("chr1", 10, 20, 5);
  gsv.add("chr1", 15, 25, 3);

  std::vector<std::pair<GenomicRegion, int>> out;
  gsv.at(GenomicRegion{"chr1", 0, 30}, out);

  REQUIRE(out.size() == 3);
  REQUIRE(region_eq(out[0].first, GenomicRegion{"chr1", 10, 15}));
  REQUIRE(out[0].second == 5);
  REQUIRE(region_eq(out[1].first, GenomicRegion{"chr1", 15, 20}));
  REQUIRE(out[1].second == 8);
  REQUIRE(region_eq(out[2].first, GenomicRegion{"chr1", 20, 25}));
  REQUIRE(out[2].second == 3);
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: GenomicRegion overload of add()") {
  // add(GenomicRegion, val) must behave identically to add(chr, start, end, val).

  GenomicStepVector<int> gsv;
  gsv.add(GenomicRegion{"chr1", 10, 20}, 5);

  std::vector<std::pair<GenomicRegion, int>> out;
  gsv.at(GenomicRegion{"chr1", 0, 30}, out);

  REQUIRE(out.size() == 1);
  REQUIRE(region_eq(out[0].first, GenomicRegion{"chr1", 10, 20}));
  REQUIRE(out[0].second == 5);
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: keep_0=false excludes zero-valued regions") {

  GenomicStepVector<int> gsv;
  gsv.add("chr1", 10, 20, 5);

  std::vector<std::pair<GenomicRegion, int>> out;
  // Default keep_0=false: regions [0,10) and [20,30) have value 0 and are excluded.
  gsv.at(GenomicRegion{"chr1", 0, 30}, out);

  REQUIRE(out.size() == 1);
  REQUIRE(region_eq(out[0].first, GenomicRegion{"chr1", 10, 20}));
  REQUIRE(out[0].second == 5);
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: keep_0=true includes zero-valued regions") {

  GenomicStepVector<int> gsv;
  gsv.add("chr1", 10, 20, 5);

  std::vector<std::pair<GenomicRegion, int>> out;
  gsv.at(GenomicRegion{"chr1", 0, 30}, out, true);

  // Expect three regions: [0,10)=0, [10,20)=5, [20,30)=0
  REQUIRE(out.size() == 3);
  REQUIRE(region_eq(out[0].first, GenomicRegion{"chr1",  0, 10}));
  REQUIRE(out[0].second == 0);
  REQUIRE(region_eq(out[1].first, GenomicRegion{"chr1", 10, 20}));
  REQUIRE(out[1].second == 5);
  REQUIRE(region_eq(out[2].first, GenomicRegion{"chr1", 20, 30}));
  REQUIRE(out[2].second == 0);
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: equal-value steps are merged in output") {
  // Two adjacent intervals with the same value should collapse into one
  // output GenomicRegion rather than two.

  GenomicStepVector<int> gsv;
  gsv.add("chr1", 10, 20, 3);
  gsv.add("chr1", 20, 30, 3);   // adjacent, same value

  std::vector<std::pair<GenomicRegion, int>> out;
  gsv.at(GenomicRegion{"chr1", 0, 40}, out);

  REQUIRE(out.size() == 1);
  REQUIRE(region_eq(out[0].first, GenomicRegion{"chr1", 10, 30}));
  REQUIRE(out[0].second == 3);
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: single-position intervals") {
  // Minimal-width interval (end = start + 1).

  GenomicStepVector<int> gsv;
  gsv.add("chr1", 10, 11, 5);

  std::vector<std::pair<GenomicRegion, int>> out;
  gsv.at(GenomicRegion{"chr1", 0, 20}, out);

  REQUIRE(out.size() == 1);
  REQUIRE(region_eq(out[0].first, GenomicRegion{"chr1", 10, 11}));
  REQUIRE(out[0].second == 5);
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: order independence") {
  // Adding intervals in two different orders must produce identical results.

  GenomicStepVector<int> gsv1, gsv2;

  gsv1.add("chr1", 10, 30, 1);
  gsv1.add("chr1", 20, 40, 2);
  gsv1.add("chr1", 15, 25, 3);

  gsv2.add("chr1", 15, 25, 3);   // reversed order
  gsv2.add("chr1", 20, 40, 2);
  gsv2.add("chr1", 10, 30, 1);

  std::vector<std::pair<GenomicRegion, int>> out1, out2;
  gsv1.at(GenomicRegion{"chr1", 0, 50}, out1);
  gsv2.at(GenomicRegion{"chr1", 0, 50}, out2);

  REQUIRE(out1.size() == out2.size());
  for (size_t i = 0; i < out1.size(); ++i) {
    REQUIRE(region_eq(out1[i].first, out2[i].first));
    REQUIRE(out1[i].second == out2[i].second);
  }
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: accumulated value matches interval count") {
  // Add N random intervals with value=1. For a set of query positions,
  // verify that at(pos) equals the brute-force count of intervals containing pos.

  const int    N     = 100;
  const size_t range = 1000;

  std::mt19937 rng(123);
  std::uniform_int_distribution<size_t> pos_dist(0, range - 1);
  std::uniform_int_distribution<size_t> len_dist(1, 100);

  struct Interval { size_t start, end; };
  std::vector<Interval> intervals(N);

  GenomicStepVector<int> gsv;
  for (auto &iv : intervals) {
    size_t start = pos_dist(rng);
    size_t end   = std::min(start + len_dist(rng), range);
    iv = {start, end};
    gsv.add("chr1", start, end, 1);
  }

  const std::vector<size_t> test_positions = {0, 50, 100, 250, 500, 750, 900, 999};
  std::vector<std::pair<GenomicRegion, int>> out;

  for (size_t pos : test_positions) {
    int expected = 0;
    for (const auto &iv : intervals)
      if (iv.start <= pos && pos < iv.end)
        ++expected;

    // Query [pos, pos+1) with keep_0=true to always get exactly one entry.
    gsv.at(GenomicRegion{"chr1", pos, pos + 1}, out, true);
    REQUIRE(out.size() == 1);
    REQUIRE(out[0].second == expected);
  }
}

// ---------------------------------------------------------------------------

TEST_CASE("GenomicStepVector: cross-chromosome isolation at scale") {
  // Add the same two overlapping intervals to 20 chromosomes.
  // Each chromosome must independently return the correct merged output.

  const int n_chrom = 20;
  GenomicStepVector<int> gsv;

  for (int c = 1; c <= n_chrom; ++c) {
    std::string chr = "chr" + std::to_string(c);
    gsv.add(chr, 10, 30, 1);
    gsv.add(chr, 20, 40, 2);
  }

  std::vector<std::pair<GenomicRegion, int>> out;
  for (int c = 1; c <= n_chrom; ++c) {
    std::string chr = "chr" + std::to_string(c);
    gsv.at(GenomicRegion{chr, 0, 50}, out);

    REQUIRE(out.size() == 3);
    REQUIRE(region_eq(out[0].first, GenomicRegion{chr, 10, 20}));
    REQUIRE(out[0].second == 1);
    REQUIRE(region_eq(out[1].first, GenomicRegion{chr, 20, 30}));
    REQUIRE(out[1].second == 3);
    REQUIRE(region_eq(out[2].first, GenomicRegion{chr, 30, 40}));
    REQUIRE(out[2].second == 2);
  }
}
