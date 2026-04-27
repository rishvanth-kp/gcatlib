/*
 * test_genomic_step_vector.cpp: Catch2 correctness tests for GenomicStepVector
 */

#include "catch_amalgamated.hpp"

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
