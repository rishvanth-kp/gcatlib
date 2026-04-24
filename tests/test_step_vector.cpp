/*
 * test_step_vector.cpp: Catch2 correctness tests for StepVector
 */

#include "catch_amalgamated.hpp"

#include <utility>
#include <vector>

#include "StepVector.hpp"

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Convenience: build a pair with size_t first so comparisons type-check.
static std::pair<size_t, int> sp(size_t pos, int val) {
  return {pos, val};
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

TEST_CASE("StepVector: empty vector") {

  StepVector<int> sv;

  SECTION("at() returns 0 for any position") {
    REQUIRE(sv.at(0)   == 0);
    REQUIRE(sv.at(100) == 0);
  }

  SECTION("at_range() returns two zero-valued boundary points") {
    std::vector<std::pair<size_t, int>> out;
    sv.at_range(10, 20, out);
    REQUIRE(out.size() == 2);
    REQUIRE(out[0] == sp(10, 0));
    REQUIRE(out[1] == sp(20, 0));
  }
}

// ---------------------------------------------------------------------------

TEST_CASE("StepVector: single interval [10, 20) with value 5") {

  StepVector<int> sv;
  sv.add(10, 20, 5);

  SECTION("at() inside interval") {
    REQUIRE(sv.at(10) == 5);   // inclusive start
    REQUIRE(sv.at(15) == 5);   // middle
    REQUIRE(sv.at(19) == 5);   // last included position
  }

  SECTION("at() outside interval") {
    REQUIRE(sv.at(9)  == 0);   // just before
    REQUIRE(sv.at(20) == 0);   // exclusive end
    REQUIRE(sv.at(25) == 0);   // well after
  }

  SECTION("add() is a no-op when start >= end") {
    sv.add(15, 10, 99);        // start > end
    REQUIRE(sv.at(12) == 5);
    sv.add(15, 15, 99);        // start == end
    REQUIRE(sv.at(15) == 5);
  }
}

// ---------------------------------------------------------------------------

TEST_CASE("StepVector: overlapping intervals accumulate") {
  //  [10 -------- 20)  val=5
  //          [15 ---------- 25)  val=3
  //
  //  Result: [10,15)=5, [15,20)=8, [20,25)=3

  StepVector<int> sv;
  sv.add(10, 20, 5);
  sv.add(15, 25, 3);

  REQUIRE(sv.at(9)  == 0);   // before both
  REQUIRE(sv.at(12) == 5);   // first only
  REQUIRE(sv.at(17) == 8);   // overlap: 5+3
  REQUIRE(sv.at(22) == 3);   // second only
  REQUIRE(sv.at(25) == 0);   // after both
}

// ---------------------------------------------------------------------------

TEST_CASE("StepVector: adjacent non-overlapping intervals") {
  // [0, 10) val=1  |  [10, 20) val=2
  // No gap, no overlap — boundary at 10 must switch cleanly.

  StepVector<int> sv;
  sv.add(0,  10, 1);
  sv.add(10, 20, 2);

  REQUIRE(sv.at(5)  == 1);
  REQUIRE(sv.at(10) == 2);   // second interval starts exactly here
  REQUIRE(sv.at(15) == 2);
  REQUIRE(sv.at(20) == 0);
}

// ---------------------------------------------------------------------------

TEST_CASE("StepVector: nested intervals") {
  // [0 ------------------- 30)  val=1
  //       [10 ---- 20)  val=2
  //
  // Result: [0,10)=1, [10,20)=3, [20,30)=1

  StepVector<int> sv;
  sv.add(0,  30, 1);
  sv.add(10, 20, 2);

  REQUIRE(sv.at(5)  == 1);   // outer only
  REQUIRE(sv.at(15) == 3);   // outer + inner: 1+2
  REQUIRE(sv.at(25) == 1);   // outer only again
}

// ---------------------------------------------------------------------------

TEST_CASE("StepVector: at_range") {
  // Build: add(10,30,1) then add(20,40,2)
  // Internal steps: {10:1, 20:3, 30:2, 40:0}

  StepVector<int> sv;
  sv.add(10, 30, 1);
  sv.add(20, 40, 2);

  std::vector<std::pair<size_t, int>> out;

  SECTION("range spanning multiple step boundaries") {
    sv.at_range(15, 35, out);
    REQUIRE(out.size() == 4);
    REQUIRE(out[0] == sp(15, 1));
    REQUIRE(out[1] == sp(20, 3));
    REQUIRE(out[2] == sp(30, 2));
    REQUIRE(out[3] == sp(35, 2));
  }

  SECTION("range within a single step") {
    sv.at_range(11, 19, out);
    REQUIRE(out.size() == 2);
    REQUIRE(out[0] == sp(11, 1));
    REQUIRE(out[1] == sp(19, 1));
  }

  SECTION("range wider than all data") {
    sv.at_range(0, 50, out);
    REQUIRE(out.front().second == 0);   // zero before any data
    REQUIRE(out.back().second  == 0);   // zero after all data
  }

  SECTION("start == end returns empty output") {
    sv.at_range(15, 15, out);
    REQUIRE(out.empty());
  }

  SECTION("query boundaries coincide with step boundaries") {
    sv.at_range(10, 40, out);
    REQUIRE(out.front() == sp(10, 1));
    REQUIRE(out.back()  == sp(40, 0));
  }
}
