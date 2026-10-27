# Performance

Performance of `GenomicStepVector`, measured with
[nanobench](https://nanobench.ankerl.com) v4.3.11. The benchmark program is
`tests/bench_genomic_step_vector.cpp`.

## Environment

| Item | Details |
|---|---|
| CPU | AMD Ryzen 9 5900X |
| Memory | 96 GB |
| Compiler flags | `-O2 -std=c++17` |
| System tuning | `pyperf system tune` (performance mode) |
| Core pinning | One core, with `taskset` |

The benchmarks were run after a reboot, with only necessary user programs
running.

## Method

- Each benchmark runs 3 warm-up runs and 10 measured runs (epochs). The
  times reported are nanobench's median, and nanobench's error estimate was
  below 5% for every measurement.
- All intervals and queries are generated before timing, with a fixed random
  seed, so random number generation is not measured.
- Unless stated otherwise, random intervals have starts uniformly distributed
  over a 10 Mbp range and lengths uniformly distributed in [1, 1,000] bp.
- Each measurement times one complete run, for example building a
  `GenomicStepVector` from N intervals. The per-operation times below are
  that total divided by the number of operations.

---

## add()

### Sequential, non-overlapping

N intervals of 5 bp, each starting 5 bp after the previous one ends, on one
chromosome. No interval overlaps an existing step, so this is the fastest
case.

| N | Total time | Time per add |
|---:|---:|---:|
| 1,000 | 130 µs | 130 ns |
| 10,000 | 2.39 ms | 239 ns |
| 100,000 | 33.2 ms | 332 ns |
| 1,000,000 | 483 ms | 483 ns |
| 10,000,000 | 6.74 s | 674 ns |

Each `add()` does a few `std::map` operations, each O(log N), so the total is
**O(N log N)** (nanobench's best fit, with 3% error).

### Random overlapping: unsorted vs sorted

N random intervals on one chromosome, added either in random order or sorted
by start position. Both columns use the same intervals.

| N | Unsorted: total | Unsorted: per add | Sorted: total | Sorted: per add | Sorted speed-up |
|---:|---:|---:|---:|---:|---:|
| 1,000 | 164 µs | 164 ns | 128 µs | 128 ns | 1.3× |
| 10,000 | 2.88 ms | 288 ns | 1.67 ms | 167 ns | 1.7× |
| 100,000 | 48.8 ms | 488 ns | 23.8 ms | 238 ns | 2.0× |
| 1,000,000 | 4.99 s | 4.99 µs | 509 ms | 509 ns | 9.8× |
| 10,000,000 | 335 s | 33.5 µs | 11.7 s | 1.17 µs | 28.6× |

An `add()` updates every step boundary inside the new interval. As more
intervals are added to the same 10 Mbp range, the boundaries get denser and
each `add()` updates more of them, so the cost per add grows with N. In
unsorted order the cost grows as **O(N²)** (nanobench's best fit, with 1.1%
error).

Sorting the input makes this much cheaper: **28.6× faster at N = 10M**. The
hardware counters suggest that the unsorted runs are limited by memory
access: their instructions per cycle (IPC) fall from 1.8 at N = 1,000 to
0.05 at N = 10M, while the sorted runs stay above 1.1. A likely reason is
that in sorted order, the boundaries an `add()` visits were created recently
and are still in cache, while in random order they are scattered across the
`std::map`.

### Multiple chromosomes

Sorted random intervals, each assigned to one of 1, 5, 10 or 20
chromosomes at random. Each chromosome covers the same 10 Mbp range.

**N intervals per chromosome** (on average). The total number of adds is
N × (number of chromosomes).

| N per chromosome | 1 chromosome | 5 chromosomes | 10 chromosomes | 20 chromosomes |
|---:|---:|---:|---:|---:|
| 1,000 | 128 µs | 645 µs | 1.30 ms | 2.65 ms |
| 10,000 | 1.67 ms | 8.39 ms | 16.8 ms | 34.3 ms |
| 100,000 | 23.8 ms | 141 ms | 292 ms | 597 ms |
| 1,000,000 | 546 ms | 2.84 s | 5.31 s | 10.8 s |
| 10,000,000 | 10.6 s | 58.5 s | 117 s | 231 s |

The total time grows linearly with the number of chromosomes: at 10M
intervals per chromosome, the time per add is 1.06 to 1.17 µs for every
chromosome count. Each chromosome has its own `StepVector`, so the chromosomes
do not slow each other down, and the chromosome lookup adds no cost that grows
with the number of chromosomes.

**N intervals in total, split across the chromosomes.** Time per add:

| Total N | 1 chromosome | 5 chromosomes | 10 chromosomes | 20 chromosomes |
|---:|---:|---:|---:|---:|
| 1,000 | 128 ns | 105 ns | 101 ns | 92 ns |
| 10,000 | 168 ns | 138 ns | 129 ns | 124 ns |
| 100,000 | 238 ns | 183 ns | 168 ns | 157 ns |
| 1,000,000 | 521 ns | 349 ns | 297 ns | 251 ns |
| 10,000,000 | 1.10 µs | 596 ns | 537 ns | 443 ns |

The same number of intervals is faster to add when it is spread over more
chromosomes, because each chromosome's map is smaller and less dense.

---

## at()

The queries are random regions, with the same distribution as the intervals.
`at()` returns every step that a query overlaps, so its cost depends on how
many steps it returns.

### Number of queries: unsorted vs sorted

A `GenomicStepVector` with 1,000,000 sorted random intervals on one
chromosome, queried with N random regions in random order or sorted by
start.

| Queries | Unsorted: total | Unsorted: per query | Sorted: total | Sorted: per query |
|---:|---:|---:|---:|---:|
| 1,000 | 2.54 ms | 2.54 µs | 2.44 ms | 2.44 µs |
| 10,000 | 46.5 ms | 4.65 µs | 35.7 ms | 3.57 µs |
| 100,000 | 556 ms | 5.56 µs | 204 ms | 2.04 µs |
| 1,000,000 | 5.51 s | 5.51 µs | 1.88 s | 1.88 µs |
| 10,000,000 | 55.1 s | 5.51 µs | 18.0 s | 1.80 µs |

The total time is **O(number of queries)** (nanobench's best fit, with 0.0%
error unsorted and 0.9% sorted). From 100,000 queries on, sorted queries are
**2.7 to 3.1× faster** than unsorted ones, probably for the same reason as
sorted adds: consecutive queries visit nearby steps that are still in cache.

### Number of entries

1,000,000 sorted random queries, against a `GenomicStepVector` built from N
sorted random intervals on one chromosome.

| Intervals in the vector (N) | Total for 1M queries | Per query |
|---:|---:|---:|
| 1,000 | 71.0 ms | 71 ns |
| 10,000 | 104 ms | 104 ns |
| 100,000 | 334 ms | 334 ns |
| 1,000,000 | 1.89 s | 1.89 µs |
| 10,000,000 | 7.78 s | 7.78 µs |

Finding the query's start and end takes two O(log N) lookups, but the
number of steps that a query overlaps grows with the density of the vector.
With 1,000 intervals over 10 Mbp, most queries overlap no interval and
return nothing (`at()` leaves out zero-valued steps by default); with 10M
intervals, a query returns hundreds of steps. The cost of `at()` is
therefore dominated by the number of steps it returns, not by the lookup.

---

## Key findings

1. **Sort intervals before adding them.** With overlapping intervals, sorted
   input was 2× faster at 100,000 intervals and 28.6× faster at 10M. In random
   order, the cost per add grows with N (O(N²) in total) and becomes limited
   by memory access.

2. **Non-overlapping adds are fast.** Sequential, non-overlapping adds cost
   130 to 674 ns each, up to 10M intervals, and scale as O(N log N).

3. **The cost of an add depends on overlap density.** Each `add()` updates
   every step boundary inside the new interval, so adds slow down as the
   intervals on a chromosome get denser.

4. **Chromosomes are independent.** The total time grows linearly with the
   number of chromosomes, and the same number of intervals is faster to add
   when it is spread over more chromosomes.

5. **Query cost follows the output size.** An `at()` query costs 71 ns in a
   sparse vector (1,000 intervals), where most queries return nothing, and
   7.78 µs in a dense one (10M intervals), where it returns hundreds of
   steps. Sorted queries are about 3× faster than unsorted ones.

---

*This document was co-authored with [Claude](https://claude.ai) (Anthropic).*
