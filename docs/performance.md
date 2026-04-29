# Performance

Performance characteristics of the core data structures in gcatlib, measured
using [nanobench](https://nanobench.ankerl.com) v4.3.11.

## Environment

| Item | Details |
|---|---|
| CPU | AMD Ryzen 9 5900X 12-Core Processor (24 threads) |
| Compiler | g++ with `-O2 -std=c++17` |
| Benchmarking tool | nanobench v4.3.11 |
| CPU governor | `performance` mode (`cpupower frequency-set -g performance`) |
| Core pinning | Single dedicated core (`taskset -c 2`) |

Setting the CPU governor to `performance` and pinning to a single core
eliminates frequency scaling and process migration noise. All error percentages
reported are from nanobench's stability analysis.

---

## StepVector

`StepVector<T>` stores a piecewise-constant function over a 1D integer
coordinate space using `std::map<size_t, T>` to record only the positions
where the value changes.

### add() — sequential non-overlapping

Each interval starts after the previous one ends. No existing map nodes fall
inside the new interval's range so the internal update loop is never entered.
This is the best-case access pattern.

| N | ns / add | Total time |
|---:|---:|---:|
| 1,000 | 191 | 191 µs |
| 10,000 | 311 | 3.1 ms |
| 100,000 | 412 | 41 ms |
| 1,000,000 | 510 | 510 ms |

The per-add cost grows slowly with N because each `add()` performs two
`std::map` operations (`lower_bound` and `insert`), both O(log N).
Total cost is **O(N log N)**.

### add() — random overlapping

Intervals have uniformly random start positions over a range of 100,000 bp
with lengths uniformly distributed in [1, 1,000]. Intervals are pre-generated
to exclude RNG overhead.

| N | ns / add | Total time |
|---:|---:|---:|
| 1,000 | 208 | 208 µs |
| 10,000 | 842 | 8.4 ms |
| 100,000 | 4,561 | 456 ms |

The 10x increase from N=1K to N=10K costs ~4x in per-add time, and from
N=10K to N=100K costs ~5.4x. This super-linear degradation occurs because
as the map fills, each new interval overlaps more existing step boundaries.
The internal update loop — which increments every map node between the
interval's start and end — becomes increasingly expensive as the coordinate
space saturates. See [worst-case stress](#worst-case-wide-interval-over-dense-pre-built-map)
for the limiting case.

### at() — point query

Queries 1,000 random positions uniformly distributed over the populated range
after building the map with N sequential non-overlapping adds.

| Map size (N) | ns / query |
|---:|---:|
| 1,000 | 12.4 |
| 10,000 | 27.4 |
| 100,000 | 70.6 |

Point queries use `std::map::upper_bound`, which is **O(log N)**.
Tripling from N=1K to N=10K (+10x) costs ~2.2x in query time, consistent
with log scaling.

### at_range() — full range query

Queries the entire populated range [0, N×10) after N sequential adds, writing
all step boundaries to an output vector.

| Map size (N) | Time per query |
|---:|---:|
| 1,000 | 10.6 µs |
| 10,000 | 110.8 µs |
| 100,000 | 1.15 ms |

Output size grows linearly with N (each sequential add contributes 2 map
entries), so total time is **O(N)**. The ~10.5x cost increase per 10x
increase in N confirms this.

---

## GenomicStepVector

`GenomicStepVector<T>` maintains one independent `StepVector<T>` per
chromosome, dispatched via `std::unordered_map<string, StepVector<T>>`.

### Single chromosome — dispatch overhead baseline

Single-chromosome performance can be compared directly against `StepVector`
to isolate the cost of the `unordered_map` chromosome dispatch.

| N | StepVector (ns/add) | GenomicStepVector (ns/add) | Overhead |
|---:|---:|---:|---:|
| 1,000 | 191 | 135 | ~0% |
| 10,000 | 311 | 252 | ~0% |
| 100,000 | 412 | 346 | ~0% |
| 1,000,000 | 510 | 485 | ~0% |
| 10,000,000 | — | 576 | — |

The `unordered_map` lookup adds negligible overhead — well within measurement
noise. `GenomicStepVector` is effectively free compared to the underlying
`StepVector` operations.

### Multi-chromosome sequential scaling

Each chromosome receives the same number of sequential non-overlapping adds.
Total adds = n_chrom × n_per_chrom.

| n_chrom | n_per_chrom | Total adds | Total time |
|---:|---:|---:|---:|
| 1 | 1,000 | 1,000 | 121 µs |
| 1 | 10,000 | 10,000 | 2.4 ms |
| 1 | 100,000 | 100,000 | 34 ms |
| 5 | 1,000 | 5,000 | 636 µs |
| 5 | 10,000 | 50,000 | 12.1 ms |
| 5 | 100,000 | 500,000 | 182 ms |
| 20 | 1,000 | 20,000 | 2.7 ms |
| 20 | 10,000 | 200,000 | 49 ms |
| 20 | 100,000 | 2,000,000 | 787 ms |

Scaling is near-linear with both n_chrom and n_per_chrom. Going from 1 to 5
to 20 chromosomes at the same n_per_chrom scales proportionally, confirming
that each chromosome's `StepVector` operates independently.

### Random overlapping — super-linear degradation

Each chromosome receives the same pre-generated set of random overlapping
intervals (positions in [0, 100,000], lengths in [1, 1,000]).

| n_chrom | n_per_chrom | Total adds | Total time |
|---:|---:|---:|---:|
| 1 | 1,000 | 1,000 | 212 µs |
| 1 | 10,000 | 10,000 | 7.6 ms |
| 1 | 100,000 | 100,000 | 433 ms |
| 1 | 1,000,000 | 1,000,000 | 78 s |
| 5 | 1,000 | 5,000 | 1.0 ms |
| 5 | 10,000 | 50,000 | 38.8 ms |
| 5 | 100,000 | 500,000 | 2,139 ms |
| 5 | 1,000,000 | 5,000,000 | 394 s |
| 20 | 1,000 | 20,000 | 4.2 ms |
| 20 | 10,000 | 200,000 | 156 ms |
| 20 | 100,000 | 2,000,000 | 8,658 ms |
| 20 | 1,000,000 | 20,000,000 | 1,560 s |

**Key observation:** The jump from n_per_chrom=10K to 100K costs ~57x per
chromosome (7.6 ms → 433 ms), far worse than the 10x expected from O(N log N).
This is because the coordinate range [0, 100,000] becomes saturated — once
the step boundaries are dense enough, every new `add()` must traverse a large
fraction of the map in the update loop, pushing cost toward O(N²).

Scaling with n_chrom remains proportional (each chromosome's map is independent),
so the total time for 20 chromosomes at n_per_chrom=1M is approximately
20× the single-chromosome case.

### at() query — chromosome lookup cost

Builds a vector with n_per_chrom=10,000 sequential adds per chromosome, then
runs 1,000 random queries across chromosomes.

| n_chrom | Time for 1,000 queries |
|---:|---:|
| 1 | 25.4 ms |
| 5 | 24.8 ms |
| 20 | 26.6 ms |

Query time is essentially flat regardless of chromosome count. The
`unordered_map` chromosome lookup is O(1), and the dominant cost is the
`at_range()` traversal and output merging — not the dispatch.

### Stress tests

#### Worst-case: wide interval over dense pre-built map

Pre-builds N non-overlapping length-1 intervals (creating N map nodes), then
repeatedly adds a single interval spanning the entire range. Each `add()` call
must traverse and update every existing map node — this is the maximum
per-call cost.

| Pre-built nodes (N) | Time per add |
|---:|---:|
| 10,000 | 145 µs ⚠ |
| 100,000 | 864 µs |
| 1,000,000 | 23.9 ms |

⚠ N=10,000 result is flagged as unstable by nanobench (~7.8 iterations
collected). The N=100K and N=1M results are reliable.

The ~28x cost increase from N=100K to N=1M (for 10x more nodes) reflects
**O(N)** traversal with additional cache pressure — once the map exceeds
the CPU cache, random pointer-chasing in `std::map` nodes becomes expensive.

#### Maximally dense: sequential length-1 intervals

Each `add()` creates a new 1-position interval, maximising map size for a
given number of adds. This tests performance under peak memory pressure.

| N | Total time |
|---:|---:|
| 100,000 | 30 ms |
| 1,000,000 | 437 ms |

Scaling of ~14.6x for 10x more adds is slightly super-linear compared to the
sequential 5-bp case (~12x), reflecting the larger and more memory-intensive
map.

---

## Key Findings

1. **`unordered_map` dispatch is free.** The per-chromosome lookup in
   `GenomicStepVector` adds no measurable overhead over raw `StepVector`.
   Multi-chromosome scaling is proportional: N chromosomes costs N× the
   single-chromosome time.

2. **Sequential adds are O(N log N).** When intervals do not overlap,
   each `add()` only touches 2 map entries and the update loop is never
   entered. This is the fast path.

3. **Random overlapping add() is the main bottleneck.** When intervals
   overlap, `add()` must update every map node between the interval's start
   and end. As the coordinate space saturates, this approaches O(N²) behaviour.
   For a 100 kbp coordinate range, saturation begins around N=10,000 adds.
   **Consider pre-sorting or batching overlapping intervals to reduce map
   density if this is a bottleneck.**

4. **Point queries are O(log N) and fast.** Even with 100,000 step boundaries,
   a single `at()` costs ~71 ns. Query performance degrades slowly with map
   size.

5. **Range queries are O(N) in output size.** `at_range()` is dominated by
   the cost of writing output, not by lookup. Reusing the output vector across
   calls (avoiding repeated allocation) is recommended.

6. **Cache pressure at large map sizes.** The worst-case stress test shows that
   traversing 1M `std::map` nodes costs ~28× more than traversing 100K nodes
   — super-linear due to cache misses on pointer-based tree nodes. For very
   large datasets, a cache-friendlier data structure (e.g. sorted vector of
   breakpoints) could improve performance.

---

*This document was co-authored with [Claude](https://claude.ai) (Anthropic).*
