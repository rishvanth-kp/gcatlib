# Core types

Basic types used throughout gcatlib.

## GenomicRegion

{cpp:class}`GenomicRegion` represents a genomic interval: a sequence
name, a start and end position, and a strand. As in BED, intervals are
0-based and half-open; see
[Coordinate conventions](readers.md#coordinate-conventions) for how the
other file formats compare. Writing a region to a stream with `<<` gives
its name, start, and end separated by tabs (BED3); the strand is not
written.

```{doxygenclass} GenomicRegion
:members:
```

## FeatureVector

{cpp:class}`FeatureVector` is a list of values that accumulates by
concatenation. Use it as the value type of a {cpp:class}`StepVector` or
{cpp:class}`GenomicStepVector` when each position should record which
features cover it, rather than a numeric sum.

```{doxygenclass} FeatureVector
:members:
```
