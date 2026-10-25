# Analysis

Classes that summarize aligned reads against genomic features:

- {cpp:class}`AlignedGenomicFeature` counts aligned bases per genomic
  feature category (e.g. CDS, UTR, intron, intergenic), for a bulk
  sample or per cell barcode.
- {cpp:class}`AlignedOverlapGenomicFeature` counts paired-end fragments
  per cell barcode against features from GTF and BED files, which may
  overlap.
- {cpp:class}`AlignmentMismatch` records mismatched bases and coverage
  at each genomic position.

See also [Metagene](metagene.md), which maps positions to normalized
positions along features for metagene plots.

## Choosing a feature counter

| | `AlignedGenomicFeature` | `AlignedOverlapGenomicFeature` |
|---|---|---|
| Input | single alignments | read pairs (fragments) |
| Counts | aligned bases | fragments, fragment bases, mean fragment length |
| Features | one category per position, from GTF | overlapping features, from GTF and BED |
| Bulk / barcoded | both; optional UMI dedup | barcoded only |
| Typical use | e.g. RNA-seq QC | e.g. ATAC-seq |

## AlignedGenomicFeature

Assigns each aligned base to exactly one feature category.

```{doxygenclass} AlignedGenomicFeature
:members:
```

## AlignedOverlapGenomicFeature

Assigns each fragment to every feature at its left end, allowing
features to overlap.

```{doxygenclass} AlignedOverlapGenomicFeature
:members:
```

## AlignmentMismatch

Records the read bases that differ from the reference, with coverage.

```{doxygenclass} AlignmentMismatch
:members:
```
