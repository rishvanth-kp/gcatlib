# Readers

Classes for reading genomic file formats:

- **SAM/BAM**: {cpp:class}`SamReader` reads alignments into
  {cpp:class}`SamEntry` records, using HTSlib. The `SamCigar`, `SamFlags`,
  and `SamTags` helpers interpret the CIGAR string, FLAG field, and
  optional tags of an entry.
- **FASTQ**: {cpp:class}`FastqReader` reads single-end or paired-end
  records into {cpp:struct}`FastqEntry`.
- **BED**: {cpp:class}`BedReader` reads regions into
  {cpp:class}`GenomicRegion`.
- **GTF**: {cpp:class}`GtfReader` reads annotation lines into
  {cpp:struct}`GtfEntry` or {cpp:struct}`GencodeGtfEntry`.

## Coordinate conventions

Each format stores coordinates differently, and the readers keep them
exactly as they are in the file:

| Format | Stored in | Convention |
|---|---|---|
| SAM/BAM | `SamEntry::pos` | 1-based |
| GTF | `GtfEntry::start`, `GtfEntry::end` | 1-based, end inclusive |
| BED | `GenomicRegion` | 0-based, half-open |

The rest of gcatlib uses {cpp:class}`GenomicRegion` (0-based, half-open),
so SAM and GTF coordinates must be converted: subtract 1 from a SAM
position or a GTF start, and keep a GTF end as it is. The `SamCigar`
helpers `reference_end_pos()` and `cigar_to_reference_regions()` already
return 0-based coordinates.

## SAM/BAM

### SamReader
```{doxygenclass} SamReader
:members:
```

### get_seq_lengths
```{doxygenfunction} get_seq_lengths
```

### SamEntry
```{doxygenclass} SamEntry
:members:
```

### SamCigar
```{doxygennamespace} SamCigar
:members:
```

### SamFlags
```{doxygennamespace} SamFlags
:members:
```

### SamTags
```{doxygennamespace} SamTags
:members:
```

## FASTQ

### FastqReader
```{doxygenclass} FastqReader
:members:
```

### FastqEntry
```{doxygenstruct} FastqEntry
:members:
```

## BED

### BedReader
```{doxygenclass} BedReader
:members:
```

## GTF

### GtfReader
```{doxygenclass} GtfReader
:members:
```

### GtfEntry
```{doxygenstruct} GtfEntry
:members:
```

### GencodeGtfEntry
```{doxygenstruct} GencodeGtfEntry
:members:
```
