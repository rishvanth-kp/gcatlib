# Metagene

{cpp:class}`Metagene` maps genomic positions to normalized positions
along features, for metagene plots. It reads features (e.g. transcripts)
from a BED6 file and divides each one into a fixed number of divisions,
numbered from the 5' end to the 3' end, so that signal from features of
different lengths can be combined on a common scale. See the class
documentation below for the BED file requirements.

See also [Analysis](analysis.md) for the other classes that summarize
aligned reads against genomic features.

## Metagene
```{doxygenclass} Metagene
:members:
```

