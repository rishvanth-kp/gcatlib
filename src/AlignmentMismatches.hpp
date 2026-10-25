/*
* AlignmentMismatches: class for keeping track of mismatches
* Copyright (C) 2022 Rishvanth Prabakar
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License along
* with this program; if not, write to the Free Software Foundation, Inc.,
* 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef ALIGNMENT_MISMATCHES_HPP
#define ALIGNMENT_MISMATCHES_HPP

#include <string>
#include <vector>

#include "SamEntry.hpp"
#include "GenomicStepVector.hpp"

/**
* \brief Keep track of mismatched bases and their location.
*
* Parses the CIGAR string and MD tag of each added SAM entry and tracks,
* per genomic position, the read bases that differ from the reference
* and the alignment coverage. For each position, the mismatched read
* bases from all added entries are joined into a single string in the
* order the entries were added (e.g. "AG" if one read had an A and a
* later read a G). Reference bases are not stored, and insertions and
* deletions are not recorded as mismatches.
*
* Coverage counts the read bases aligned to a position by M, = or X
* CIGAR operations; deletions (D) and skipped regions (N) do not add
* coverage.
*
* @note No filtering is done: every added entry contributes. Strand is
*   not considered, so mismatches from forward and reverse reads are
*   pooled. Secondary, supplementary, duplicate, and low-MAPQ alignments
*   and base qualities are not checked, and overlapping mates of a pair
*   are both counted. Filter entries before calling add() as needed.
*/
class AlignmentMismatch {
public:
  /**
  * Default constructor.
  */
  AlignmentMismatch();
  /**
  * Default destructor.
  */
  ~AlignmentMismatch();

  /**
  * Adds the coverage and mismatches of a SAM entry.
  * The CIGAR string is parsed, and all matching and mismatching bases
  * are added to coverage. Then the MD tag is parsed to find the
  * mismatches, and the read base at each mismatch is recorded. The
  * reference base is not recorded.
  *
  * The entry must be a mapped alignment with an MD tag.
  *
  * @param [in] sam_entry SAM entry to be added.
  * @throws std::invalid_argument if the entry has no MD tag. This
  *   includes all unmapped reads, so filter them out before calling.
  *   Coverage from the entry has already been added when this is thrown.
  */
  void add(const SamEntry &sam_entry);

  /**
  * After all the SAM entries are added, this method is used to retrieve
  * the mismatched bases and the coverage at each mismatched position.
  * Only positions with at least one mismatch are returned, in order of
  * position. All output vectors are cleared first, and are left empty
  * if the chromosome of in_region has no recorded mismatches.
  *
  * @param [in] in_region Region to get the mismatches from.
  * @param [out] out_region One-base-wide region for each mismatched
  *   position in the input region. 'aln_count' and 'mismatches' are
  *   parallel arrays corresponding to the locations in 'out_region'.
  * @param [out] aln_count Coverage at the position.
  * @param [out] mismatches A string with all the mismatched bases at the
  *   position. The length of this string is equal to the number of
  *   mismatched bases.
  * @warning in_region.start must be less than in_region.end. As with
  *   GenomicStepVector::at(), a query with start >= end on a chromosome
  *   that has mismatches reads past the end of an internal vector
  *   (undefined behaviour, typically a crash).
  */
  void at(const GenomicRegion &in_region,
          vector<GenomicRegion> &out_region,
          vector<size_t> &aln_count,
          vector<string> &mismatches) const;

private:
  GenomicStepVector<string> mismatches;
  GenomicStepVector<size_t> coverage;
};

#endif
