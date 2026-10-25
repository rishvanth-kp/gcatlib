/*
* GenomicRegion: class to store genomic intervals
* Copyright (C) 2021 Rishvanth Prabakar
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

#ifndef GENOMIC_REGION_HPP
#define GENOMIC_REGION_HPP

#include <iostream>
#include <string>
#include <limits>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

/**
* \brief A genomic interval on a named reference sequence.
*
* GenomicRegion stores a sequence (chromosome) name, a start and end
* position, and a strand. Following the BED convention, intervals are
* 0-based and half-open: [start, end) covers positions start through
* end - 1, and its length is end - start.
*
* Other formats use different conventions and must be converted before
* constructing a GenomicRegion: SAM positions (\ref SamEntry::pos) are
* 1-based, and GTF coordinates (\ref GtfEntry) are 1-based and inclusive.
*/
class GenomicRegion {
public:
  /**
  * Constructs a region from its components. All arguments have
  * defaults, so GenomicRegion() creates an unnamed region spanning
  * [0, SIZE_MAX) with unknown strand. An end of SIZE_MAX denotes
  * "to the end of the sequence", and is used by \ref GenomicStepVector
  * to query an entire chromosome.
  *
  * @param [in] n  Sequence name (e.g. "chr1"). Defaults to "".
  * @param [in] s  Start position (inclusive, 0-based). Defaults to 0.
  * @param [in] e  End position (exclusive, 0-based). Defaults to SIZE_MAX.
  * @param [in] st Strand: '+', '-', or '.' if unknown. Defaults to '.'.
  */
  explicit GenomicRegion(string n="", size_t s=0,
                size_t e=std::numeric_limits<size_t>::max(),
                char st='.') {
    name = n;
    start = s;
    end = e;
    strand = st;
  }

  string name;  /**< Sequence (chromosome) name. */
  size_t start; /**< Start position (inclusive, 0-based). */
  size_t end;   /**< End position (exclusive, 0-based). */
  char strand;  /**< Strand: '+', '-', or '.' if unknown. */

};

/**
* Writes a GenomicRegion to an output stream as tab-separated name,
* start, and end (i.e. BED3 format). The strand is not written.
*
* @param [in,out] out Output stream to write to.
* @param [in]     g   GenomicRegion to write.
* @return Reference to the output stream.
*/
std::ostream&
operator<<(std::ostream &out, const GenomicRegion &g);

#endif
