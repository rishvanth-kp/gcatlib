/*
* SamEntry: class to store SAM entry
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

#ifndef SAM_ENTRY_HPP
#define SAM_ENTRY_HPP

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <cstdint>

#include "GenomicRegion.hpp"

using std::string;
using std::string;
using std::vector;
using std::pair;

/**
* \brief A class for holing values in a SAM entry.
*
* Takes an entry read from a SAM/BAM file as a string and parsers
* the mandatory and optional fields.
*/
class SamEntry {
public:

  /**
  * Default constructor. Initializes SAM fields to default values.
  */
  SamEntry(): flag{0}, pos{0}, mapq{255}, pnext{0}, tlen{0} {};
  /**
  * Constructs from a SAM/BAM string by calling parse_entry.
  */
  SamEntry(const string &line);

  /**
  * Default destructor.
  */
  ~SamEntry();

  
  /**
  * Parses a string into SAM/BAM fields.
  * The required fields are stored in the appropriate member.
  * The optional tags are stored as a TAG:TYPE:VALUE string in
  * a vector.
  *
  * @param [in] line SAM/BAM line.
  */
  void parse_entry(const string &line);

  string qname; /**< QNAME. Query template name. */
  uint16_t flag; /**< FLAG. Bitwise flag. */
  string rname; /**< RNAME. Reference sequence name. */
  uint32_t pos; /**< POS. 1-based leftmost mapping position. */
  uint16_t mapq; /**< MAPQ. Mapping quality. */
  string cigar; /**< CIGAR. CIGAR string. */
  string rnext; /**< RNEXT. Reference name of mate. */
  uint32_t pnext; /**< PNEXT. Position of mate. */
  int tlen; /**< TLEN. Observed template length. */
  string seq; /**< SEQ. Segment sequence. */
  string qual; /**< QUAL. ASCII if Phred-scaled base quality. */

  vector<string> tags /**< TAGS. Optional tags. */;
};

std::ostream& operator<< (std::ostream&, const SamEntry&);

/**
* \brief Helper functions for working with SAM CIGAR.
*
*
*/
namespace SamCigar {

  enum class Cigar : char {
    aln_match = 'M',
    ref_insert = 'I',
    ref_del = 'D',
    ref_skip = 'N',
    soft_clip = 'S',
    hard_clip = 'H',
    padding = 'P',
    seq_match = '=',
    seq_mismatch = 'X'
  };

  using CigarTuples = vector<pair<Cigar, size_t>>;
  using CigarRegions = vector<pair<Cigar, GenomicRegion>>;

  /**
  * Takes a SamEntry and converts to cigar stirng into CigarTuples.
  * Each cigar operation is coverted into a pair of operation and length,
  * and inserted into CigarTuples. The lenght of CiagrTuples is equal to the
  * number of operations in the cigar string.
  *
  * @param [in] e SamEntry containin the cigar string to be parsed.
  * @param [out] tuples CigarTuples containing the parsed cigar operations.
  */
  void
  string_to_tuple(const SamEntry &e, CigarTuples &tuples);

  void
  move_in_reference(const CigarTuples &tuples, const size_t len,
    size_t &ref_pos, size_t &query_pos, bool skipN = true);

  void
  move_in_query(const CigarTuples &tuples, const size_t len,
    size_t &ref_pos, size_t &query_pos);

  size_t
  reference_end_pos(const SamEntry &e);

  /**
  * Takes a SamEntry with a cigar stirng, and converts
  * them to a vector of GenomicRegions for each cigar operation.
  * Each region corresponds to the reference regions for the respecive
  * cigar operation.
  *
  * @param [in] e SamEntry containing the cigar string to be parsed.
  * @param [out] ref_regions A vector of a pair of cigar op and GenomicRegions
  *   with a regions corresponding to each cigar operation. The regions are
  *   represented as half open intervals with 0-based offsets.
  */
  void
  cigar_to_reference_regions(const SamEntry &e,
    CigarRegions &ref_regions);
}

namespace SamFlags {

  enum class Flag : uint16_t {
    read_paired = 0x0001,
    proper_pair = 0x0002,
    read_unmapped = 0x0004,
    mate_unmapped = 0x0008,
    read_reverse = 0x0010,
    mate_reverse = 0x0020,
    first_in_pair = 0x0040,
    second_in_pair = 0x0080,
    not_primary_aln = 0x0100,
    fail_qc = 0x0200,
    pcr_duplicate = 0x0400,
    supplementary_aln = 0x0800
  };

  constexpr bool is_any_set(const uint16_t flag, const uint16_t check) {
    return (flag & check);
  }

  constexpr bool is_all_set(const uint16_t flag, const uint16_t check) {
    return ((flag & check) == check);
  }

  constexpr bool is_set(const uint16_t flag, Flag f) {
    return (static_cast<uint16_t>(f) & flag);
  }

  constexpr uint16_t set(const uint16_t flag, Flag f) {
    return (static_cast<uint16_t>(f) | flag);
  }

}

namespace SamTags {
bool
get_tag(const vector<string> &tags, const string &tag, string &value);

void
md_to_tuple(const string &md_tag, vector<pair<size_t, string>> &tuples);

}

#endif
