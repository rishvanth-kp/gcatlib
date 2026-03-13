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
* \brief A class for holding values in a SAM entry.
*
* Takes an entry read from a SAM/BAM file as a string and parses
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
  *
  * @param [in] line SAM/BAM line to parse.
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
  uint16_t flag; /**< FLAG. Bitwise flag. See \ref SamFlags for interpretation. */
  string rname; /**< RNAME. Reference sequence name. */
  uint32_t pos; /**< POS. 1-based leftmost mapping position. */
  uint16_t mapq; /**< MAPQ. Mapping quality. 255 indicates not available. */
  string cigar; /**< CIGAR. CIGAR string. */
  string rnext; /**< RNEXT. Reference name of mate. */
  uint32_t pnext; /**< PNEXT. Position of mate. */
  int tlen; /**< TLEN. Observed template length. */
  string seq; /**< SEQ. Segment sequence. */
  string qual; /**< QUAL. ASCII of Phred-scaled base quality. */

  vector<string> tags /**< TAGS. Optional tags. */;
};

/**
* Writes a SamEntry to an output stream in SAM format, with fields
* tab-separated in the standard SAM column order.
*
* @param [in,out] os Output stream to write to.
* @param [in]     e  SamEntry to write.
* @return Reference to the output stream.
*/
std::ostream& operator<< (std::ostream&, const SamEntry&);

/**
* \brief Helper functions for working with SAM CIGAR strings.
*
* CIGAR (Compact Idiosyncratic Gapped Alignment Report) strings encode
* the alignment of a read to a reference sequence as a series of
* operations (e.g. match, insertion, deletion). These utilities parse
* and interpret CIGAR strings from \ref SamEntry objects.
*/
namespace SamCigar {

  /**
  * CIGAR operation codes as defined in the SAM specification.
  * Each value corresponds to the character used in the CIGAR string.
  */
  enum class Cigar : char {
    aln_match    = 'M', /**< Alignment match (can be match or mismatch). */
    ref_insert   = 'I', /**< Insertion relative to the reference. */
    ref_del      = 'D', /**< Deletion from the reference. */
    ref_skip     = 'N', /**< Skipped region from the reference (e.g. intron). */
    soft_clip    = 'S', /**< Soft clipping (clipped bases present in SEQ). */
    hard_clip    = 'H', /**< Hard clipping (clipped bases absent from SEQ). */
    padding      = 'P', /**< Padding (silent deletion from padded reference). */
    seq_match    = '=', /**< Sequence match. */
    seq_mismatch = 'X'  /**< Sequence mismatch. */
  };

  /** A vector of (Cigar operation, length) pairs representing a parsed
  *  CIGAR string. */
  using CigarTuples = vector<pair<Cigar, size_t>>;

  /** A vector of (Cigar operation, GenomicRegion) pairs mapping each
  *  CIGAR operation to its reference-coordinate span. */
  using CigarRegions = vector<pair<Cigar, GenomicRegion>>;

  /**
  * Takes a SamEntry and converts the cigar string into CigarTuples.
  * Each cigar operation is converted into a pair of operation and length,
  * and inserted into CigarTuples. The length of CigarTuples is equal to
  * the number of operations in the cigar string.
  *
  * @param [in]  e      SamEntry containing the cigar string to be parsed.
  * @param [out] tuples CigarTuples containing the parsed cigar operations.
  *   Cleared before populating.
  */
  void
  string_to_tuple(const SamEntry &e, CigarTuples &tuples);

  /**
  * Advances ref_pos and query_pos by len bases consumed in the reference.
  * Alignment match (M, =, X) and deletion (D) operations consume the
  * reference and count toward len. Insertion (I) and soft-clip (S)
  * consume only the query and do not count toward len. Intron/skip (N)
  * operations consume the reference but, when skipN is true (default),
  * do not count toward len.
  *
  * @param [in]     tuples    Parsed CIGAR operations from string_to_tuple().
  * @param [in]     len       Number of reference bases to advance.
  * @param [in,out] ref_pos   Reference position to advance (0-based).
  * @param [in,out] query_pos Query position to advance (0-based).
  * @param [in]     skipN     If true (default), N operations advance the
  *   reference position but do not count toward len.
  */
  void
  move_in_reference(const CigarTuples &tuples, const size_t len,
    size_t &ref_pos, size_t &query_pos, bool skipN = true);

  /**
  * Advances ref_pos and query_pos by len bases consumed in the query.
  *
  * @param [in]     tuples    Parsed CIGAR operations from string_to_tuple().
  * @param [in]     len       Number of query bases to advance.
  * @param [in,out] ref_pos   Reference position to advance (0-based).
  * @param [in,out] query_pos Query position to advance (0-based).
  *
  * @note Not yet implemented.
  */
  void
  move_in_query(const CigarTuples &tuples, const size_t len,
    size_t &ref_pos, size_t &query_pos);

  /**
  * Returns the 0-based exclusive end position of the alignment in the
  * reference (i.e., one past the last aligned base). Accounts for all
  * CIGAR operations that consume the reference (M, =, X, D, N).
  *
  * @param [in] e SamEntry containing the alignment.
  * @return 0-based exclusive reference end position.
  */
  size_t
  reference_end_pos(const SamEntry &e);

  /**
  * Takes a SamEntry with a cigar string and converts it to a vector of
  * GenomicRegions for each cigar operation. Each region corresponds to
  * the reference span of the respective cigar operation. Operations that
  * do not consume the reference (I, S, H, P) produce a zero-length region
  * where start equals end.
  *
  * @param [in]  e           SamEntry containing the cigar string to be parsed.
  * @param [out] ref_regions A vector of pairs of cigar operation and
  *   GenomicRegion, one per cigar operation. Regions are half-open
  *   intervals with 0-based offsets. Cleared before populating.
  */
  void
  cigar_to_reference_regions(const SamEntry &e,
    CigarRegions &ref_regions);
}

/**
* \brief Helper functions for working with the SAM FLAG field.
*
* The SAM FLAG field is a 16-bit integer where each bit encodes a
* property of the alignment. The \ref Flag enum defines the individual
* bit values, and the functions in this namespace provide convenient
* ways to test and set them.
*/
namespace SamFlags {

  /**
  * Bit flags for the SAM FLAG field, as defined in the SAM specification.
  */
  enum class Flag : uint16_t {
    read_paired       = 0x0001, /**< Read is paired. */
    proper_pair       = 0x0002, /**< Read is mapped in a proper pair. */
    read_unmapped     = 0x0004, /**< Read is unmapped. */
    mate_unmapped     = 0x0008, /**< Mate is unmapped. */
    read_reverse      = 0x0010, /**< Read is mapped to the reverse strand. */
    mate_reverse      = 0x0020, /**< Mate is mapped to the reverse strand. */
    first_in_pair     = 0x0040, /**< Read is the first in the pair. */
    second_in_pair    = 0x0080, /**< Read is the second in the pair. */
    not_primary_aln   = 0x0100, /**< Not a primary alignment. */
    fail_qc           = 0x0200, /**< Read fails platform/vendor quality checks. */
    pcr_duplicate     = 0x0400, /**< Read is a PCR or optical duplicate. */
    supplementary_aln = 0x0800  /**< Supplementary alignment. */
  };

  /**
  * Returns true if any of the bits in check are set in flag.
  *
  * @param [in] flag  The FLAG value from a \ref SamEntry.
  * @param [in] check Bitmask of flags to test.
  * @return True if any bit in check is set in flag.
  */
  constexpr bool is_any_set(const uint16_t flag, const uint16_t check) {
    return (flag & check);
  }

  /**
  * Returns true if all of the bits in check are set in flag.
  *
  * @param [in] flag  The FLAG value from a \ref SamEntry.
  * @param [in] check Bitmask of flags to test.
  * @return True if every bit in check is set in flag.
  */
  constexpr bool is_all_set(const uint16_t flag, const uint16_t check) {
    return ((flag & check) == check);
  }

  /**
  * Returns true if the given Flag bit is set in flag.
  *
  * @param [in] flag The FLAG value from a \ref SamEntry.
  * @param [in] f    The \ref Flag bit to test.
  * @return True if flag f is set.
  */
  constexpr bool is_set(const uint16_t flag, Flag f) {
    return (static_cast<uint16_t>(f) & flag);
  }

  /**
  * Returns a new flag value with the given Flag bit set.
  *
  * @param [in] flag The original FLAG value.
  * @param [in] f    The \ref Flag bit to set.
  * @return New FLAG value with bit f set.
  */
  constexpr uint16_t set(const uint16_t flag, Flag f) {
    return (static_cast<uint16_t>(f) | flag);
  }

}

/**
* \brief Helper functions for working with SAM optional tags.
*
* SAM optional tags are extra fields appended to each alignment record
* in TAG:TYPE:VALUE format. This namespace provides utilities for
* retrieving tag values and parsing specific tags such as the MD string.
*/
namespace SamTags {

/**
* Searches for a tag in a vector of SAM optional tags and returns its
* value. Tags are expected to be in TAG:TYPE:VALUE format, where TAG is
* 2 characters and TYPE is 1 character.
*
* @param [in]  tags  Vector of optional tag strings from a \ref SamEntry.
* @param [in]  tag   Two-character tag name to search for (e.g. "MD").
* @param [out] value The value portion of the tag if found. Cleared
*   before populating.
* @return True if the tag was found. False otherwise.
*/
bool
get_tag(const vector<string> &tags, const string &tag, string &value);

/**
* Parses an MD tag string into a vector of (offset, base) tuples.
* Each tuple contains the number of matching reference bases preceding
* an event, and a string for the event: a mismatched reference base
* (e.g. "A") or a deletion prefixed with '^' (e.g. "^ACG"). The final
* tuple always has an empty string for the event.
*
* @param [in]  md_tag The MD tag value string.
* @param [out] tuples Parsed tuples of (reference_offset, event_string).
*   Cleared before populating.
*/
void
md_to_tuple(const string &md_tag, vector<pair<size_t, string>> &tuples);

}

#endif
