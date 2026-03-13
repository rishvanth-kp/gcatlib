/*
* SamReader: class to read SAM files
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

#ifndef SAM_READER_HPP
#define SAM_READER_HPP

#include <iostream>
#include <string>
#include <unordered_map>

#include "SamEntry.hpp"
#include "GenomicRegion.hpp"

#include "htslib/hts.h"
#include "htslib/sam.h"

using std::string;
using std::unordered_map;

/**
* \brief SAM/BAM file reader.
*
* Uses htslib to read SAM/BAM header and entries. This class only reads
* the SAM header and entries, it uses \ref SamEntry to parse them.
*
* @see SamEntry
*/
class SamReader {
public:
  /**
  * Opens a SAM/BAM file, verifies that it is a sequence data file,
  * requires a header, and initializes htslib handlers.
  *
  * @param [in] in_file SAM/BAM file name.
  * @throws std::runtime_error if the file cannot be opened, is not a
  *   sequence data file, has no header, or htslib initialization fails.
  */
  SamReader(const string &in_file);
  /**
  * Closes the file and destroys htslib handlers.
  */
  ~SamReader();

  /**
  * Reads one entry and populates 'entry' with the help of \ref SamEntry.
  *
  * @param [out] entry SamEntry to populate with the parsed fields.
  *   Contents are not valid if false is returned.
  * @return True on successfully reading an entry. False if end of file
  *   is reached.
  * @throws std::runtime_error if the entry cannot be parsed.
  */
  bool read_sam_line(SamEntry &entry);

  /**
  * Reads a pair of entries from a paired-end SAM/BAM file, and populates
  * 'entry1' and 'entry2' in the order the pair appears in the SAM file.
  * This uses a hash table to keep track of SAM entries whose mates have
  * not yet been read. SAM entries are read using 'read_sam_line', and
  * checks if another entry with the same QNAME exists in the hash table.
  * If it exists, the current read entry and the entry in the hash table
  * are returned. If it does not exist, the current entry is stored in the
  * hash table.
  * It is assumed that each read pair contains exactly two entries in the
  * SAM file. Having more than two entries present for a read pair (for
  * example, supplementary alignments) will result in returning more than one
  * read pair for the same QNAME (if there are even entries) or having dangling
  * entries in the hash table (if there are odd entries).
  *
  * @param [out] entry1 SamEntry with the first read
  * @param [out] entry2 SamEntry with the second read
  * @return True on successfully reading a pair. False if end of file is
  * reached.
  */
  bool read_pe_sam(SamEntry &entry1, SamEntry &entry2);

  /**
  * Returns the entire SAM header as a string, including all header
  * lines (@HD, @SQ, @RG, @PG, etc.).
  *
  * @param [out] hdr String to hold the header.
  */
  void read_sam_header(string &hdr);

private:
  htsFile *hts;
  sam_hdr_t *header;
  bam1_t *b;
  kstring_t sam_entry;
  bool eof;

  unordered_map<string, SamEntry> lonley_mates;
};

/**
* Parses a SAM header string and extracts reference sequence names
* and lengths from @SQ header lines. Each sequence is returned as a
* \ref GenomicRegion with name set to the sequence name (SN tag) and
* end set to the sequence length (LN tag).
*
* @param [in]  hdr SAM header string, as returned by read_sam_header().
* @param [out] out Vector of GenomicRegions, one per reference sequence.
*   Cleared before populating.
*/
void
get_seq_lengths(const string &hdr, vector<GenomicRegion> &out);

#endif
