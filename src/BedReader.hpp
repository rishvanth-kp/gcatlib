/*
* BedReader: class to read bed files
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

#ifndef BED_READER_HPP
#define BED_READER_HPP

#include <vector>
#include <string>
#include <fstream>

#include "GenomicRegion.hpp"


/**
* \brief BED file reader.
*
* A class to read and parse BED files. Supports reading BED3 entries
* (chrom, start, end) and variable-column BED entries where additional
* columns are returned as a vector of strings. Coordinates are stored
* as read, so they follow the BED convention used by \ref GenomicRegion
* (0-based, half-open).
*
* @note Every line is treated as a data line with at least three
*   columns. Blank lines and track, browser, or '#' comment lines are
*   not skipped, and malformed lines are not reported: they are returned
*   as records with bogus values (e.g. a non-numeric start or end is
*   read as 0, or left unchanged). Remove such lines before reading.
*/
class BedReader {
public:
  /**
  * Opens a BED file.
  *
  * @param [in] in_file BED file name.
  * @throws std::runtime_error if the file cannot be opened.
  */
  BedReader(const string &in_file);
  /**
  * Closes the BED file.
  */
  ~BedReader();

  /**
  * Reads the next BED3 (chrom, start, end) line from the file
  * and populates a \ref GenomicRegion. Columns must be separated by
  * tabs; any columns after the third are ignored. Only g.name,
  * g.start, and g.end are set; g.strand keeps its previous value.
  *
  * @param [out] g GenomicRegion to populate.
  * @return True on successfully reading a BED3 entry. False if
  *   end of file is reached.
  * @warning A line with fewer than three tab-separated columns,
  *   including a blank or space-separated line, causes undefined
  *   behaviour (an out-of-range read of the split columns).
  */
  bool read_bed3_line(GenomicRegion &g);

  /**
  * Reads the next BED line and populates a \ref GenomicRegion with
  * the first 3 columns (chrom, start, end). Any remaining columns
  * are returned in a vector of strings. If the line has only 3
  * columns, the fields vector will be empty.
  *
  * Columns are split on any whitespace (tabs or spaces), so a field
  * that itself contains spaces is split into several entries. Only
  * g.name, g.start, and g.end are set; g.strand keeps its previous
  * value. For a BED6 line, the extra columns are fields[0] = name,
  * fields[1] = score, and fields[2] = strand.
  *
  * @param [out] g      GenomicRegion to populate with chrom, start, end.
  * @param [out] fields Vector of strings containing any additional
  *   columns beyond the first three. Cleared before populating.
  * @return True on successfully reading a BED entry. False if end of
  *   file is reached.
  */
  bool read_bed_line(GenomicRegion &g, std::vector<std::string> &fields);

  /**
  * Reads the rest of the BED file and appends a \ref GenomicRegion for
  * each line, with the BED3 (chrom, start, end) fields. Existing
  * contents of g are kept. Reading starts from the current position, so
  * if read_bed3_line() or read_bed_line() has already been called, only
  * the remaining lines are read. Columns must be separated by tabs; any
  * columns after the third are ignored. Appended regions have strand '.'.
  *
  * @param [in,out] g Vector of GenomicRegions to append to.
  * @warning As with read_bed3_line(), a line with fewer than three
  *   tab-separated columns, including a blank or space-separated line,
  *   causes undefined behaviour.
  */
  void read_bed3_file(std::vector<GenomicRegion> &g);

private:
  std::ifstream in;

};

#endif
