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
* columns are returned as a vector of strings.
*/
class BedReader {
public:
  /**
  * Opens a BED file. Throws a runtime error if the file cannot
  * be opened.
  *
  * @param [in] in_file BED file name.
  */
  BedReader(const string &in_file);
  /**
  * Closes the BED file.
  */
  ~BedReader();

  /**
  * Reads the next BED3 (chrom, start, end) line from the file
  * and populates a \ref GenomicRegion.
  *
  * @param [out] g GenomicRegion to populate.
  * @return True on successfully reading a BED3 entry. False if
  *   end of file is reached.
  */
  bool read_bed3_line(GenomicRegion &g);

  /**
  * Reads the next BED line and populates a \ref GenomicRegion with
  * the first 3 columns (chrom, start, end). Any remaining columns
  * are returned in a vector of strings. If the line has only 3
  * columns, the fields vector will be empty.
  *
  * @param [out] g      GenomicRegion to populate with chrom, start, end.
  * @param [out] fields Vector of strings containing any additional
  *   columns beyond the first three.
  * @return True on successfully reading a BED entry. False if end of
  *   file is reached.
  */
  bool read_bed_line(GenomicRegion &g, std::vector<std::string> &fields);

  /**
  * Reads the entire BED file and populates a vector of \ref GenomicRegion
  * with the BED3 (chrom, start, end) fields.
  *
  * @param [out] g Vector of GenomicRegions to populate.
  */
  void read_bed3_file(std::vector<GenomicRegion> &g);

private:
  std::ifstream in;

};

#endif
