/*
* FastqReader: class to read fastq files
* Copyright (C) 2023 Rishvanth Prabakar
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

#ifndef FASTQ_READER_HPP
#define FASTQ_READER_HPP

#include <iostream>
#include <fstream>
#include <string>

/**
* \brief A structure for holing values in a FASTQ entry.
*
* Stores the four fields of a FASTQ entry as strings.
*/
struct FastqEntry {
  std::string name;
  std::string seq;
  std::string separator;
  std::string quality;
};


/**
* \brief FASTQ file reader
*
* A class to read and parse single and paired end FASTQ files
*/
class FastqReader {
public:
  /**
  * Opens a single-end FASTQ file. Throws a runtime error if the file
  * cannot be opened.
  *
  * @param [in] in_file FASTQ file name
  */
  FastqReader(const std::string &in_file);
  /**
  * Opens a paired-end FASTQ file. Throws a runtime error if either file
  * cannot be opened.
  *
  * @param [in] in_file_1 first FASTQ file name
  * @param [in] in_file_2 second FASTQ file name
  */
  FastqReader(const std::string &in_file_1, const std::string &in_file_2);
  
  /**
  * Closes any open files
  */
  ~FastqReader();

  /**
  * Read an entry from a single-end FASTQ file and populates
  * a @see FastqEntry.
  * 
  * @param [out] e FastqEntry to populate with the FASTQ entry
  * @return True on successfully reading a FASTQ entry. Flase if
  *   end of file is reached. 
  */
  bool read_se_entry(FastqEntry &e);
  /**
  * Read a pair of entries from paired-end FASTQ files and populates
  * two @see FastqEntry.
  * 
  * @param [out] e1 FastqEntry to populate with the first FASTQ entry
  * @param [out] e2 FastqEntry to populate with the second FASTQ entry
  * @return True on successfully reading a FASTQ entry. Flase if
  *   end of file is reached. 
  */
  bool read_pe_entry(FastqEntry &e1, FastqEntry &e2);


private:
  std::ifstream in_1;
  std::ifstream in_2;

  bool is_pe;
};

std::ostream&
operator<<(std::ostream &out, const FastqEntry &e);

#endif
