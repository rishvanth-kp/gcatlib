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
* \brief A structure for holding values in a FASTQ entry.
*
* Stores the four lines of a FASTQ entry as strings, exactly as read
* from the file, without the trailing newline. Only '\\n' is stripped,
* so files with Windows (CRLF) line endings leave a trailing '\\r' on
* each field.
*/
struct FastqEntry {
  std::string name;      /**< Header line, including the leading '@'. */
  std::string seq;       /**< Sequence line. */
  std::string separator; /**< Separator line, including the leading '+'. */
  std::string quality;   /**< Quality line. */
};


/**
* \brief FASTQ file reader.
*
* A class to read single-end and paired-end FASTQ files. Each record is
* read as exactly four lines into a \ref FastqEntry.
*
* @note Only uncompressed FASTQ is supported; gzipped files
*   (.fastq.gz) are not decompressed. Records are not validated: the
*   reader does not check for the leading '@' or '+', or that the
*   sequence and quality have the same length, and multi-line FASTQ
*   records are not supported.
*/
class FastqReader {
public:
  /**
  * Opens a single-end FASTQ file.
  *
  * @param [in] in_file FASTQ file name.
  * @throws std::runtime_error if the file cannot be opened.
  */
  FastqReader(const std::string &in_file);
  /**
  * Opens a pair of paired-end FASTQ files (e.g. R1 and R2).
  *
  * @param [in] in_file_1 First FASTQ file name.
  * @param [in] in_file_2 Second FASTQ file name.
  * @throws std::runtime_error if either file cannot be opened.
  */
  FastqReader(const std::string &in_file_1, const std::string &in_file_2);

  /**
  * Closes any open files.
  */
  ~FastqReader();

  /**
  * Reads the next entry from a single-end FASTQ file and populates
  * a \ref FastqEntry. If the reader was opened with two files, only
  * the first file is read.
  *
  * @param [out] e FastqEntry to populate with the FASTQ entry.
  *   Contents are undefined if false is returned.
  * @return True on successfully reading a FASTQ entry. False if
  *   end of file is reached.
  */
  bool read_se_entry(FastqEntry &e);
  /**
  * Reads the next pair of entries from paired-end FASTQ files and
  * populates two \ref FastqEntry objects. The files are read in
  * lockstep, so e1 and e2 are the n-th records of the first and second
  * files; read names are not checked to match. Requires the two-file
  * constructor; on a single-end reader this always returns false.
  *
  * @param [out] e1 FastqEntry to populate with the entry from the first
  *   file. Contents are undefined if false is returned.
  * @param [out] e2 FastqEntry to populate with the entry from the
  *   second file. Contents are undefined if false is returned.
  * @return True on successfully reading a pair of entries. False if
  *   the end of either file is reached.
  * @warning If the two files have different numbers of records, the
  *   extra records in the longer file are silently ignored. This is
  *   indistinguishable from a normal end of file.
  */
  bool read_pe_entry(FastqEntry &e1, FastqEntry &e2);


private:
  std::ifstream in_1;
  std::ifstream in_2;

  bool is_pe;
};

/**
* Writes a FastqEntry to an output stream as four newline-terminated
* lines (name, seq, separator, quality). Because the fields hold the
* raw lines, an entry read by \ref FastqReader is written back unchanged.
*
* @param [in,out] out Output stream to write to.
* @param [in]     e   FastqEntry to write.
* @return Reference to the output stream.
*/
std::ostream&
operator<<(std::ostream &out, const FastqEntry &e);

#endif
