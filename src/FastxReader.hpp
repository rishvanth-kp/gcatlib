/*
* FastxReader: class to read fastq/a files
* Copyright (C) 2025 Rishvanth Prabakar
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

#ifndef FASTX_READER_HPP
#define FASTX_READER_HPP

#include <iostream>
#include <string>
#include <unordered_map>

#include "htslib/hts.h"
#include "htslib/sam.h"


/**
* \brief FASTQ/FASTA file reader.
* 
* Used htslib to read FASTX entries. This class reads the entries
* and parses them into a FastxEntry. 
*
*/
class FastxReader {
public:
  /**
  * Open a fastx file, verifies that it is in the right format,
  * and initialized the htslib handlers.
  * 
  * @param [in] in_file FASTQ/FASTA (compressed) file name. 
  */
  FastxReader(const std::string &in_file);
  /**
  * Closes the file and destroys htslib handlers.
  */
  ~FastxReader();

private:
  htsFile *hts;
  sam_hdr_t *header;
  bam1_t *fastx_entry;
};


#endif
