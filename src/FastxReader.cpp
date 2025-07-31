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

#include "FastxReader.hpp"
#include <sstream>

using std::cout;
using std::cerr;
using std::endl;
using std::string;

FastxReader::FastxReader(const string &in_file) {

  // open the file
  if (!(hts = hts_open(in_file.c_str(), "r")))
    throw std::runtime_error("Cannot open FASTX file: " + in_file);
 
  // check format
  if (hts_get_format(hts)->format == fasta_format)
    is_fasta = true;
  else if (hts_get_format(hts)->format == fastq_format)
    is_fasta = false;
  else
    throw std::runtime_error("Invalid FASTX file format: " + in_file);

  // need to read header even for a fastx file
  if (!(header = sam_hdr_read(hts)))
    throw std::runtime_error("No header in FASTX file!: " + in_file);

  // initialize data structure for reading 
  if (!(fastx_data = bam_init1())) 
    throw std::runtime_error("Failed to initialize fastx data");  

  // initialize the kstring
  ks_initialize(&fastx_kstr);  

}

FastxReader::~FastxReader() {
  hts_close(hts);
  sam_hdr_destroy(header);
  bam_destroy1(fastx_data);
  ks_free(&fastx_kstr);
}

void 
FastxReader::set_2nd_name_column() {
  hts_set_opt(hts, FASTQ_OPT_NAME2); 
}

bool
FastxReader::read_fastx_entry(FastxEntry &e) {
  
  if (sam_read1(hts, header, fastx_data) >= 0) {
    if (sam_format1(header, fastx_data, &fastx_kstr) >= 0) {
      std::istringstream iss(fastx_kstr.s);
      string dummy;
      iss >> e.name >> dummy >> dummy >> dummy >> dummy
          >> dummy >> dummy >> dummy >> dummy >> e.seq >> e.qual;
      cout << e.name << endl;
      cout << e.seq << endl;
      cout << e.qual << endl;
      return true;
    }
  }
  return false;
}
