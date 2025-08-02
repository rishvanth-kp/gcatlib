/*
* FastxWriter: class to write fastq/a files
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

#include <iostream>

#include "FastxWriter.hpp"
  
using std::cout;
using std::cerr;
using std::endl;
using std::string;

FastxWriter::FastxWriter(const std::string &out_file) {

  // determine the write mode from file format
  char mode[4] = "w";
  if (sam_open_mode(mode + 1, out_file.c_str(), NULL) < 0)
    throw std::runtime_error("Invalid FASTX file format: " + out_file);
  
  cout << mode << endl;

  // open the file
  if (!(hts = hts_open(out_file.c_str(), mode)))
    throw std::runtime_error("Cannot open FASTX file: " + out_file);

  // check format
  if ((hts_get_format(hts)->format != fasta_format) &&
      (hts_get_format(hts)->format != fastq_format)) 
    throw std::runtime_error("Invalid FASTX file format: " + out_file); 

  // initialize data structure for writing
  if (!(fastx_data = bam_init1())) 
    throw std::runtime_error("Failed to initialize fastx data for writing");
}

FastxWriter::~FastxWriter() {
  hts_close(hts);
  bam_destroy1(fastx_data);
}

  
bool 
FastxWriter::write_fastx_entry(FastxEntry &e) {

  // need to subtract 33 as sam_format1 adds 33 
  for (size_t i = 0; i < e.qual.length(); ++i)
    e.qual[i] -= 33;

  if (bam_set1(fastx_data, e.name.length(), e.name.c_str(), 
                BAM_FUNMAP, -1, -1, 0, 0, NULL, -1, -1, 
                0, e.seq.length(), e.seq.c_str(), e.qual.c_str(), 0) < 0) {
    return false;
  }
  if (sam_write1(hts, NULL, fastx_data) < 0) {
    return false;
  }

  return true;
}
