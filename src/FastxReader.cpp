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

using std::cout;
using std::cerr;
using std::endl;
using std::string;

FastxReader::FastxReader(const string &in_file) {

  // open the file
  if (!(hts = hts_open(in_file.c_str(), "r")))
    throw std::runtime_error("Cannot open FASTX file: " + in_file);
 
  // check format

  // need to read header even for a fastx file

  // initialize data structure for reading 
   
}

FastxReader::~FastxReader() {
  hts_close(hts);
}
