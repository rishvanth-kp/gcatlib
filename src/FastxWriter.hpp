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

#ifndef FASTX_WRITER_HPP
#define FASTX_WRITER_HPP

#include <string>

#include "htslib/hts.h"
#include "htslib/sam.h"

#include "FastxEntry.hpp" 

class FastxWriter {
public:

  FastxWriter(const std::string &out_file);

  ~FastxWriter();

  bool write_fastx_entry(FastxEntry &e);

  

private:
  htsFile *hts;
  bam1_t *fastx_data;
};  

#endif
