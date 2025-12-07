/*
* PreprocessGff: parse features from GFF file 
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

#ifndef PREPROCESS_GFF_HPP
#define PREPROCESS_GFF_HPP

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdint>
#include <unordered_map>

#include "GtfReader.hpp"

using std::string;
using std::vector;
using std::unordered_map;

class PreprocessGff {
public:
  PreprocessGff(const string &chrom_size_file, bool verbose=false);
  ~PreprocessGff();

  void parse_genome_features(const string &gff_file);


private:
  unordered_map<string, size_t> chrom_sizes;
  bool VERBOSE; 

  void read_chrom_sizes(const string &chrom_size_file);
};

#endif
