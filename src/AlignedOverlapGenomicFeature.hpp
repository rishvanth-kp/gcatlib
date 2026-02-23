/*
* AlignedOverlapGenomicFeature: Determine feature of aligned reads
*   allowing a read to align to multiple overlapping features
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

#ifndef ALIGNED_OVERLAP_GENOMIC_FEATURE_HPP
#define ALIGNED_OVERLAP_GENOMIC_FEATURE_HPP

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include "SamEntry.hpp"
#include "FeatureVector.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStepVector.hpp"


/**
* \brief Assign genomic region(s) to aligned reads.
*
* A class to assign reads to genomic regions processed from GTF 
* and BED files. The aligned reads can be assigned to multiple
* overlapping regions. 
*/
class AlignedOverlapGenomicFeature {
public:
  AlignedOverlapGenomicFeature();
  ~AlignedOverlapGenomicFeature();

  void add_gtf_features(const std::string& gtf_file);
  void add_bed_features(const std::string& bed_file,
                        const std::string feature_name);

  void set_min_frag_len(const size_t min_len);
  void set_max_frag_len(const size_t max_len);
  
  void set_shift_left(const int left);
  void set_shift_right(const int right);

  void process_barcodes(const std::string& bc_file);

  void add(const SamEntry &e1, const SamEntry &e2, const string &bc);

  void feature_counts_to_file(const std::string &file_prefix) const;

private:
  GenomicStepVector<FeatureVector<std::string>> genomic_features;

  std::unordered_map<string, size_t> feature_index;
  size_t feature_counter;

  std::unordered_map<string, size_t> bc_index;
  size_t bc_counter;

  size_t min_frag_len;
  size_t max_frag_len;
    
  int shift_left;
  int shift_right;

  std::vector<size_t> counted_bases;
  std::vector<size_t> counted_frags;

  vector<vector<size_t>> feature_counts;
  vector<vector<size_t>> feature_align_len;
  vector<vector<size_t>> feature_frag_len;
};

#endif
