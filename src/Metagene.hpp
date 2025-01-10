/*
* Metagene: Determine the normazlied gene locations for metagene plots 
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

#ifndef METAGENE_HPP
#define METAGENE_HPP

#include <iostream>
#include <string>

#include "GenomicRegion.hpp"
#include "FeatureVector.hpp"
#include "GenomicStepVector.hpp"


class Metagene {

public:
  Metagene(const std::string& bed_file, const size_t divisions = 100);
  // Metagene(const std::string& bed_file) : Metagene(bed_file, 100) {};
  ~Metagene();


  void at(const GenomicRegion g) const;

private:
  GenomicStepVector<FeatureVector<pair<string, size_t>>> metagene;
  
  size_t n_features;
  size_t n_divisions;
  
  
  struct FeatureRegions {
    string chrom;
    size_t start;
    size_t end;
    string name;
    bool dir; // 1: pos, 0:neg 
  }; 

  void process_feature(const vector<FeatureRegions>& feature);
  void add_features(const string& bed_file);
  void add_region(const GenomicRegion& g, const vector<string>& fields, 
           vector<FeatureRegions>& feature);


};

# endif
