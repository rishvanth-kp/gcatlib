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

#include <fstream>
#include <cmath>

#include "Metagene.hpp"
#include "BedReader.hpp"

using std::string;
using std::vector;



void
Metagene::add_region(const GenomicRegion& g, const vector<string>& fields, 
           vector<FeatureRegions>& feature) {

  feature.push_back(FeatureRegions()); 
  feature.back().chrom = g.name; 
  feature.back().start = g.start;
  feature.back().end = g.end;
  feature.back().name = fields[0];
  if (fields[2][0] == '-')
    feature.back().dir = false;
  else
    feature.back().dir = true;

}

void
Metagene::process_feature(const vector<FeatureRegions>& feature) {

  size_t feature_size = 0;
  cout << "Number of regions: " << feature.size() << endl;
  for (size_t i = 0; i < feature.size(); ++i) {
    cout << feature[i].chrom << "\t" 
         << feature[i].start << "\t" << feature[i].end << "\t"
         << feature[i].name << "\t" << feature[i].dir << endl;

    feature_size += (feature[i].end - feature[i].start);
  }
  
  cout << "feature size: " << feature_size << endl;
  float base_pct = static_cast<float>(n_divisions) / 
                    static_cast<float>(feature_size);
  cout << "base_size: " << base_pct << endl;

  float base_counter = 0;
  if (!feature[0].dir)
    base_counter = 100;
    

  for (size_t i = 0; i < feature.size(); ++i) {
    for (size_t j = feature[i].start; j < feature[i].end; ++j) {
      metagene.add(feature[i].chrom, j, j + 1,
        FeatureVector<pair<string, size_t>>
        (make_pair(feature[i].name, static_cast<size_t>(round(base_counter)) )) );
      cout << j << "\t" << base_counter << "\t" << static_cast<size_t>(round(base_counter)) << endl;
      if (feature[i].dir) {
        base_counter += base_pct;
      }
      else {
        base_counter -= base_pct;
      }
    }    
    metagene.add(feature[i].chrom, feature[i].end, feature[i].end + 1,
      FeatureVector<pair<string, size_t>>
      (make_pair(feature[i].name, static_cast<size_t>(round(base_counter)) )) );
    cout << feature[i].end << "\t" << base_counter << "\t" << static_cast<size_t>(round(base_counter)) << endl;
  } 

}

void 
Metagene::add_features(const string& bed_file) {
  BedReader reader(bed_file);
  GenomicRegion g;
  vector<string> fields;
  cout << bed_file << endl;
 
  vector<FeatureRegions> feature;
  reader.read_bed_line(g, fields);
  add_region(g, fields, feature);

  while(reader.read_bed_line(g, fields)) {
    
    // check if the regions belongs to a new feature.
    // if it does, process the previous feature
    if (feature.back().name != fields[0]) {
      process_feature(feature);
      ++n_features;    
      // get rid of the old feature
      feature.clear();
    }
    
    // if the region is part of an old feature, just keep track
    add_region(g, fields, feature);
  }

  // process the very last feature
  process_feature(feature);
  ++n_features;
}


Metagene::Metagene(const string& bed_file, const size_t divisions) {
  
  cout << "Metagene construct" << endl; 
  n_features = 0;
  n_divisions = divisions;
  add_features(bed_file);
}

Metagene::~Metagene() {

}



void 
Metagene::at(const GenomicRegion g) const {

}
