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

#include <unordered_map>
#include <fstream>
#include <cmath>

#include "Metagene.hpp"
#include "BedReader.hpp"

using std::string;
using std::vector;
using std::unordered_map;


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

  // keep track of feature name
  feature_names.push_back(feature[0].name);

  size_t feature_size = 0;
  // cout << "Number of regions: " << feature.size() << endl;
  for (size_t i = 0; i < feature.size(); ++i) {
    // cout << feature[i].chrom << "\t"
    //      << feature[i].start << "\t" << feature[i].end << "\t"
    //      << feature[i].name << "\t" << feature[i].dir << endl;

    feature_size += (feature[i].end - feature[i].start);
  }

  // cout << "feature size: " << feature_size << endl;
  double base_pct = static_cast<double>(n_divisions) /
                    static_cast<double>(feature_size);
  // cout << "base_size: " << base_pct << endl;

  double base_counter = 0;



  for (size_t i = 0; i < feature.size(); ++i) {
    size_t prev_counter = static_cast<size_t>(floor(base_counter));
    size_t prev_start = feature[i].start;
    for (size_t j = feature[i].start; j < feature[i].end; ++j) {

      /*
      metagene.add(feature[i].chrom, j, j + 1,
        FeatureVector<pair<string, size_t>>
        (make_pair(feature[i].name, static_cast<size_t>(floor(base_counter)) )) );
      cout << j << "\t" << base_counter << "\t"
           << static_cast<size_t>(floor(base_counter)) << endl;
      */
      // cout << j << "\t" << base_counter << "\t"
      //      << static_cast<size_t>(floor(base_counter)) << endl;

      size_t curr_counter = static_cast<size_t>(floor(base_counter));
      if (curr_counter != prev_counter) {
        if (feature[i].dir) {
          metagene.add(feature[i].chrom, prev_start, j,
            FeatureVector<pair<string, size_t>>
            (make_pair(feature[i].name, prev_counter )) );
          // cout << "Adding: " << prev_start << " " << j << " "
          //      << prev_counter << endl;
        }
        else {
          metagene.add(feature[i].chrom, prev_start, j,
            FeatureVector<pair<string, size_t>>
            (make_pair(feature[i].name, n_divisions - prev_counter - 1)) );
          // cout << "Adding: " << prev_start << " " << j << " "
          //      << n_divisions - prev_counter - 1 << endl;

        }

        prev_counter = curr_counter;
        prev_start = j;
      }

      base_counter += base_pct;

/*
      size_t curr_counter = static_cast<size_t>(floor(base_counter));
      if (curr_counter != prev_counter) {
        metagene.add(feature[i].chrom, prev_start, j + 1,
          FeatureVector<pair<string, size_t>>
          (make_pair(feature[i].name, prev_counter )) );
        cout << "Adding: " << prev_start << " " << j + 1 << " "
             << prev_counter << endl;
        prev_counter = curr_counter;
        prev_start = j + 1;
      }
*/

    }
    /*
    metagene.add(feature[i].chrom, feature[i].end, feature[i].end + 1,
      FeatureVector<pair<string, size_t>>
      (make_pair(feature[i].name, static_cast<size_t>(floor(base_counter)) )) );
    cout << feature[i].end << "\t" << base_counter << "\t"
         << static_cast<size_t>(floor(base_counter)) << endl;
    */
    // cout << feature[i].end << "\t" << base_counter << "\t"
    //      << static_cast<size_t>(floor(base_counter)) << endl;

    if (feature[i].dir) {
      metagene.add(feature[i].chrom, prev_start, feature[i].end,
        FeatureVector<pair<string, size_t>>
        (make_pair(feature[i].name, prev_counter )) );
      // cout << "Adding e: " << prev_start << " " << feature[i].end
      //      << " " << prev_counter << endl;
    }
    else {
      metagene.add(feature[i].chrom, prev_start, feature[i].end,
        FeatureVector<pair<string, size_t>>
        (make_pair(feature[i].name, n_divisions - prev_counter - 1 )) );
      // cout << "Adding e: " << prev_start << " " << feature[i].end
      //      << " " << n_divisions - prev_counter - 1 << endl;

    }
  }

}

void
Metagene::add_features(const string& bed_file) {
  BedReader reader(bed_file);
  GenomicRegion g;
  vector<string> fields;

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

  n_features = 0;
  n_divisions = divisions;
  add_features(bed_file);
}

Metagene::~Metagene() {

}

size_t
Metagene::get_n_features() const {
  return n_features;
}


void
Metagene::get_feature_names(vector<string> &names) const {
  names.clear();
  names = feature_names;
}


void
Metagene::at(const GenomicRegion &in,
             vector<pair<GenomicRegion,
                         FeatureVector<pair<string, size_t>>>> &out) const {
  metagene.at(in, out);
}

void
Metagene::at_ends(const GenomicRegion &in,
                  vector<string> &feature,
                  vector<size_t> &first, vector<size_t> &last) const {

  feature.clear();
  first.clear();
  last.clear();

  unordered_map<string, size_t> regions;
  size_t regions_counter = 0;

  // cout  << "querying " << in << endl;

  vector<pair<GenomicRegion, FeatureVector<pair<string, size_t>>>>
    feature_out;
  metagene.at(in, feature_out);

  for (size_t i = 0; i < feature_out.size(); ++i) {
    // cout << feature_out[i].first << "\t";
    for (size_t j = 0; j < feature_out[i].second.size(); ++j) {
      // cout << feature_out[i].second.at(j).first << "\t";
      // cout << feature_out[i].second.at(j).second << "\t";

      unordered_map<string, size_t>::const_iterator it;
      it = regions.find(feature_out[i].second.at(j).first);
      // check if it is a new region
      if (it == regions.end()) {
        // keep track of region
        regions[feature_out[i].second.at(j).first] = regions_counter++;
        // add the region to output
        feature.push_back(feature_out[i].second.at(j).first);
        first.push_back(feature_out[i].second.at(j).second);
        last.push_back(feature_out[i].second.at(j).second);
      }
      // if the region already exists just upate the min and max
      else {
        // update the min and the max
        // (only one will be updated depending on the direction).
        if (feature_out[i].second.at(j).second < first[it->second]) {
          first[it->second] = feature_out[i].second.at(j).second;
        }
        if (feature_out[i].second.at(j).second > last[it->second]) {
          last[it->second] = feature_out[i].second.at(j).second;
        }

      }
    }
    // cout << endl;
  }
}
