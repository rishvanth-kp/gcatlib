/*
* Metagene: Determine the normalized gene locations for metagene plots
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
#include <vector>

#include "GenomicRegion.hpp"
#include "FeatureVector.hpp"
#include "GenomicStepVector.hpp"


/**
* \brief Map genomic positions to normalized positions along features,
* for metagene plots.
*
* Metagene reads a set of features (e.g. transcripts) from a BED file
* and divides the spliced length of each feature into a fixed number
* of equal divisions, numbered 0 to divisions - 1 from the 5' end to
* the 3' end. For a feature on the minus strand, division 0 is
* therefore at its highest genomic coordinate. Queries return, for each
* genomic position, which features cover it and in which division, so
* that signal from features of different lengths can be combined on a
* common scale.
*
* For example, with 4 divisions, a plus-strand feature with two 64-bp
* regions, chr1:100-164 and chr1:200-264 (128 bp spliced), is divided
* into 32-bp parts:
* \code
* chr1:100-132 -> 0    chr1:200-232 -> 2
* chr1:132-164 -> 1    chr1:232-264 -> 3
* \endcode
* On the minus strand, the same regions would be numbered 3, 2, 1, 0.
*
* The BED file must meet these requirements:
* - At least 6 columns (BED6): the name (column 4) identifies the
*   feature, and the strand (column 6) is '-' for the minus strand;
*   any other value, including '.', is treated as plus.
* - All regions of a feature (e.g. the exons of a transcript) share
*   its name and are on consecutive lines. If the same name appears
*   again later, those regions are treated as a separate feature.
* - The regions of a feature are listed in ascending genomic order,
*   even on the minus strand, and do not overlap each other.
*
* Different features may overlap; a position covered by several
* features is reported for each of them.
*
* @note Division boundaries can be off by one base. The division of
*   each base is computed by repeatedly adding divisions / length to a
*   floating-point value, so rounding error can place a boundary one
*   base late. Features shorter than the number of divisions skip
*   divisions, because each base advances by more than one division.
*/
class Metagene {

public:
  /**
  * Reads the features from a BED file and computes their divisions.
  *
  * @param [in] bed_file  BED6 file of features; see the class
  *   documentation for its requirements.
  * @param [in] divisions Number of divisions per feature. Defaults to
  *   100.
  * @throws std::runtime_error if the file cannot be opened.
  * @warning A file with fewer than 6 columns (e.g. BED3), or an empty
  *   file, causes undefined behaviour (typically a crash).
  */
  Metagene(const std::string &bed_file, const size_t divisions = 100);
  // Metagene(const std::string& bed_file) : Metagene(bed_file, 100) {};
  /**
  * Default destructor.
  */
  ~Metagene();


  /**
  * Returns the features and divisions covering the genomic region in.
  * Each output pair is a sub-region of in, together with the list of
  * (feature name, division) pairs covering it; adjacent bases with the
  * same list are merged into one sub-region. Bases not covered by any
  * feature are omitted.
  *
  * @param [in]  in  Genomic region to query.
  * @param [out] out (sub-region, list of (feature name, division))
  *   pairs, in order of position. Cleared before populating, and left
  *   empty if the chromosome of in has no features.
  * @warning in.start must be less than in.end. As with
  *   GenomicStepVector::at(), a query with start >= end on a chromosome
  *   that has features reads past the end of an internal vector
  *   (undefined behaviour, typically a crash).
  */
  void at(const GenomicRegion &in,
          vector<pair<GenomicRegion,
                      FeatureVector<pair<string, size_t>>>> &out) const;
  /**
  * For each feature overlapping the genomic region in, returns the
  * lowest (most 5') and highest (most 3') division covered by in. For
  * example, for a read this gives the span of each feature, in
  * divisions, that the read covers. The three output vectors are
  * parallel, with one entry per feature, in order of first appearance
  * along in. Features with the same name are combined, even if they
  * were read as separate features (see get_feature_names()).
  *
  * @param [in]  in      Genomic region to query.
  * @param [out] feature Names of the overlapping features. Cleared
  *   before populating.
  * @param [out] first   Lowest division of each feature covered by in.
  *   Cleared before populating.
  * @param [out] last    Highest division of each feature covered by in.
  *   Cleared before populating.
  * @warning As with at(), in.start must be less than in.end.
  */
  void at_ends(const GenomicRegion &in,
               vector<string> &feature,
               vector<size_t> &first, vector<size_t> &last) const;

  /**
  * Returns the number of features read from the BED file. A name that
  * appears in more than one group of consecutive lines is counted once
  * per group.
  *
  * @return Number of features.
  */
  size_t get_n_features() const;
  /**
  * Returns the names of the features, in the order they appear in the
  * BED file. A name that appears in more than one group of consecutive
  * lines is listed once per group.
  *
  * @param [out] names Feature names. Replaced with the new list.
  */
  void get_feature_names(vector<string> &names) const;

private:
  GenomicStepVector<FeatureVector<pair<string, size_t>>> metagene;

  size_t n_features;
  size_t n_divisions;
  vector<string> feature_names;

  struct FeatureRegions {
    string chrom;
    size_t start;
    size_t end;
    string name;
    bool dir; // 1: pos, 0:neg
  };

  void process_feature(const vector<FeatureRegions> &feature);
  void add_features(const string &bed_file);
  void add_region(const GenomicRegion& g, const vector<string> &fields,
           vector<FeatureRegions> &feature);


};

# endif
