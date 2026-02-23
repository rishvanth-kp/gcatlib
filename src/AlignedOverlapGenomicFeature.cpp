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

#include <sstream>
#include <fstream>
#include <cassert>
#include <unordered_set>

#include "GtfReader.hpp"
#include "BedReader.hpp"
#include "AlignedOverlapGenomicFeature.hpp"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::unordered_set;

AlignedOverlapGenomicFeature::AlignedOverlapGenomicFeature() {
  bc_counter = 0;
  feature_counter = 0;
  min_frag_len = 0;
  max_frag_len = std::numeric_limits<size_t>::max();
  shift_left = 4;
  shift_right = -5;
}

AlignedOverlapGenomicFeature::~AlignedOverlapGenomicFeature() {

}


void
AlignedOverlapGenomicFeature::add_gtf_features(const string& gtf_file) {

  GtfReader gtf_reader(gtf_file);
  GencodeGtfEntry entry;

  GenomicStepVector<FeatureVector<string>> gtf_features;
  unordered_set<string> chroms;

  while (gtf_reader.read_gencode_gtf_line(entry)) {
    // keep track of all the chromosomes
    chroms.insert(entry.name);

    // keep track of genes, UTR, CDS. Need to keep genes to figure out
    // introns later.
    if (entry.feature == "gene") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{entry.gene_name});
    }
    else if (entry.feature == "UTR") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{"UTR"});
    }
    else if (entry.feature == "CDS") {
      gtf_features.add(entry.name, entry.start-1, entry.end,
        FeatureVector<string>{"CDS"});
    }
  }
  
  // don't add chrM features
  chroms.erase("chrM");

  for (auto chrom_it = chroms.begin(); chrom_it != chroms.end(); ++chrom_it) {

    vector<pair<GenomicRegion, FeatureVector<string>>> out;
    gtf_features.at(*chrom_it, out);
    for (auto it = out.begin(); it != out.end(); ++it) {
      bool cds = false;
      bool utr = false;
      bool gene = false;

      for (size_t j = 0; j < it->second.size(); ++j) {
        if (it->second.at(j) == "CDS") {
          cds = true;
        }
        else if (it->second.at(j) == "UTR") {
          utr = true;
        }
        else {
          gene = true;
        }
      }

      // build reference
      // A region can be associated with multiple features.
      // If a region covers an UTR or CDS, it is assined as that.
      // If a region covers a gene, but does not have an UTR or CDS,
      // it is assigned as intronic.

      if (cds) {
        genomic_features.add(it->first, FeatureVector<string>{"CDS"});
      }
      if (utr) {
        genomic_features.add(it->first, FeatureVector<string>{"UTR"});
      }
      if (!cds & !utr & gene) {
        genomic_features.add(it->first, FeatureVector<string>{"intron"});
      }

    }
  }
  // keep track of the feature types
  feature_index["CDS"] = feature_counter++;
  feature_index["UTR"] = feature_counter++;
  feature_index["intron"] = feature_counter++;
  feature_index["intergene"] = feature_counter++;

  // add "chrM" as one large feature so that any read aligning to chrM is 
  // assigned to it
  genomic_features.add("chrM", 0, 10000000, FeatureVector<string>{"chrM"});
  feature_index["chrM"] = feature_counter++;

}

void
AlignedOverlapGenomicFeature::add_bed_features(const std::string& bed_file,
                                               const std::string feature_name) {

  BedReader bed_reader(bed_file);
  GenomicRegion bed_region;
  while (bed_reader.read_bed3_line(bed_region)) {
    genomic_features.add(bed_region, FeatureVector<string>{feature_name});
  }

  // added to index
  feature_index[feature_name] = feature_counter++;

}


void 
AlignedOverlapGenomicFeature::set_min_frag_len(const size_t min_len) {
  min_frag_len = min_len;
}

void 
AlignedOverlapGenomicFeature::set_max_frag_len(const size_t max_len) {
  max_frag_len = max_len;
}


void 
AlignedOverlapGenomicFeature::set_shift_left(const int left) {
  shift_left = left;
}

void 
AlignedOverlapGenomicFeature::set_shift_right(const int right) {
  shift_right = right;
}


void
AlignedOverlapGenomicFeature::process_barcodes(const std::string& bc_file) {

  std::ifstream bc_in(bc_file);
  if (!bc_in) {
    throw std::runtime_error("Cannot open " + bc_file);
  }
  string line;

  // keep crate an index for all barcodes
  while (getline(bc_in, line)) {
    string barcode;
    std::istringstream iss(line);
    iss >> barcode;

    bc_index[barcode] = bc_counter++;
  }
  bc_in.close();

  // initialize count matrices
  feature_counts.resize(bc_index.size());
  feature_frag_len.resize(bc_index.size());
  feature_align_len.resize(bc_index.size());
  for (size_t i = 0; i < feature_counts.size(); ++i) {
    feature_counts[i].resize(feature_index.size());
    feature_frag_len[i].resize(feature_index.size());
    feature_align_len[i].resize(feature_index.size());
  }

  // initialize counted_bases and counted_frags vector
  counted_bases.resize(bc_index.size());
  counted_frags.resize(bc_index.size());
}

void
AlignedOverlapGenomicFeature::add(const SamEntry &e1,
                                  const SamEntry &e2,
                                  const string &bc) {

  unordered_map<string, size_t>::iterator bc_it;
  // search for valid barcode
  bc_it = bc_index.find(bc);
  if (bc_it != bc_index.end()) {
    size_t entry_bc_index = bc_it->second;

    // find the start and end location of the reads
    size_t e1_start = e1.pos - 1;
    size_t e1_end = SamCigar::reference_end_pos(e1);
    size_t e2_start = e2.pos - 1;
    size_t e2_end = SamCigar::reference_end_pos(e2);

    // find the stand and end location of the fragments
    size_t frag_start = e1_start;
    if (e2_start < e1_start)
      frag_start = e2_start;

    size_t frag_end = e1_end;
    if (e2_end > e1_end)
      frag_end = e2_end;

    // add the ATAC fragment shift
    frag_start += shift_left;
    frag_end += shift_right;

    size_t frag_len = frag_end - frag_start;

    /*
    if (frag_len != std::abs(e1.tlen)) {
      cerr << e1.qname << "\t" << frag_len << "\t"
           << e1.tlen << "\t" << e1.tlen << endl;
    }
    */

    // count a fragment if it falls withing the fragment length range
    if ((frag_len >= min_frag_len) && (frag_len <= max_frag_len)) {
      // find overlapping regions
      vector<pair<GenomicRegion, FeatureVector<string>>> out;
      genomic_features.at(GenomicRegion(e1.rname, frag_start, frag_end),
                          out, true);

      // count for the start of fragment
      // size if 0 if the alignment is to a chr that is not present in the 
      // features
      if (out.size() > 0) {
        ++counted_frags[entry_bc_index];
        auto first_it = out.begin();
        if (first_it->second.size() == 0) {
          size_t index = feature_index["intergene"];
          ++feature_counts[entry_bc_index][index];
          feature_frag_len[entry_bc_index][index] += frag_len;
        }

        for (size_t j = 0; j < first_it->second.size(); ++j) {
          size_t index = feature_index[first_it->second.at(j)];
          ++feature_counts[entry_bc_index][index];
          feature_frag_len[entry_bc_index][index] += frag_len;
        }
      }

      // count for the entire fragment
      for (auto it = out.begin(); it != out.end(); ++it) {
        const size_t match_len = it->first.end - it->first.start;
        // add match length for the sample.
        counted_bases[entry_bc_index] += match_len;

        // intergenic, if not aligned to any other feature
        if (it->second.size() == 0) {
          size_t index = feature_index["intergene"];
          feature_align_len[entry_bc_index][index] += match_len;
        }

        for (size_t j = 0; j < it->second.size(); ++j) {
          size_t index = feature_index[it->second.at(j)];
          feature_align_len[entry_bc_index][index] += match_len;
        }

      }
    }
    
  }

}

void AlignedOverlapGenomicFeature::feature_counts_to_file(
      const std::string &file_prefix) const {

  std::ofstream out_counts(file_prefix + "_feature_counts.txt");
  std::ofstream out_align_len(file_prefix + "_feature_align_len.txt");
  std::ofstream out_frag_len(file_prefix + "_feature_frag_len.txt");

  // write header
  out_counts << "barcode" << "\tcount";
  out_align_len << "barcode" << "\tbases";
  out_frag_len << "barcode";
  vector<pair<string, size_t>> features;
  for (auto it = feature_index.begin(); it != feature_index.end(); ++it) {
    features.push_back(std::make_pair(it->first, it->second));
    out_counts << "\t" << it->first;
    out_align_len << "\t" << it->first;
    out_frag_len << "\t" << it->first;
  }
  out_counts << endl;
  out_align_len << endl;
  out_frag_len << endl;



  for (auto it = bc_index.begin(); it != bc_index.end(); ++it) {
    // barcode names
    out_counts << it->first;
    out_align_len << it->first;
    out_frag_len << it->first;

    // base or frag count
    out_counts << "\t" << counted_frags[it->second];
    out_align_len << "\t" << counted_bases[it->second];

    // feature counts
    for (size_t j = 0; j < features.size(); ++j) {

      size_t counts = feature_counts[it->second][features[j].second];
      out_counts << "\t" << counts;
      out_align_len << "\t"
                    << feature_align_len[it->second][features[j].second];
      if (counts) {
        float mean_frag_len = static_cast<float>
          (feature_frag_len[it->second][features[j].second]) /
          static_cast<float>(counts);
        out_frag_len << "\t" << mean_frag_len;
      }
      else {
        out_frag_len << "\t0";
      }
    }
    out_counts << endl;
    out_align_len << endl;
    out_frag_len << endl;
  }

  out_counts.close();
  out_align_len.close();
  out_frag_len.close();
}
