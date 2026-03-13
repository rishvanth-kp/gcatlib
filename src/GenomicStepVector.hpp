/*
* GenomcRegion: class to store genomic intervals
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

#ifndef GENOMIC_STEP_VECTOR_HPP
#define GENOMIC_STEP_VECTOR_HPP

#include <iostream>
#include <string>
#include <unordered_map>

#include "StepVector.hpp"
#include "GenomicRegion.hpp"

using std::pair;
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::make_pair;
using std::unordered_map;

/**
* \brief A chromosome-aware step function for genomic data.
*
* GenomicStepVector extends \ref StepVector to genomic coordinates by
* maintaining an independent \ref StepVector per chromosome. Chromosomes
* are created on demand the first time an interval is added. The
* accumulation semantics are the same as \ref StepVector: adding to an
* interval that already has a value sums the new value on top of it.
*
* The main query method at(GenomicRegion, out) returns results as a
* vector of (GenomicRegion, T) pairs where consecutive step boundaries
* with equal values are merged into a single region, and zero-valued
* intervals can optionally be filtered out.
*
* The type T must support the += and + operators, equality comparison,
* and must be default-constructible.
*
* \tparam T The value type stored at each step.
*/
template<typename T>
class GenomicStepVector {
public:
  /**
  * Default constructor. Initializes an empty genomic step vector
  * with no chromosomes.
  */
  GenomicStepVector();

  /**
  * Accumulates val over the half-open interval [start, end) on the
  * given chromosome. If the chromosome does not yet exist, it is
  * created. Accumulation semantics are inherited from \ref StepVector:
  * overlapping intervals sum their values.
  *
  * @param [in] chr   Chromosome name.
  * @param [in] start Start of the interval (inclusive, 0-based).
  * @param [in] end   End of the interval (exclusive, 0-based).
  * @param [in] val   Value to accumulate over the interval.
  */
  void add(const string chr, const size_t start, const size_t end, const T &val);

  /**
  * Convenience overload of add() using a \ref GenomicRegion.
  * Delegates to add(chr, start, end, val) using g.name, g.start, g.end.
  *
  * @param [in] g   GenomicRegion defining the chromosome and interval.
  * @param [in] val Value to accumulate over the interval.
  */
  void add(const GenomicRegion &g, const T &val);

  /**
  * Prints the step boundaries in [start, end) on the given chromosome
  * to stdout as tab-separated position-value pairs. If the chromosome
  * does not exist, prints a message to stderr. Intended for debugging.
  *
  * @param [in] chr   Chromosome name.
  * @param [in] start Start of the query interval (inclusive, 0-based).
  * @param [in] end   End of the query interval (exclusive, 0-based).
  */
  void at(const string chr, const size_t start, const size_t end) const;

  /**
  * Returns the accumulated values within the genomic region g as a
  * vector of (GenomicRegion, T) pairs. Each pair represents a contiguous
  * sub-interval where the value is constant. Adjacent step boundaries
  * with the same value are merged into a single entry.
  *
  * By default (keep_0 = false), intervals where the value equals T{}
  * (the default-constructed zero value) are excluded from the output.
  * Set keep_0 = true to include them.
  *
  * If the chromosome in g does not exist in the vector, out is cleared
  * and left empty.
  *
  * @param [in]  g      GenomicRegion defining the chromosome and interval
  *   to query.
  * @param [out] out    Vector of (GenomicRegion, T) pairs representing
  *   the merged constant-value sub-intervals within g. Cleared before
  *   populating.
  * @param [in]  keep_0 If false (default), intervals with value T{}
  *   are excluded. If true, all intervals are returned.
  */
  void at(const GenomicRegion &g,
          vector<pair<GenomicRegion, T>> &out,
          bool keep_0 = false) const;

  /**
  * Queries the entire chromosome and returns all accumulated intervals
  * as a vector of (GenomicRegion, T) pairs. Equivalent to calling
  * at(GenomicRegion{chr, 0, SIZE_MAX}, out). Zero-valued intervals are
  * excluded (keep_0 = false).
  *
  * @param [in]  chr Chromosome name to query.
  * @param [out] out Vector of (GenomicRegion, T) pairs for the entire
  *   chromosome. Cleared before populating.
  */
  void at(const string chr, vector<pair<GenomicRegion, T>> &out);

  /**
  * Returns the number of distinct chromosomes that have been added.
  *
  * @return Number of chromosomes.
  */
  size_t chrom_count() const { return n_chrom; }

  /**
  * Returns the total number of add() calls made across all chromosomes.
  * This counts calls, not unique genomic intervals.
  *
  * @return Total number of add() calls.
  */
  size_t entry_count() const { return n_entry; }

private:
  size_t n_chrom{};
  size_t n_entry{};

  typedef unordered_map<string, StepVector<T>> GenomicStepVectorType;
  GenomicStepVectorType genomic_vector;
};


template<typename T>
GenomicStepVector<T>::GenomicStepVector() {
}

template<typename T>
void
GenomicStepVector<T>::add(const string chr, const size_t start,
                          const size_t end, const T &val) {

  typename GenomicStepVectorType::iterator chr_map;
  chr_map = genomic_vector.find(chr);
  if (chr_map == genomic_vector.end()) {
    genomic_vector[chr].add(start, end, val);
    ++n_chrom;
    ++n_entry;
    // genomic_vector[chr].print_elements();
  }
  else {
    chr_map->second.add(start, end, val);
    ++n_entry;
    // chr_map->second.print_elements();
  }
}

template<typename T>
void
GenomicStepVector<T>::add(const GenomicRegion &g, const T &val) {
  add(g.name, g.start, g.end, val);
}

template<typename T>
void
GenomicStepVector<T>::at(const string chr, const size_t start,
                         const size_t end) const {

  typename GenomicStepVectorType::const_iterator chr_map;
  chr_map = genomic_vector.find(chr);
  if (chr_map != genomic_vector.end()) {
    vector<pair<size_t, T>> out;
    chr_map->second.at_range(start, end, out);

    for (auto it = out.begin(); it != out.end(); ++it) {
      cout << it->first << "\t" << it->second << endl;
    }
  }
  else {
    cerr << "chromosome does not exist" << endl;
  }
}

template<typename T>
void
GenomicStepVector<T>::at(const GenomicRegion &g,
                         vector<pair<GenomicRegion, T>> &out,
                         bool keep_0) const {

  out.clear();

  typename GenomicStepVectorType::const_iterator chr_map;
  chr_map = genomic_vector.find(g.name);
  if (chr_map != genomic_vector.end()) {
    vector<pair<size_t, T>> step_out;
    chr_map->second.at_range(g.start, g.end, step_out);

    GenomicRegion region;
    size_t i = 0;
    while(i < step_out.size() - 1) {
      size_t j = i + 1;
      if (keep_0 || (!keep_0 && step_out[i].second != T{})) {
        while (j < step_out.size() - 1 &&
               step_out[i].second == step_out[j].second) {
          ++j;
        }
        region.name = g.name;
        region.start = step_out[i].first;
        region.end = step_out[j].first;
        out.push_back(make_pair(region, step_out[i].second));
      }
      i = j;
    }
  }
}


template<typename T>
void 
GenomicStepVector<T>::at(const string chr, 
                         vector<pair<GenomicRegion, T>> &out) {
  at(GenomicRegion{chr, 0, std::numeric_limits<size_t>::max()}, out);
}

#endif
