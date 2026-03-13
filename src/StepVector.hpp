/*
* StepVector:
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

#ifndef STEP_VECTOR_HPP
#define STEP_VECTOR_HPP

#include <iostream>
#include <vector>
#include <map>

using std::map;
using std::cout;
using std::cerr;
using std::endl;
using std::pair;
using std::vector;

/**
* \brief A memory-efficient step function over integer positions.
*
* StepVector stores a piecewise-constant function over a 1D integer
* coordinate space. Rather than storing a value at every position,
* it stores only the positions where the value changes (step boundaries),
* making it efficient for sparse data.
*
* Values are added over half-open intervals [start, end) and are
* cumulative: adding to an interval that already has a value accumulates
* the new value on top of the existing one. Positions that have never
* been written to return the default-constructed value T{}.
*
* The type T must support the += and + operators and must be
* default-constructible.
*
* \tparam T The value type stored at each step.
*/
template<typename T>
class StepVector {
public:
  /**
  * Default constructor.
  */
  StepVector();

  /**
  * Accumulates val over the half-open interval [start, end).
  * If the interval overlaps with previously added intervals, the
  * values are summed. If start >= end, the call is a no-op.
  *
  * @param [in] start Start of the interval (inclusive, 0-based).
  * @param [in] end   End of the interval (exclusive, 0-based).
  * @param [in] val   Value to accumulate over the interval.
  */
  void add(const size_t start, const size_t end, const T &val);

  /**
  * Returns the step boundaries and their values for the half-open
  * interval [start, end).
  *
  * The output is a vector of (position, value) pairs representing
  * the step function within the queried range. Consecutive entries
  * out[i] and out[i+1] define the sub-interval [out[i].first,
  * out[i+1].first) with constant value out[i].second. The first
  * entry is always (start, value_at_start) and the last entry is
  * always (end, value_at_end), bounding the output. Any step
  * boundaries strictly inside (start, end) are included between them.
  *
  * If start >= end, out is cleared and left empty.
  *
  * @param [in]  start Start of the query interval (inclusive, 0-based).
  * @param [in]  end   End of the query interval (exclusive, 0-based).
  * @param [out] out   Vector of (position, value) pairs representing
  *   the step boundaries within [start, end). Cleared before populating.
  */
  void at_range(const size_t start, const size_t end,
                vector<pair<size_t, T>>& out) const;

  /**
  * Returns the value at a single position.
  * Returns the default-constructed value T{} if no value has been
  * accumulated at or before pos.
  *
  * @param [in] pos Query position (0-based).
  * @return The accumulated value at pos.
  */
  T at(const size_t pos) const;

  /**
  * Prints all internal step boundaries to stdout as tab-separated
  * position-value pairs, one per line. Intended for debugging.
  */
  void print_elements();


private:

  map<size_t, T> step_vec;
  bool VERBOSE;
};

template<typename T>
StepVector<T>::StepVector() {
}





template<typename T>
void
StepVector<T>::add(const size_t start,
                const size_t end,
                const T &val) {

  if (start < end) {

    // Inserting at the end
    // if the end postion is an an existing element, there is nothing to
    // be done
    typename map<size_t, T>::iterator end_it = step_vec.lower_bound(end);
    // ending before or after the end of all existing entries,
    // add default value.
    if (end_it == step_vec.begin() || end_it == step_vec.end()) {
      // cout << "End next element at end" << endl;
      T default_val{};
      end_it = step_vec.insert(pair<size_t, T>(end, default_val)).first;
    }
    // ending before the end of exisiting entries,
    else if (!(end_it->first == end)) {
      T prev_val = (--end_it)->second;
      // cout << "insering end prev val: " << prev_val << endl;
      end_it = step_vec.insert(pair<size_t, T>(end, prev_val)).first;
    }
    // else {
    //   cout << "ending at an exitsting element" << endl;
    // }

    // Inserting at the start
    typename map<size_t, T>::iterator start_it = step_vec.lower_bound(start);
    if (start_it->first == start) {
      // If the element exists, update it's value.
      // cout << "start exists" << endl;
      start_it->second += val;

      // Iterator points to the start element. Increment to point to
      // next element
      ++start_it;
    }
    else {
      // If an element does not exist, insert a new element with
      // the current value added to the previous element's value.
      // cout << "start does not exist" << endl;
      if (start_it == step_vec.begin()) {
        start_it =
          step_vec.insert(pair<size_t, T>(start, val)).first;
      }
      else {
        T prev_val = (--start_it)->second;
        start_it =
          step_vec.insert(pair<size_t, T>(start, prev_val + val)).first;
      }

      // Iterator points to the inserted element
      // Increment to point to next element
     ++start_it;
    }

    // Modifying values in [start + 1. end)
    // cout << "inserting in the middle" << endl;
    // cout << start_it->first << "\t" << end_it->first << endl;
    for (auto it = start_it; it != end_it; ++it) {
      // T prev_val = it->second;
      // cout << "prev val: " << prev_val << endl;
      it->second += val;
    }

  }
}


template<typename T>
T
StepVector<T>::at(const size_t pos) const {
  typename map<size_t, T>::const_iterator it = step_vec.upper_bound(pos);
  if (it == step_vec.begin())
    return T{};
  else
    return (--it)->second;

}

template<typename T>
void
StepVector<T>::at_range(const size_t start, const size_t end,
                        vector<pair<size_t, T>>& out) const {

  out.clear();

  if (start < end) {

    typename map<size_t, T>::const_iterator start_it, end_it;
    start_it = step_vec.upper_bound(start);
    if (start_it == step_vec.begin()) {
      out.push_back(std::make_pair(start, T{}));
    }
    else {
      out.push_back(std::make_pair(start, (--start_it)->second));
      ++start_it;
    }

    end_it = step_vec.upper_bound(end);
    for (auto it = start_it; it != end_it; ++it) {
      out.push_back(std::make_pair(it->first, it->second));
    }

    if (end_it == step_vec.begin()) {
      out.push_back(std::make_pair(end, T{}));
    }
    else {
      --end_it;
      if (end_it->first != end)
        out.push_back(std::make_pair(end, end_it->second));
    }

  }
}

template<typename T>
void
StepVector<T>::print_elements() {

  for (auto it = step_vec.begin(); it != step_vec.end(); ++it) {
    cout << it->first << "\t" << it->second << endl;
  }
  cout << endl;
}


#endif
