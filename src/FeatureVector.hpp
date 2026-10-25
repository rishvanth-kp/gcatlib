/*
* FeatureVector: Vector compatible with GenomicStepVector
* Copyright (C) 2021 Rishvanth Prabakar
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

#ifndef FEATURE_VECTOR_HPP
#define FEATURE_VECTOR_HPP

#include <iostream>
#include <vector>

using std::vector;

/**
* \brief An ordered list of values that accumulates by concatenation.
*
* FeatureVector is designed to be used as the value type T of a
* \ref StepVector or \ref GenomicStepVector when each position should
* record which features cover it, rather than a numeric sum. Its + and
* += operators concatenate lists instead of adding numbers, so where
* intervals overlap, the overlapping region holds the features from all
* of them. A default-constructed (empty) FeatureVector plays the role of
* zero: positions not covered by any feature hold an empty list.
*
* \code
* GenomicStepVector<FeatureVector<string>> features;
* features.add("chr1", 100, 200, FeatureVector<string>{"geneA"});
* features.add("chr1", 150, 300, FeatureVector<string>{"geneB"});
* // [100, 150) -> {geneA}
* // [150, 200) -> {geneA, geneB}
* // [200, 300) -> {geneB}
* \endcode
*
* Elements are kept in the order they were added and are not
* deduplicated: adding the same value twice over a region stores it
* twice. Equality is order-sensitive, so {A, B} and {B, A} compare
* unequal, and \ref GenomicStepVector will not merge adjacent regions
* holding them.
*
* \tparam T The element type. Must be copyable and support ==.
*/
template<typename T>
class FeatureVector {
public:
  /**
  * Constructs an empty FeatureVector.
  */
  FeatureVector() = default;
  /**
  * Constructs a FeatureVector holding a single element.
  *
  * @param [in] in The element to store.
  */
  explicit FeatureVector(const T &in);
  /**
  * Copy constructor.
  *
  * @param [in] in FeatureVector to copy.
  */
  FeatureVector(const FeatureVector<T> &in);

  /**
  * Returns a new FeatureVector containing the elements of this vector
  * followed by the elements of in. Neither operand is modified.
  *
  * @param [in] in FeatureVector whose elements are appended.
  * @return The concatenated FeatureVector.
  */
  FeatureVector<T> operator+(const FeatureVector<T> &in);
  /**
  * Appends all elements of in to the end of this vector.
  *
  * @param [in] in FeatureVector whose elements are appended.
  */
  void operator+=(const FeatureVector<T> in);
  /**
  * Returns true if both vectors have the same size and equal elements
  * in the same order.
  *
  * @param [in] in FeatureVector to compare against.
  * @return True if the vectors are equal.
  */
  bool operator==(const FeatureVector &in) const;
  /**
  * Returns the negation of operator==.
  *
  * @param [in] in FeatureVector to compare against.
  * @return True if the vectors differ.
  */
  bool operator!=(const FeatureVector &in) const;

  /**
  * Appends an element to the end of the vector.
  *
  * @param [in] in Element to append.
  */
  void push_back(const T in);
  /**
  * Returns a copy of the element at index i.
  *
  * @param [in] i Index of the element (0-based).
  * @return The element at index i.
  * @warning No bounds checking is performed, unlike std::vector::at.
  *   i must be less than size().
  */
  T at(const size_t i) const;
  /**
  * Returns the number of elements.
  *
  * @return Number of elements in the vector.
  */
  size_t size() const;

private:
  vector<T> features;
};



template<typename T>
FeatureVector<T>::FeatureVector(const T &in) {
  features.push_back(in);
}


template<typename T>
FeatureVector<T>::FeatureVector(const FeatureVector<T> &in) {
  features.clear();
  for (size_t i = 0; i < in.size(); ++i) {
    features.push_back(in.at(i));
  }
}


template<typename T>
FeatureVector<T>
FeatureVector<T>::operator+(const FeatureVector<T> &in) {
  FeatureVector<T> out;
  for (size_t i = 0; i < features.size(); ++i) {
    out.push_back(features[i]);
  }
  for (size_t i = 0; i < in.size(); ++i) {
    out.push_back(in.at(i));
  }
  return out;
}


template<typename T>
void
FeatureVector<T>::operator+=(const FeatureVector<T> in) {
  for (size_t i = 0; i < in.size(); ++i) {
    features.push_back(in.at(i));
  }
}


template<typename T>
bool
FeatureVector<T>::operator==(const FeatureVector &in) const {
  bool match = false;
  if (in.size() == features.size()) {
    size_t i = 0;
    while ((i < features.size()) && (in.at(i) == features[i]))
      ++i;
    if (i == features.size())
      match = true;
  }
  return match;
}


template<typename T>
bool
FeatureVector<T>::operator!=(const FeatureVector<T> &in) const {
  return !(*this == in);
}


template<typename T>
void
FeatureVector<T>::push_back(const T in) {
  features.push_back(in);
}

template<typename T>
T
FeatureVector<T>::at(const size_t i) const {
  return features[i];
}

template<typename T>
size_t
FeatureVector<T>::size() const {
  return features.size();
}

#endif
