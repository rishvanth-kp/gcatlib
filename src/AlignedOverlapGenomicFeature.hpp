/*
* AlignedOverlapGenomicFeature: Determine features of aligned fragments
*   allowing a fragment to overlap multiple features
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
* \brief Assign aligned fragments to overlapping genomic features.
*
* A class to assign paired-end fragments (e.g. from ATAC-seq) to
* genomic features loaded from GTF and BED files, and count them per
* cell barcode. Features may overlap, so a fragment can be assigned to
* several features at once.
*
* Features from add_gtf_features():
* - CDS and UTR: regions covered by a CDS or UTR line. A region covered
*   by both is assigned to both.
* - intron: in a gene but in no CDS or UTR. This includes the exons of
*   non-coding genes (e.g. lncRNA) and the non-coding exons of
*   protein-coding genes; exon lines are not used.
* - intergene: not in any feature.
* - chrM: all of chrM, as a single feature.
*
* add_bed_features() adds each region of a BED file as a feature with
* a given name (e.g. peaks).
*
* For each read pair, the fragment runs from the leftmost start to the
* rightmost end of the two alignments. Its ends are then shifted (by
* default +4 at the left end and -5 at the right end, the usual Tn5
* offset for ATAC-seq), and the fragment is kept only if its shifted
* length is within the limits set by set_min_frag_len() and
* set_max_frag_len(). For each kept fragment:
* - It is counted once under each feature at its shifted left end.
*   The right end is not considered.
* - Every base of the fragment interval, including the unsequenced gap
*   between the mates, is added to the features it overlaps.
* - Its length is added to each feature it was counted under, to give
*   the mean fragment length per feature.
*
* All features must be loaded before process_barcodes(), and
* process_barcodes() must be called before add():
* \code
* AlignedOverlapGenomicFeature aogf;
* aogf.add_gtf_features("gencode.gtf");        // 1. features, GTF first
* aogf.add_bed_features("peaks.bed", "peaks");
* aogf.set_min_frag_len(1);                    //    optional settings
* aogf.set_max_frag_len(1000);
* aogf.process_barcodes("barcodes.txt");       // 2. barcodes
* SamReader reader("aln.bam");
* SamEntry e1, e2;
* string cb;
* while (reader.read_pe_sam(e1, e2)) {         // 3. add read pairs
*   if (SamTags::get_tag(e1.tags, "CB", cb))
*     aogf.add(e1, e2, cb);
* }
* aogf.feature_counts_to_file("sample");       // 4. write the counts
* \endcode
*
* @note The GTF must follow GENCODE conventions: UTRs must be UTR lines
*   (Ensembl's five_prime_utr and three_prime_utr lines are not
*   recognised, so those regions would be labelled intron), and the
*   mitochondrial chromosome must be named chrM. Counting is per
*   barcode only; there is no bulk mode, so use a single barcode for
*   bulk data. No filtering is done: strand is not considered, and
*   secondary, supplementary, duplicate, and low-MAPQ alignments are
*   counted. Both reads of a pair must be mapped to the same chromosome.
*/
class AlignedOverlapGenomicFeature {
public:
  /**
  * Constructs an empty object with the default settings: no fragment
  * length limits (0 to SIZE_MAX), and shifts of +4 (left) and -5
  * (right).
  */
  AlignedOverlapGenomicFeature();
  /**
  * Default destructor.
  */
  ~AlignedOverlapGenomicFeature();

  /**
  * Reads a GENCODE GTF file and adds the CDS, UTR, intron, intergene,
  * and chrM features described in the class documentation. GTF
  * coordinates are converted to 0-based, half-open. Features on chrM in
  * the GTF are not used; instead, chrM is added as a single feature
  * covering positions 0 to 10,000,000.
  *
  * @param [in] gtf_file GTF file name.
  * @throws std::runtime_error if the file cannot be opened.
  * @warning This must be called even if only BED features are wanted,
  *   because it defines the intergene feature. Without it, fragments
  *   that are in no feature are counted under the first BED feature,
  *   and the intergene column in the output repeats that feature's
  *   counts. Call it only once, before process_barcodes(); calling it
  *   again leads to out-of-bounds writes in add().
  */
  void add_gtf_features(const std::string& gtf_file);
  /**
  * Reads a BED file and adds each region as a feature named
  * feature_name. Only the first three columns are used (see
  * BedReader::read_bed3_line()).
  *
  * @param [in] bed_file     BED file name.
  * @param [in] feature_name Name of the feature, used as its column
  *   name in the output files.
  * @throws std::runtime_error if the file cannot be opened.
  * @warning feature_name must be unique: each call must use a new name
  *   that is also not one of the GTF feature names (CDS, UTR, intron,
  *   intergene, chrM). Reusing a name, or calling this after
  *   process_barcodes(), leads to out-of-bounds writes in add().
  *   Regions in the file should not overlap each other: where they do,
  *   a fragment is counted once per overlapping region. Merge
  *   overlapping regions first if needed.
  */
  void add_bed_features(const std::string& bed_file,
                        const std::string feature_name);

  /**
  * Sets the minimum fragment length, after shifting. Shorter fragments
  * are not counted. The default is 0.
  *
  * @param [in] min_len Minimum fragment length (inclusive).
  */
  void set_min_frag_len(const size_t min_len);
  /**
  * Sets the maximum fragment length, after shifting. Longer fragments
  * are not counted. The default is SIZE_MAX (no limit).
  *
  * @param [in] max_len Maximum fragment length (inclusive).
  * @note Setting a finite maximum, together with set_min_frag_len(1),
  *   also prevents the crash on very short fragments described in
  *   add().
  */
  void set_max_frag_len(const size_t max_len);

  /**
  * Sets the shift applied to the left (start) end of each fragment.
  * The default is +4, the usual Tn5 offset for ATAC-seq. Use 0 for no
  * shift.
  *
  * @param [in] left Number of bases to add to the fragment start (may
  *   be negative).
  */
  void set_shift_left(const int left);
  /**
  * Sets the shift applied to the right (end) end of each fragment.
  * The default is -5, the usual Tn5 offset for ATAC-seq. Use 0 for no
  * shift.
  *
  * @param [in] right Number of bases to add to the fragment end (may
  *   be negative).
  */
  void set_shift_right(const int right);

  /**
  * Reads a list of cell barcodes and sets up the per-barcode counts.
  * Only read pairs whose barcode is in this list are counted. The file
  * has one barcode per line; only the first whitespace-separated token
  * of each line is used.
  *
  * @param [in] bc_file Barcode file name.
  * @throws std::runtime_error if the file cannot be opened.
  * @warning Must be called after all add_gtf_features() and
  *   add_bed_features() calls, because the count tables are sized by
  *   the number of features at this point. The file must not contain
  *   duplicate barcodes: a repeated barcode gets an index beyond the end
  *   of the count tables, leading to out-of-bounds writes in add(). A
  *   blank line adds an empty barcode, so more than one blank line
  *   counts as a duplicate.
  */
  void process_barcodes(const std::string& bc_file);

  /**
  * Adds a read pair (one fragment) to the counts of cell barcode bc, as
  * described in the class documentation. The pair is ignored if bc is
  * not in the list from process_barcodes(), if the shifted fragment
  * length is outside the set limits, or if its chromosome has no
  * features at all. The chromosome is taken from e1; e2 is assumed to
  * be on the same chromosome.
  *
  * @param [in] e1 First read of the pair.
  * @param [in] e2 Second read of the pair (the mate of e1).
  * @param [in] bc Cell barcode of the pair (e.g. the CB tag).
  * @warning If shifting leaves the fragment with start >= end (with the
  *   default shifts, any fragment of 9 bp or less) and the length
  *   limits do not exclude it, the call reads past the end of an
  *   internal vector (undefined behaviour, typically a crash). The
  *   default limits never exclude such fragments, because their length
  *   wraps around to a very large value. Use set_min_frag_len(1) and a
  *   finite set_max_frag_len() to exclude them.
  */
  void add(const SamEntry &e1, const SamEntry &e2, const string &bc);

  /**
  * Writes the per-barcode counts to three tab-separated files, each
  * with a header line and one line per barcode in the list:
  * - `<prefix>_feature_counts.txt`: barcode, total number of counted
  *   fragments, then the number of fragments counted under each
  *   feature.
  * - `<prefix>_feature_align_len.txt`: barcode, total number of
  *   fragment bases, then the number of fragment bases overlapping each
  *   feature.
  * - `<prefix>_feature_frag_len.txt`: barcode, then the mean length of
  *   the fragments counted under each feature (0 if there are none).
  *
  * Because features can overlap, the per-feature columns of a row can
  * sum to more than its total. Barcodes with no counts are included.
  * The order of barcode rows and feature columns is arbitrary and does
  * not follow the input files, but it is the same in all three files.
  *
  * @param [in] file_prefix Prefix of the output file names. Existing
  *   files are overwritten, and failure to open a file is not reported.
  */
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
