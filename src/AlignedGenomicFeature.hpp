/*
* AlignedGenomicFeature: Determine feature of aligned reads
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

#ifndef ALIGNED_GENOMIC_FEATURE_HPP
#define ALIGNED_GENOMIC_FEATURE_HPP

#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "SamEntry.hpp"
#include "GtfReader.hpp"
#include "FeatureVector.hpp"
#include "GenomicRegion.hpp"
#include "GenomicStepVector.hpp"

/**
* \brief Count aligned bases per genomic feature category.
*
* Builds a map of feature categories from a GENCODE GTF file and, for
* each added alignment, counts how many aligned bases fall into each
* category. Counts can be kept for a bulk sample, or per cell barcode
* for single-cell data, optionally with UMI deduplication. Only bases
* aligned by M, = or X CIGAR operations are counted; insertions,
* deletions, clipped bases, and skipped regions (N) are not.
*
* Each position of the genome is assigned exactly one category:
* - chrM_\<gene_type\>: any gene on chrM, named by its gene_type
*   (e.g. chrM_protein_coding). Overlapping genes on chrM are not
*   reported as multi_gene.
* - multi_gene: covered by more than one gene.
* - For a position covered by a single protein-coding gene:
*   - pc_cds: in a CDS but not a UTR.
*   - pc_utr: in a UTR but not a CDS.
*   - pc_both: in both a CDS and a UTR (e.g. from different
*     transcripts).
*   - pc_none: in an exon that is neither CDS nor UTR.
*   - pc_intron: in the gene but not in any exon.
* - \<gene_type\>: covered by a single gene of any other type (e.g.
*   lncRNA). Exons and introns are not distinguished for these genes.
* - intergenic: on a chromosome in the GTF, but not in any gene.
*
* Only categories that occur in the GTF file, plus intergenic, are
* reported.
*
* Methods must be called in this order:
* \code
* AlignedGenomicFeature agf;
* agf.preprocess_gff("gencode.gtf");        // 1. build the feature map
* agf.process_barcodes("barcodes.txt");     // 2. single-cell only
* SamReader reader("aln.bam");
* SamEntry e;
* string cb, umi;
* while (reader.read_sam_line(e)) {         // 3. add alignments
*   if (SamTags::get_tag(e.tags, "CB", cb) &&
*       SamTags::get_tag(e.tags, "UB", umi))
*     agf.add(e, cb, umi);
* }
* agf.feature_count_to_file("counts.txt");  // 4. write the counts
* \endcode
*
* @note The GTF must follow GENCODE conventions: gene types are read
*   from the gene_type attribute (Ensembl's gene_biotype is not
*   recognised, so with an Ensembl GTF every base would be counted as
*   intergenic), and the mitochondrial chromosome must be named chrM.
*   No filtering is done: strand is not considered; secondary,
*   supplementary, duplicate, and low-MAPQ alignments are counted; and
*   overlapping mates of a pair are both counted. Filter alignments
*   before calling add() as needed.
*/
class AlignedGenomicFeature {
public:
  /**
  * Constructs an empty counter. Call preprocess_gff() before adding
  * alignments.
  */
  AlignedGenomicFeature();
  /**
  * Default destructor.
  */
  ~AlignedGenomicFeature();

  /**
  * Reads a GENCODE GTF file and builds the map of feature categories
  * described in the class documentation, from its gene, exon, CDS, and
  * UTR lines. GTF coordinates are converted to 0-based, half-open.
  *
  * @param [in] gff_file GTF file name.
  * @throws std::runtime_error if the file cannot be opened.
  * @warning Call this only once per object. Calling it again adds the
  *   new features on top of the existing ones, and overlapping category
  *   names are concatenated (e.g. "pc_cdspc_cds"), corrupting the
  *   counts.
  */
  void preprocess_gff(const std::string& gff_file);

  /**
  * Reads a list of cell barcodes and sets up per-barcode counts. Only
  * alignments whose barcode is in this list are counted by the
  * barcoded add() overloads. The file has one barcode per line; only
  * the first whitespace-separated token of each line is used, so extra
  * columns are ignored.
  *
  * After this is called, feature_count_to_file() writes the per-barcode
  * counts instead of the bulk counts.
  *
  * @param [in] bc_file Barcode file name.
  * @throws std::runtime_error if the file cannot be opened.
  * @warning Must be called after preprocess_gff(), because the count
  *   tables are sized by the number of feature categories at this
  *   point; calling it first leads to out-of-bounds writes in add().
  *   The file must not contain duplicate barcodes: a repeated barcode
  *   gets an index beyond the end of the count tables, which also leads
  *   to out-of-bounds writes. A blank line adds an empty barcode, so
  *   more than one blank line counts as a duplicate.
  */
  void process_barcodes(const std::string& bc_file);

  /**
  * Adds an alignment to the bulk counts: each aligned base is counted
  * once under the category of its position. Bases on chromosomes that
  * are not in the GTF file are ignored entirely; they are not counted
  * as intergenic and are not included in the total. Unmapped reads
  * have no aligned bases and add nothing.
  *
  * @param [in] e SAM entry to add.
  */
  void add(const SamEntry &e);

  /**
  * Adds an alignment to the counts of cell barcode bc, in the same way
  * as add(e). Alignments whose barcode is not in the list from
  * process_barcodes() are silently ignored. Barcoded counts are kept
  * separately from the bulk counts of add(e).
  *
  * @param [in] e  SAM entry to add.
  * @param [in] bc Cell barcode of the alignment (e.g. the CB tag).
  */
  void add(const SamEntry &e, const std::string &bc);

  /**
  * Adds an alignment to the counts of cell barcode bc, counting each
  * molecule only once. A molecule is identified by its (barcode, UMI)
  * pair: the first alignment seen for a pair is counted as in
  * add(e, bc), and later alignments with the same pair are skipped,
  * wherever they align. This is deliberate: each molecule contributes
  * its bases once, regardless of which feature categories its reads
  * fall in.
  *
  * As a result, the counted alignment for a molecule is the first one
  * in input order, even if it contributes no bases (e.g. an unmapped
  * read), and two different molecules in the same cell that share a
  * UMI by chance are counted as one. The set of seen pairs is kept for
  * the lifetime of the object and is not reset by clear_counts().
  *
  * @param [in] e   SAM entry to add.
  * @param [in] bc  Cell barcode of the alignment (e.g. the CB tag).
  * @param [in] umi UMI of the alignment (e.g. the UB tag).
  * @note The barcode and UMI are joined without a separator to form the
  *   key, so pairs such as ("AC", "GT") and ("A", "CGT") collide. This
  *   only matters if barcode or UMI lengths vary.
  */
  void add(const SamEntry &e, const std::string &bc,
           const std::string &umi);

  /**
  * Resets the bulk counts and total from add(e) to zero. The feature
  * map built by preprocess_gff() is kept.
  *
  * @warning Only the bulk counts are reset. The per-barcode counts and
  *   the set of seen (barcode, UMI) pairs are not cleared.
  */
  void clear_counts();

  /**
  * Returns the bulk counts from add(e): the number of aligned bases in
  * each feature category, sorted by category name. Categories with zero
  * bases are included. Per-barcode counts are not available through
  * this method; use feature_count_to_file().
  *
  * @param [out] counts  (category, bases) pairs. Cleared before
  *   populating.
  * @param [out] n_bases Total number of aligned bases counted.
  */
  void get_feature_counts(std::vector<pair<std::string, size_t>> &counts,
                          size_t &n_bases);
  /**
  * Writes the counts to a tab-separated text file.
  *
  * Without barcodes: one line per category, sorted by category name,
  * with no header. Each line has the category, its bases, and its
  * percentage of all counted bases. The percentage is nan if no bases
  * were counted.
  *
  * After process_barcodes(): a header line (barcode, bases, then one
  * column per category), followed by one line per barcode in the list
  * with the barcode, its total counted bases, and its bases in each
  * category. Barcodes with no counts are included. Bulk counts from
  * add(e) are not written in this mode. The order of barcode rows and
  * category columns is arbitrary and does not follow the barcode file,
  * but the columns of every row match the header.
  *
  * @param [in] file Output file name. An existing file is overwritten.
  *   Failure to open the file is not reported.
  */
  void feature_count_to_file(const std::string& file) const;


private:
  GenomicStepVector<std::string> genomic_features;
  
  // for bulk samples that do not have barcodes
  // This is slightly ineffiecnet since these structures
  // are created even the sample is barcoded. 
  std::unordered_map<std::string, size_t> feature_count;
  size_t match_bases;

  // for barcoded samples
  std::unordered_map<std::string, size_t> feature_index;
  std::unordered_map<std::string, size_t> bc_index;
  size_t bc_counter;
  bool is_barcoded;

  std::vector<std::vector<size_t>> bc_feature_count;
  std::vector<size_t> bc_counted_bases;

  // for barcoded and UMI samples
  std::unordered_set<std::string> seen_cb_umi;
};


#endif
