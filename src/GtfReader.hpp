/*
* GtfReader: class to read GTF files
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

#ifndef GTF_READER_HPP
#define GTF_READER_HPP

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <unordered_map>

using std::string;
using std::vector;
using std::unordered_map;

/**
* \brief A structure for holding values in a GTF entry.
*
* GTF is a 9-column, tab-delimited format. The first 8 columns are
* mandatory; the 9th holds a semicolon-separated list of attributes,
* where each tag and value are separated by a space.
*
* Coordinates are stored as in the file: 1-based, with an inclusive
* end. To convert to a \ref GenomicRegion (0-based, half-open), use
* start - 1 as the start and end unchanged.
*/
struct GtfEntry {
  string name;    /**< Sequence name (column 1), e.g. chr1. */
  string source;  /**< Source (column 2): the program or database that
                       generated the feature. */
  string feature; /**< Feature type (column 3), e.g. gene, transcript,
                       exon, CDS. */
  size_t start;   /**< Start position (column 4), 1-based. */
  size_t end;     /**< End position (column 5), 1-based and inclusive. */
  double score;   /**< Score (column 6). If the column is '.', this is set
                       to std::numeric_limits<float>::min() (about 1.18e-38,
                       the smallest positive float, not a negative value). */
  char strand;    /**< Strand (column 7): '+', '-', or '.'. */
  size_t frame;   /**< Frame (column 8). Set to SIZE_MAX if the column is
                       '.'. Numeric frames (0, 1, 2) are not parsed, so for
                       CDS lines this field is left unset. */
  unordered_map<string, string> attribute; /**< Attributes (column 9),
                       mapping each tag to its value. Values are stored
                       exactly as in the file, including surrounding
                       double quotes (e.g. gene_id is stored as "ENSG1",
                       with the quotes). If a tag appears more than once
                       (e.g. GENCODE's tag), only the last value is
                       kept. */
};


/**
* \brief A GTF entry with the common GENCODE attributes extracted.
*
* Holds the same columns as \ref GtfEntry, except score, which is not
* kept. Commonly used GENCODE attributes are stored in named fields
* instead of a map. It is populated by
* GtfReader::read_gencode_gtf_line() and
* GtfReader::read_gencode_gtf_file().
*
* The extraction assumes GENCODE formatting:
* - String attributes are assumed to be in double quotes, and their
*   first and last characters are removed. An unquoted value would
*   lose real characters.
* - exon_number is assumed to be unquoted, as in GENCODE. A quoted
*   value (as in Ensembl GTFs) is read as 0.
* - A missing attribute gives an empty string, or 0 for exon_number.
*
* Coordinates are 1-based with an inclusive end, as in \ref GtfEntry.
*/
struct GencodeGtfEntry {
  string name;    /**< Sequence name (column 1). */
  string source;  /**< Source (column 2). */
  string feature; /**< Feature type (column 3). */
  size_t start;   /**< Start position (column 4), 1-based. */
  size_t end;     /**< End position (column 5), 1-based and inclusive. */
  char strand;    /**< Strand (column 7): '+', '-', or '.'. */
  size_t frame;   /**< Frame (column 8). SIZE_MAX if the column is '.'.
                       Numeric frames are not parsed, so this is
                       uninitialized for CDS lines. */

  //attributes
  string gene_id;         /**< gene_id attribute, without quotes. */
  string transcript_id;   /**< transcript_id attribute, without quotes. */
  string gene_type;       /**< gene_type attribute, without quotes. */
  string gene_name;       /**< gene_name attribute, without quotes. */
  string transcript_type; /**< transcript_type attribute, without quotes. */
  string transcript_name; /**< transcript_name attribute, without quotes. */
  size_t exon_number;     /**< exon_number attribute. 0 if missing or
                               quoted. */
  string exon_id;         /**< exon_id attribute, without quotes. */
};

/**
* \brief GTF reader.
*
* A class for reading and parsing GTF files, line by line or a whole
* file at a time, into \ref GtfEntry or \ref GencodeGtfEntry records.
* All read methods share the same file handle, so each call continues
* from where the previous one stopped.
*
* @note Comment lines starting with '#' are skipped only at the start
*   of the file. Comment, blank, or malformed lines later in the file
*   are not skipped or reported: they are returned as records with
*   bogus values.
*/
class GtfReader {
public:
  /**
  * Opens a GTF file and skips any leading comment lines that start
  * with '#'.
  *
  * @param [in] in_file GTF file name.
  * @throws std::runtime_error if the file cannot be opened.
  */
  GtfReader(const string &in_file);
  /**
  * Closes the GTF file.
  */
  ~GtfReader();

  /**
  * Reads the next line from the GTF file, parses it, and populates a
  * \ref GtfEntry.
  *
  * @param [out] g GtfEntry to populate.
  * @return True on successfully reading a GTF entry. False if end of
  *   file is reached.
  * @warning g.attribute is not cleared before parsing. If the same
  *   GtfEntry is reused across calls, attributes from earlier lines
  *   remain unless the current line overwrites them. Clear g.attribute,
  *   or use a fresh GtfEntry, before each call.
  */
  bool read_gtf_line(GtfEntry &g);
  /**
  * Reads the rest of the GTF file, parses each line, and appends a
  * \ref GtfEntry for each line to g. Existing contents of g are kept.
  * Because all read methods share the same file handle, if this is
  * called after read_gtf_line(), it reads from the next line to the
  * end of the file.
  *
  * @param [in,out] g Vector of GtfEntry to append to.
  */
  void read_gtf_file(vector<GtfEntry> &g);

  /**
  * Reads the next line from the GTF file, parses it, and populates a
  * \ref GencodeGtfEntry with the GENCODE attributes extracted. All
  * fields of g are overwritten on each call.
  *
  * @param [out] g GencodeGtfEntry to populate.
  * @return True on successfully reading a GTF entry. False if end of
  *   file is reached.
  */
  bool read_gencode_gtf_line(GencodeGtfEntry &g);
  /**
  * Reads the rest of the GTF file and appends a \ref GencodeGtfEntry
  * for each line to g. Existing contents of g are kept. Like
  * read_gtf_file(), it starts from the current position in the file.
  *
  * @param [in,out] g Vector of GencodeGtfEntry to append to.
  */
  void read_gencode_gtf_file(vector<GencodeGtfEntry> &g);

private:
  std::ifstream in;

  void parse_gtf_line(const string &in, GtfEntry &g);

};


#endif
