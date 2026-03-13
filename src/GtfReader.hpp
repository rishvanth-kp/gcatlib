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

/*
// GTF format
// GTF uses a 9 column tab delimated format.
// The first 8 columans are mandatroy
// 1: sequence name. E.g. Chr1.
// 2: source: porgram/source that generated the frature
// 3: feature: E.g. gene, transcript
// 4: start: 1-based start position of feature
// 5: end: inclusive end position of feature
// 6: score:
// 7: strand: strand of feature. '+', '-', or '.'
// 8: frame:
// 9: attributes: semicolon seperated list of tag-value paris
//      tag and value are separated with a space
*/
struct GtfEntry {
  string name;
  string source;
  string feature;
  size_t start;
  size_t end;
  double score;
  char strand;
  size_t frame;
  unordered_map<string, string> attribute;
};


struct GencodeGtfEntry {
  string name;
  string source;
  string feature;
  size_t start;
  size_t end;
  char strand;
  size_t frame;

  //attributes
  string gene_id;
  string transcript_id;
  string gene_type;
  string gene_name;
  string transcript_type;
  string transcript_name;
  size_t exon_number;
  string exon_id;
};

/**
* \brief GTF reader.
*
* A class for reading and parsing GTF files.
*
*/
class GtfReader {
public:
  /**
  * Opens a GTF file. Throws a runtime error if the file cannot
  * be opened.
  *
  * \param [in] in_file GTF file name
  */
  GtfReader(const string &in_file);
  /**
  * Closes the GTF file.
  */
  ~GtfReader();

  /**
  * Reads a line from a GTF file, parses it, and populates a @see GtfEntry.
` * 
  * @params [out] g GtfEntry to populate.
  * @return true on successfully reading a gtf entry, and false if end of 
  *   file is reached. 
  */
  bool read_gtf_line(GtfEntry &g);
  /**
  * Reads a GTF file, parses each line, and populates a vector of 
  * @see GtfEntry, Note that read_gtf_file and read_gtf_line uses
  * the same file handler, and so if this is called after calls to
  * read_gtf_line, then it will read from the next line to the end
  * of the file. 
  * 
  * @param [out] g a vector of @see GtfEntry to populate. 
  */
  void read_gtf_file(vector<GtfEntry> &g);

  bool read_gencode_gtf_line(GencodeGtfEntry &g);
  void read_gencode_gtf_file(vector<GencodeGtfEntry> &g);

private:
  std::ifstream in;

  void parse_gtf_line(const string &in, GtfEntry &g);

};


#endif
