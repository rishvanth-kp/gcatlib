/*
 * test_gtf_reader.cpp: Catch2 tests for GTF column 8 (frame) parsing
 */

#include "catch_amalgamated.hpp"

#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <unistd.h>

#include "GtfReader.hpp"

static std::string
write_temp_gtf(const std::string &body) {
  char path[] = "/tmp/gcatlib_gtf_XXXXXX";
  int fd = mkstemp(path);
  REQUIRE(fd != -1);
  FILE *f = fdopen(fd, "w");
  REQUIRE(f != nullptr);
  REQUIRE(fwrite(body.data(), 1, body.size(), f) == body.size());
  REQUIRE(fclose(f) == 0);
  return std::string(path);
}

TEST_CASE("parse numeric and missing GTF frames") {
  const char *gtf =
      "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tgene_id \"g1\";\n"
      "chr1\tsrc\tCDS\t4\t6\t.\t+\t1\tgene_id \"g1\";\n"
      "chr1\tsrc\tCDS\t7\t9\t.\t+\t2\tgene_id \"g1\";\n"
      "chr1\tsrc\texon\t1\t9\t.\t+\t.\tgene_id \"g1\";\n";

  std::string path = write_temp_gtf(gtf);
  GtfReader reader(path);

  GtfEntry a, b, c, d;
  REQUIRE(reader.read_gtf_line(a));
  REQUIRE(reader.read_gtf_line(b));
  REQUIRE(reader.read_gtf_line(c));
  REQUIRE(reader.read_gtf_line(d));
  REQUIRE_FALSE(reader.read_gtf_line(a));

  REQUIRE(a.frame == 0);
  REQUIRE(b.frame == 1);
  REQUIRE(c.frame == 2);
  REQUIRE(d.frame == std::numeric_limits<size_t>::max());

  // Re-read a missing frame after a numeric one must not leak the old value.
  // The entries above already cover that: d follows c (frame 2).

  std::remove(path.c_str());
}

TEST_CASE("GENCODE reader copies parsed numeric frames") {
  const char *gtf =
      "chr1\tHAVANA\tCDS\t10\t12\t.\t-\t1\t"
      "gene_id \"ENSG1\"; transcript_id \"ENST1\"; gene_type \"protein_coding\"; "
      "gene_name \"X\"; transcript_type \"protein_coding\"; transcript_name \"X-1\"; "
      "exon_number 1; exon_id \"ENSE1\";\n";

  std::string path = write_temp_gtf(gtf);
  GtfReader reader(path);

  GencodeGtfEntry g;
  REQUIRE(reader.read_gencode_gtf_line(g));
  REQUIRE(g.frame == 1);
  REQUIRE(g.strand == '-');

  std::remove(path.c_str());
}
