# Correctness tests (Catch2 v3)
check_PROGRAMS = test_step_vector test_genomic_step_vector
TESTS = test_step_vector test_genomic_step_vector

test_step_vector_SOURCES = test_step_vector.cpp \
                           catch2/catch_amalgamated.cpp
test_step_vector_CXXFLAGS = -Wall -std=c++14 -g \
                             -I$(top_srcdir)/src \
                             -I$(srcdir)/catch2
test_step_vector_LDADD = $(top_builddir)/src/libgcat.a

test_genomic_step_vector_SOURCES = test_genomic_step_vector.cpp \
                                   catch2/catch_amalgamated.cpp
test_genomic_step_vector_CXXFLAGS = -Wall -std=c++14 -g \
                                     -I$(top_srcdir)/src \
                                     -I$(srcdir)/catch2
test_genomic_step_vector_LDADD = $(top_builddir)/src/libgcat.a

# Performance benchmarks (nanobench)
noinst_PROGRAMS = bench_step_vector bench_genomic_step_vector

bench_step_vector_SOURCES = bench_step_vector.cpp
bench_step_vector_CXXFLAGS = -Wall -std=c++17 -O2 \
                              -I$(top_srcdir)/src \
                              -I$(srcdir)/nanobench
bench_step_vector_LDADD = $(top_builddir)/src/libgcat.a

bench_genomic_step_vector_SOURCES = bench_genomic_step_vector.cpp
bench_genomic_step_vector_CXXFLAGS = -Wall -std=c++17 -O2 \
                                      -I$(top_srcdir)/src \
                                      -I$(srcdir)/nanobench
bench_genomic_step_vector_LDADD = $(top_builddir)/src/libgcat.a
