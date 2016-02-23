 ## Process this file with automake to produce Makefile.in

bin_PROGRAMS = \
	bin/ace \
	bin/qls \
	bin/qgt

bin_ace_SOURCES = \
	src/mainACE.cpp \
	src/algorithm.cpp \
	src/algorithm.h \
	src/dataStructures.h \
	src/inverse.cpp \
	src/inverse.h \
	src/io.cpp \
	src/io.h \
	src/monteCarlo.cpp \
	src/monteCarlo.h \
	src/mtrnd.h \
	src/sparseInverse.cpp \
	src/sparseInverse.h \
	src/tools.cpp \
	src/tools.h

bin_qls_SOURCES = \
	src/mainQLS.cpp \
	src/io.cpp \
	src/io.h \
	src/monteCarlo.cpp \
	src/monteCarlo.h \
	src/mtrnd.h \
	src/qLearnSparse.cpp \
	src/qLearnSparse.h \
	src/tools.cpp \
	src/tools.h

bin_qgt_SOURCES = \
	src/mainQGT.cpp \
	src/io.cpp \
	src/io.h \
	src/monteCarlo.cpp \
	src/monteCarlo.h \
	src/mtrnd.h \
	src/qGenTest.cpp \
	src/qGenTest.h \
	src/tools.cpp \
	src/tools.h


GEN_TESTS = \
	tests/tools_innerProduct \
	tests/tools_L1 \
	tests/tools_L2 \
	tests/tools_LInfinity \
	tests/tools_insertInPlace \
	tests/tools_eraseInPlace \
	tests/inverse_getZ

TESTS = $(GEN_TESTS)
check_PROGRAMS = $(GEN_TESTS)
	
tests_tools_innerProduct_SOURCES = tests/tools_innerProduct.cpp src/tools.h src/tools.cpp
tests_tools_L1_SOURCES = tests/tools_L1.cpp src/tools.h src/tools.cpp
tests_tools_L2_SOURCES = tests/tools_L2.cpp src/tools.h src/tools.cpp
tests_tools_LInfinity_SOURCES = tests/tools_LInfinity.cpp src/tools.h src/tools.cpp
tests_tools_insertInPlace_SOURCES = tests/tools_insertInPlace.cpp src/tools.h src/tools.cpp
tests_tools_eraseInPlace_SOURCES = tests/tools_eraseInPlace.cpp src/tools.h src/tools.cpp
tests_inverse_getZ_SOURCES = tests/inverse_getZ.cpp src/tools.h src/tools.cpp src/inverse.h src/inverse.cpp src/sparseInverse.h src/sparseInverse.cpp

