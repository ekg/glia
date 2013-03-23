HEADERS = dump.h \
	examples.h \
	fastq.h \
	getSeeds.h \
	ghash.h \
	gliamodels.h \
	gsw.h \
	matrices.h \
	main.h \
	nodealign.h \
	show.h \
	traceback.h \
	parameters.h \
	seqtools.h \
	split.h \
	construct.h \
	cigar.h \
	utility.h

SOURCES = dump.cpp \
	examples.cpp \
	fastq.cpp \
	getSeeds.cpp \
	ghash.cpp \
	gliamodels.cpp \
	gsw.cpp \
	matrices.cpp \
	main.cpp \
	nodealign.cpp \
	show.cpp \
	traceback.cpp \
	parameters.cpp \
	seqtools.cpp \
	split.cpp \
	construct.cpp \
	cigar.cpp \
	utility.cpp

OBJECTS= $(SOURCES:.cpp=.o)

BINS = clia

all: $(OBJECTS) $(BINS)

BAMTOOLS_ROOT=bamtools
BAMTOOLS_LIB_DIR=bamtools/lib

CXX = g++
CXXFLAGS = -O3 -D_FILE_OFFSET_BITS=64
INCLUDES = -I$(BAMTOOLS_ROOT)/include
LDFLAGS =
LIBS = -lz -lm -L./ -Lvcflib/tabixpp/ -lbamtools -ltabix

FASTAHACK = fastahack/Fasta.o
VCFLIB = vcflib/tabixpp/tabix.o \
	vcflib/tabixpp/bgzf.o \
	vcflib/smithwaterman/SmithWatermanGotoh.o \
	vcflib/smithwaterman/disorder.c \
	vcflib/smithwaterman/LeftAlign.o \
	vcflib/smithwaterman/Repeats.o \
	vcflib/smithwaterman/IndelAllele.o \
	vcflib/Variant.o

#SSW = ssw.o ssw_cpp.o

#ssw.o: ssw.h
#ssw_cpp.o:ssw_cpp.h

# profiling

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

# libraries

# builds bamtools static lib, and copies into root
libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)
	cp bamtools/lib/libbamtools.a ./

$(FASTAHACK):
	cd fastahack && $(MAKE)

$(VCFLIB):
	cd vcflib && $(MAKE)

# clia build

%.o: %.cpp %.h
	$(CXX) -c $(*F).cpp -o $@ $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(BINS): $(BIN_SOURCES) $(OBJECTS) $(SOURCES) $(HEADERS) libbamtools.a $(FASTAHACK) jsoncpp.o $(VCFLIB)
	$(CXX) -o $@ $(INCLUDES) $(FASTAHACK) $(VCFLIB) $(OBJECTS) $(LDFLAGS) $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(BINS) $(OBJECTS)
	cd fastahack && $(MAKE) clean
	cd bamtools/build && $(MAKE) clean
	rm libbamtools.a

clean-clia:
	rm -f $(BINS) $(OBJECTS)

.PHONY: clean all
