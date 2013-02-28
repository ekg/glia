HEADERS = dump.h \
	examples.h \
	fastq.h \
	getSeeds.h \
	ghash.h \
	gliamodels.h \
	gsw.h \
	matrices.h \
	nodealign.h \
	show.h \
	traceback.h \
	parameters.h \
	seqtools.h \
	split.h

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
	split.cpp

OBJECTS= $(SOURCES:.cpp=.o)

BINS = clia

all: $(OBJECTS) $(BINS)

CXX = g++
CXXFLAGS = -O3 -D_FILE_OFFSET_BITS=64
INCLUDES = 
LDFLAGS =
LIBS = -lz -lm -L./ -L../vcflib/tabixpp/ -L$(BAMTOOLS_ROOT)/lib -ltabix

BAMTOOLS_ROOT=bamtools
BAMTOOLS_LIB_DIR=bamtools/lib
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

#$(OBJECTS): $(SOURCES) $(HEADERS)

%.o: %.cpp %.h
	$(CXX) -c -o $@ $(*F).cpp $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(BINS): $(BIN_SOURCES) $(OBJECTS) $(SOURCES) $(HEADERS) libbamtools.a $(FASTAHACK) jsoncpp.o $(VCFLIB)
	$(CXX) $(OBJECTS) $(VCFLIB) -o $@ $(INCLUDES) $(FASTAHACK) $(LDFLAGS) $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(BINS) $(OBJECTS)
	cd fastahack && $(MAKE) clean
	cd bamtools/build && $(MAKE) clean
	rm libbamtools.a

.PHONY: clean all
