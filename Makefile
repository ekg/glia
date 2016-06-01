HEADERS = dump.h \
	examples.h \
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
	utility.h \
	alignmentstats.h

SOURCES = dump.cpp \
	examples.cpp \
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
	utility.cpp \
	alignmentstats.cpp

OBJECTS= $(SOURCES:.cpp=.o)

BINS = glia

all: $(OBJECTS) $(BINS)

BAMTOOLS_ROOT=bamtools
BAMTOOLS_LIB_DIR=bamtools/lib

CXX = g++
CXXFLAGS = -O3 -D_FILE_OFFSET_BITS=64 -fpermissive #-pg
INCLUDES = -I$(BAMTOOLS_ROOT)/include -Ivcflib/src -Ivcflib/ -Igssw/src
LDFLAGS =
LIBS = -L./ -Lvcflib/ -Lvcflib/tabixpp/ -lbamtools -ltabix -lz -lm

FASTAHACK = fastahack/Fasta.o
VCFLIB = 	vcflib/src/Variant.o
VCFLIB_OBJS = vcflib/tabixpp/tabix.o \
			vcflib/tabixpp/bgzf.o \
			vcflib/smithwaterman/SmithWatermanGotoh.o \
			vcflib/smithwaterman/disorder.c \
			vcflib/smithwaterman/LeftAlign.o \
			vcflib/smithwaterman/Repeats.o \
			vcflib/smithwaterman/IndelAllele.o \

GSSW = gssw/obj/gssw.o


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

$(VCFLIB): libbamtools.a
	cd vcflib && $(MAKE)

$(GSSW): gssw/obj/gssw.o
	cd gssw && $(MAKE)

# glia build

%.o: %.cpp %.h $(VCFLIB)
	$(CXX) -c $(*F).cpp -o $@ $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(BINS): $(OBJECTS) $(SOURCES) $(HEADERS) $(FASTAHACK) jsoncpp.o $(VCFLIB) libbamtools.a $(GSSW)
	$(CXX) -o $@ $(INCLUDES) $(FASTAHACK) $(VCFLIB) $(VCFLIB_OBJS) $(OBJECTS) $(GSSW) $(LDFLAGS) $(CXXFLAGS) $(LIBS)

clean:
	rm -f $(BINS) $(OBJECTS)
	cd fastahack && $(MAKE) clean
	rm -f libbamtools.a
	cd bamtools/build && $(MAKE) clean
	cd gssw/src && $(MAKE) clean

clean-glia:
	rm -f $(BINS) $(OBJECTS)

.PHONY: clean all
