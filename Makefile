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
	json/json.h \
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
	jsoncpp.cpp \
	split.cpp

OBJECTS= $(SOURCES:.cpp=.o)

BINS = clia

all: $(OBJECTS) $(BINS)

CXX = g++
CXXFLAGS = -O3 -D_FILE_OFFSET_BITS=64
INCLUDES = 
LDFLAGS =
LIBS =  -lm -L. -lbamtools -lz

BAMTOOLS_ROOT=bamtools
BAMTOOLS_LIB_DIR=bamtools/lib
FASTAHACK = fastahack/Fasta.o


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

split.o: split.h split.cpp
	$(CXX) $(CFLAGS) -c split.cpp

# clia build

$(OBJECTS): $(SOURCES) $(HEADERS)
	$(CXX) -c -o $@ $(*F).cpp $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(BINS): $(BIN_SOURCES) $(OBJECTS) libbamtools.a $(FASTAHACK)
	$(CXX) $(OBJECTS) -o $@ $(INCLUDES) $(FASTAHACK) $(LDFLAGS) $(CXXFLAGS)

clean:
	rm -f $(BINS) $(OBJECTS)
	cd fastahack && $(MAKE) clean
	cd bamtools/build && $(MAKE) clean

.PHONY: clean all
