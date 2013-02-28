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
	json/json.h

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
	jsoncpp.cpp

OBJECTS= $(SOURCES:.cpp=.o)

BINS = clia

all: $(OBJECTS) $(BINS)

CXX = g++
CXXFLAGS = -O3 -D_FILE_OFFSET_BITS=64
INCLUDES = 
LDFLAGS =

#SSW = ssw.o ssw_cpp.o

#ssw.o: ssw.h
#ssw_cpp.o:ssw_cpp.h

profiling:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -g" all

gprof:
	$(MAKE) CXXFLAGS="$(CXXFLAGS) -pg" all

$(OBJECTS): $(SOURCES) $(HEADERS)
	$(CXX) -c -o $@ $(*F).cpp $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

$(BINS): $(BIN_SOURCES) $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ $(INCLUDES) $(LDFLAGS) $(CXXFLAGS)

clean:
	rm -f $(BINS) $(OBJECTS)

.PHONY: clean all
