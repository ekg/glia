#include "cigar.h"
#include "utility.h"
#include "api/BamAlignment.h"

using namespace std;
using namespace BamTools;

struct AlignmentStats {
    int mismatches;
    int gaps;
    int gapslen;
    int softclips;
    int mismatch_qsum;
    int softclip_qsum;
AlignmentStats(void)
      : mismatches(0),
        gaps(0),
        gapslen(0),
        softclips(0),
        mismatch_qsum(0),
        softclip_qsum(0)
        { }
};

void countMismatchesAndGaps(
    BamAlignment& alignment,
    Cigar& cigar,
    string& referenceSequence,
    long int& referenceStart,
    AlignmentStats& stats,
    bool debug);

