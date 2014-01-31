#include "alignmentstats.h"


void countMismatchesAndGaps(
    BamAlignment& alignment,
    Cigar& cigar,
    string& referenceSequence,
    long int& referenceStart,
    AlignmentStats& stats
    ) {

    int sp = alignment.Position - referenceStart + 1;
    int rp = 0;

    /*
    cerr << "ref  : " << referenceSequence << endl;
    cerr << "read : " << alignment.QueryBases << endl;
    cerr << "quals: " << alignment.Qualities << endl;
    cerr << "cigar: ";
    for (Cigar::const_iterator c = cigar.begin();
         c != cigar.end(); ++c) {
        int l = c->length;
        char t = c->type;
        cerr << l << t;
    }
    cerr << endl;
    cerr << "starts at " << sp << " of cached ref (" << referenceSequence.size() << "bp long)" << endl;
    */

    for (Cigar::const_iterator c = cigar.begin();
         c != cigar.end(); ++c) {
        int l = c->length;
        char t = c->type;
        //cerr << l << t << " " << sp << " " << rp << endl;
        if (t == 'M') { // match or mismatch
            for (int i = 0; i < l; ++i) {
                if (sp >= referenceSequence.length() ||
                      alignment.QueryBases.at(rp) != referenceSequence.at(sp)) {
                    ++stats.mismatches;
                    stats.mismatch_qsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
                }
                ++sp;
                ++rp;
            }
        } else if (t == 'D') { // deletion
            stats.gaps++;
            stats.gapslen += l;
            sp += l;  // update reference sequence position
        } else if (t == 'I') { // insertion
            stats.gaps++;
            stats.gapslen += l;
            rp += l;  // update read position
        } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
            stats.softclips += l;
            for (int i = 0; i < l; ++i) {
                stats.softclip_qsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
                ++rp;
            }
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            sp += l;
        }
    }

}

