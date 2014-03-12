#include "alignmentstats.h"


void countMismatchesAndGaps(
    BamAlignment& alignment,
    Cigar& cigar,
    string& referenceSequence,
    long int& referenceStart,
    AlignmentStats& stats,
    bool debug
    ) {

    int sp = alignment.Position - referenceStart;
    int rp = 0;

    int softclip_start = (cigar.begin()->type == 'S' ? cigar.begin()->length : 0);
    if (debug) {
        cerr << "pos  : " << alignment.Position << endl;
        cerr << "ref  : " << referenceSequence << endl;
        cerr << "ref  : " << string(softclip_start, ' ') << referenceSequence.substr(sp, cigar.refLen()) << endl;
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
    }

    if (alignment.Qualities.size() != alignment.QueryBases.size()) {
        cerr << "glia error: alignment has mismatch between length of aligned bases ("
             << alignment.QueryBases.size() << ") and length of qualities ("
             << alignment.Qualities.size() << ")" << endl;
        exit(1);
    }

    for (Cigar::const_iterator c = cigar.begin();
         c != cigar.end(); ++c) {
        int l = c->length;
        char t = c->type;
        //cerr << l << t << " " << sp << " " << rp << endl;
        if (t == 'M') { // match or mismatch
            for (int i = 0; i < l; ++i) {
                if (alignment.QueryBases.at(rp) != referenceSequence.at(sp)) {
                    ++stats.mismatches;
                    //cerr << "trying to access qualities at " << rp << " they are " << alignment.Qualities.size() << " long" << endl;
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
                //cerr << "trying to access qualities at " << rp << " they are " << alignment.Qualities.size() << " long" << endl;
                stats.softclip_qsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
                ++rp;
            }
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            sp += l;
        }
    }

}

