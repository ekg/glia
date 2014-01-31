#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "gliamodels.h"
#include "nodealign.h"
#include "traceback.h"
#include "gsw.h"
#include "examples.h"
#include "show.h"
#include "seqtools.h"
#include "parameters.h"
#include "construct.h"
#include "utility.h"
#include "alignmentstats.h"
#include "Variant.h"
#include "fastahack/Fasta.h"

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"

using namespace std;
//using namespace vcf;
using namespace BamTools;


void gswalign(vector<sn*>& nlist,
              string& read,
              string& qualities,
              Parameters& params,
              bt& backtrace,
              mbt& trace_report,
              int& score,
              string& strand) {

    sn* result_F;
    sn* result_R;
    int score_F=0;
    int score_R=0;
    mbt trace_report_F;
    mbt trace_report_R;
    bt backtrace_F;
    bt backtrace_R;

    if (params.debug) cerr << "aligning forward" << endl;

    result_F = gsw(read, nlist,
                   params.match, params.mism, params.gap);
    score_F = result_F->top_score.score;
    backtrace_F = master_backtrack(result_F, trace_report_F, read, qualities);
    vector<sn*>& nodes_F = trace_report_F.node_list;

    if (params.display_backtrace) {
        cout << "==== forward alignment ====" << endl;
        for (vector<sn*>::iterator n = nodes_F.begin();
             n != nodes_F.end(); ++n) {
            displayAlignment(*n, read);
        }
    } else if (params.display_all_nodes) {
        cout << "==== forward alignment ====" << endl;
        for (vector<sn*>::iterator n = nlist.begin();
             n != nlist.end(); ++n) {
            displayAlignment(*n, read);
        }
    }

    // check if the reverse complement provides a better alignment
    if (params.alignReverse) {
        if (params.debug) cerr << "aligning reverse" << endl;
        string readrc = reverseComplement(read);
        string qualitiesrc = qualities;
        reverse(qualitiesrc.begin(), qualitiesrc.end());
        result_R = gsw(readrc, nlist,
                       params.match, params.mism, params.gap);
        score_R = result_R->top_score.score;
        backtrace_R = master_backtrack(result_R, trace_report_R, readrc, qualitiesrc);
        vector<sn*>& nodes_R = trace_report_R.node_list;
        if (params.display_backtrace) {
            cout << "==== reverse alignment ====" << endl;
            for (vector<sn*>::iterator n = nodes_R.begin();
                 n != nodes_R.end(); ++n) {
                displayAlignment(*n, readrc);
            }
        } else if (params.display_all_nodes) {
            cout << "==== reverse alignment ====" << endl;
            for (vector<sn*>::iterator n = nlist.begin();
                 n != nlist.end(); ++n) {
                displayAlignment(*n, readrc);
            }
        }
    }

    //displayNode(result);
    //displayAlignment(result);
    //displayAlignment(nlist[0]);

    if (score_F > score_R) {
        backtrace = backtrace_F;
        trace_report = trace_report_F;
        score = score_F;
        strand = "+";
    } else {
        backtrace = backtrace_R;
        trace_report = trace_report_R;
        score = score_R;
        strand = "-";
    }

    //cout << "x: " << trace_report.x << " y: " << trace_report.y << endl;

}

// one-off
void construct_dag_and_align_single_sequence(Parameters& params) {

    if (params.debug) {
        cout << "read: " << params.read_input << endl;
        //cout << "fastq file:" << params.fastq_file << endl;
        cout << "fasta reference:" << params.fasta_reference << endl;
        cout << "vcf file " << params.vcf_file << endl;
        cout << "target " << params.target << endl;
        cout << endl;
    }

    // get sequence of target
    FastaReference reference;
    reference.open(params.fasta_reference);
    FastaRegion target(params.target);
    string targetSequence = reference.getTargetSubSequence(target);

    // get variants in target
    vector<vcf::Variant> variants;
    vcf::VariantCallFile vcffile;

    if (!params.vcf_file.empty()) {
        vcffile.open(params.vcf_file);
        vcf::Variant var(vcffile);
    
        vcffile.setRegion(params.target);
        while (vcffile.getNextVariant(var)) {
            if (var.position + var.ref.length() <= target.stopPos) {
                variants.push_back(var);
            }
        }
    }

    long offset = max(target.startPos, 1); // start is -1 when coordinates are not specified

    // Declare the target DAG to align against.
    vector<sn*> nlist;
    constructDAG(nlist, targetSequence, target.startSeq, variants, offset);

    if (params.display_dag) {
        cout << "DAG generated from input variants:" << endl;
        //displayDAG(nlist.back());
        for (vector<sn*>::iterator n = nlist.begin(); n != nlist.end(); ++n) {
            cout << *n << endl;
        }
        cout << endl;
    }


    // run the alignment

    string read = params.read_input;
    string qualities(read.size(), shortInt2QualityChar(30));
    bt backtrace;
    mbt trace_report;
    int score;
    string strand;
    gswalign(nlist,
             read,
             qualities,
             params,
             backtrace,
             trace_report,
             score,
             strand);

    cout << score << " " << strand << " "
         << (trace_report.node->position - 1) + trace_report.x << " "
         << trace_report.fcigar
         << " seq:" << trace_report.x << " read:" << trace_report.y
         << " " << trace_report.gcigar << " " << trace_report.fcigar << endl;

    if (params.display_alignment) {
        string refseq;
        for (vector<sn*>::iterator n = trace_report.node_list.begin();
             n != trace_report.node_list.end(); ++n) {
            refseq.append((*n)->sequence);
        }
        refseq = refseq.substr(trace_report.x, read.size());
        cout << refseq << endl;
        if (strand == "+") {
            cout << read << endl;
        } else {
            cout << reverseComplement(read) << endl;
        }
    }
}

bool shouldRealign(BamAlignment& alignment,
                   string& ref,
                   long int offset,
                   Parameters& params,
                   AlignmentStats& stats) {

    if (allN(alignment.QueryBases)) {
        cerr << "not realigning because query is all Ns! " << alignment.Name << endl;
        return false;
    }
    if (!alignment.IsMapped()) {
        if (params.debug) {
            cerr << "realigning because read " << alignment.Name << " is not mapped " << endl;
        }
        return true;
    }

    Cigar cigar(alignment.CigarData);
    countMismatchesAndGaps(alignment, cigar, ref, offset, stats);
    if (stats.mismatch_qsum > params.mismatch_qsum_threshold
        || stats.softclip_qsum > params.softclip_qsum_threshold) {
        if (params.debug) {
            cerr << "realigning because read " << alignment.Name << " meets mismatch (" << stats.mismatch_qsum << " vs. "
                 << params.mismatch_qsum_threshold << ") or softclip ("
                 << stats.softclip_qsum << " vs. " << params.softclip_qsum_threshold << ") thresholds" << endl;
        }
        return true;
    } else {
        return false;
    }
}

bool acceptRealignment(BamAlignment& alignment,
                       string& ref,
                       long int offset,
                       Parameters& params,
                       AlignmentStats& stats) {

    Cigar cigar(alignment.CigarData);
    countMismatchesAndGaps(alignment, cigar, ref, offset, stats);
    if (stats.mismatch_qsum > params.mismatch_qsum_max
        || stats.softclip_qsum > params.softclip_qsum_max
        || stats.gaps > params.gap_count_max) {
        return false;
    } else {
        return true;
    }
}

void realign_bam(Parameters& params) {

    FastaReference reference;
    reference.open(params.fasta_reference);

    bool suppress_output = false;

    int dag_window_size = params.dag_window_size;
    
    // open BAM file
    BamReader reader;
    if (!reader.Open("stdin")) {
        cerr << "could not open stdin for reading" << endl;
        exit(1);
    }

    BamWriter writer;
    if (!params.dry_run && !writer.Open("stdout", reader.GetHeaderText(), reader.GetReferenceData())) {
        cerr << "could not open stdout for writing" << endl;
        exit(1);
    }

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<RefData> referenceSequences = reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    vcf::VariantCallFile vcffile;
    if (!params.vcf_file.empty()) {
        vcffile.open(params.vcf_file);
    } else {
        cerr << "realignment requires VCF file" << endl;
        exit(1);
    }
    vcf::Variant var(vcffile);

    BamAlignment alignment;
    map<long int, vector<BamAlignment> > alignmentSortQueue;

    // get alignment
    // assemble DAG in region around alignment
    // loop for each alignment in BAM:
    //     update DAG when current alignment gets close to edge of assembled DAG
    //     attempt to realign if read has a certain number of mismatches + gaps or softclips, weighted by basequal
    //     if alignment to DAG has fewer mismatches and gaps than original alignment, use it
    //         flatten read into reference space (for now just output alleles from VCF un-spanned insertions)
    //     write read to queue for streaming re-sorting (some positional change will occur)

    long int dag_start_position = 1;
    string currentSeqname;
    string ref;
    vector<sn*> nlist; // the DAG container

    int total_reads = 0;
    int total_realigned = 0;
    int total_improved = 0;
    bool emptyDAG = false; // if the dag is constructed over empty sequence
                           // such as when realigning reads mapped to all-N sequence
    if (params.debug) {
        cerr << "about to start processing alignments" << endl;
    }

    while (reader.GetNextAlignment(alignment)) {
        if (params.debug) {
            cerr << "processing alignment " << alignment.Name << ":" << alignment.Position << endl;
        }

        ++total_reads;

        BamAlignment originalAlignment = alignment;
        long unsigned int initialAlignmentPosition = alignment.Position;
        //if (dag_start_position == 1) {
        //    dag_start_position = max(1, (int)initialAlignmentPosition - dag_window_size/2);
        //}
        string& seqname = referenceIDToName[alignment.RefID];

        // should we construct a new DAG?  do so when 3/4 of the way through the current one
        // center on current position + 1/2 dag window
        // TODO check this scheme using some scribbles on paper
        if (seqname != currentSeqname
            || (alignment.Position + (alignment.QueryBases.size()/2)
                > (3*dag_window_size/4) + dag_start_position)) {

            // reset current sequence name, if different
            currentSeqname = seqname;

            // recenter DAG
            if (!nlist.empty()) {
                dag_start_position = dag_start_position + dag_window_size/2;
                dag_start_position = max(dag_start_position,
                                         (long int) (alignment.GetEndPosition() - dag_window_size/2));
            } else {
                dag_start_position = alignment.Position - dag_window_size/2;
            }

            // TODO get sequence length and use to bound noted window size (edge case)
            ref = reference.getSubSequence(seqname,
                                           max((long int) 0, dag_start_position - 1),
                                           dag_window_size); // 0/1 conversion

            // get variants for new DAG
            vector<vcf::Variant> variants;
            if (vcffile.setRegion(currentSeqname,
                                   dag_start_position,
                                   dag_start_position + ref.size())) {
                while (vcffile.getNextVariant(var)) {
                    if (params.debug) cerr << "getting variant at " << var.sequenceName << ":" << var.position << endl;
                    if (var.position + var.ref.length() <= dag_start_position + ref.size()
                        && var.position >= dag_start_position) {
                        variants.push_back(var);
                    }
                }
            }

            // clear nlist
            for (vector<sn*>::iterator s = nlist.begin(); s != nlist.end(); ++s) {
                delete *s;
            }
            nlist.clear();


            if (params.debug) { cerr << "constructing DAG" << endl; }
            // and build the DAG
            constructDAG(nlist,
                         ref,
                         currentSeqname,
                         variants,
                         dag_start_position);

            if (params.display_dag || params.debug) {
                cerr << "DAG generated from input variants over "
                     << seqname << ":" << dag_start_position << "-" << dag_window_size
                     << endl;
                //displayDAG(nlist.back());
                for (vector<sn*>::iterator n = nlist.begin(); n != nlist.end(); ++n) {
                    cerr << *n << endl;
                }
                cerr << endl;
            }

            if (nlist.size() == 1 && allN(nlist.front()->sequence)) {
                emptyDAG = true;
            } else {
                emptyDAG = false;
            }

        }

        AlignmentStats stats_before;
        bool was_mapped = alignment.IsMapped();
        bool has_realigned = false;

        if (!emptyDAG && shouldRealign(alignment, ref, dag_start_position, params, stats_before)) {

            ++total_realigned;

            if (params.debug) {
                cerr << "realigning: " << alignment.Name << " " << alignment.QueryBases << endl
                     << "to variant graph over "
                     << seqname
                     << ":" << dag_start_position
                     << "-" << dag_start_position + dag_window_size << endl;
            }

            try {

                Cigar cigar;
                string read = alignment.QueryBases;
                string qualities = alignment.Qualities;
                bt backtrace;
                mbt trace_report;
                int score;
                string strand;
                gswalign(nlist,
                         read,
                         qualities,
                         params,
                         backtrace,
                         trace_report,
                         score,
                         strand);

                if (params.dry_run) {

                    //if (strand == "-" && !alignment.IsMapped()) {
                    //    read = reverseComplement(trace_report.read);
                    //}
                    cout << read << endl
                         << trace_report.read << endl;
                    cout << score << " " << strand << " "
                         << (trace_report.node->position - 1) + trace_report.x << " "
                         << trace_report.fcigar
                         << " seq:" << trace_report.x << " read:" << trace_report.y
                         << " " << trace_report.gcigar << " " << trace_report.fcigar << endl;

                } else {

                    /*
                    if (strand == "-") {
                        read = reverseComplement(trace_report.read);
                    }
                      cerr << read << endl
                      << trace_report.read << endl;
                      cerr << score << " " << strand
                      << " seq:" << trace_report.x << " read:" << trace_report.y
                      << " " << trace_report.gcigar << " " << trace_report.fcigar << endl;
                    */
                    // TODO the qualities are not on the right side of the read
                    if (strand == "-" && alignment.IsMapped()) {
                        // if we're realigning, this is always true unless we swapped strands
                        alignment.SetIsReverseStrand(true);
                        //alignment.QueryBases = reverseComplement(trace_report.read);
                        alignment.QueryBases = trace_report.read;
                        alignment.Qualities = trace_report.qualities;
                        //reverse(alignment.Qualities.begin(), alignment.Qualities.end()); // reverse qualities
                    } else {
                        // nothing to do for forward strand---
                        // BAM is already 5'-3' except for reverse strand flag
                        alignment.QueryBases = trace_report.read;
                        alignment.Qualities = trace_report.qualities;
                    }
                    alignment.Position = (trace_report.node->position - 1) + trace_report.x;
                    alignment.SetIsMapped(true);
                    if (!alignment.MapQuality) {
                        alignment.MapQuality = 20; // horrible hack...
                    }

                    // check if somehow we've ended up with an indel at the ends
                    // if so, grab the reference sequence right beyond it and add
                    // a single match to the cigar, allowing variant detection methods
                    // to run on the results without internal modification
                    Cigar& cigar = trace_report.fcigar;

                    int flankSize = params.flatten_flank;
                    if (cigar.front().isIndel() ||
                        (cigar.front().isSoftclip() && cigar.at(1).isIndel())) {
                        alignment.Position -= flankSize;
                        string refBase = reference.getSubSequence(seqname, alignment.Position, flankSize);
                        if (cigar.front().isSoftclip()) {
                            alignment.QueryBases.erase(alignment.QueryBases.begin(),
                                                       alignment.QueryBases.begin()+cigar.front().length);
                            alignment.Qualities.erase(alignment.Qualities.begin(),
                                                       alignment.Qualities.begin()+cigar.front().length);
                            cigar.erase(cigar.begin());
                        }
                        alignment.QueryBases.insert(0, refBase);
                        alignment.Qualities.insert(0, string(flankSize, shortInt2QualityChar(30)));
                        Cigar newCigar; newCigar.push_back(CigarElement(flankSize, 'M'));
                        newCigar.append(trace_report.fcigar);
                        trace_report.fcigar = newCigar;
                    }
                    if (cigar.back().isIndel() ||
                        (cigar.back().isSoftclip() && cigar.at(cigar.size()-2).isIndel())) {
                        string refBase = reference.getSubSequence(seqname,
                                                                  alignment.Position
                                                                  + trace_report.fcigar.refLen(),
                                                                  flankSize);
                        if (cigar.back().isSoftclip()) {
                            alignment.QueryBases.erase(alignment.QueryBases.end()-cigar.back().length,
                                                       alignment.QueryBases.end());
                            alignment.Qualities.erase(alignment.Qualities.end()-cigar.back().length,
                                                      alignment.Qualities.end());
                            cigar.pop_back();
                        }
                        Cigar newCigar; newCigar.push_back(CigarElement(flankSize, 'M'));
                        trace_report.fcigar.append(newCigar);
                        alignment.QueryBases.append(refBase);
                        alignment.Qualities.append(string(flankSize, shortInt2QualityChar(30)));
                    }

                    trace_report.fcigar.toCigarData(alignment.CigarData);

                    AlignmentStats stats_after;
                    countMismatchesAndGaps(alignment, trace_report.fcigar, ref, dag_start_position, stats_after);
                    if ((!was_mapped || stats_before.softclip_qsum >= stats_after.softclip_qsum)
                        && acceptRealignment(alignment, ref, dag_start_position, params, stats_after)) {

                        // keep the alignment
                        // TODO require threshold of softclips to keep alignment (or count of gaps, mismatches,...)
                        ++total_improved;
                        has_realigned = true;
                    } else {
                        // reset to old version of alignment
                        has_realigned = false;
                        alignment = originalAlignment;
                    }

                }

            } catch (...) {
                cerr << "exception when realigning " << alignment.Name
                     << " at position " << referenceIDToName[alignment.RefID]
                     << ":" << alignment.Position
                     << " " << alignment.QueryBases << endl;
                // reset to original alignment
                has_realigned = false;
                alignment = originalAlignment;
            }
        }

        if (!params.dry_run) {
            if (!params.only_realigned || has_realigned) {
                alignmentSortQueue[alignment.Position].push_back(alignment);
            }
            // ensure correct order if alignments move
            long int maxOutputPos = initialAlignmentPosition - dag_window_size;
            map<long int, vector<BamAlignment> >::iterator p = alignmentSortQueue.begin();
            for ( ; p != alignmentSortQueue.end(); ++p) {
                // except if we are running in unsorted mode, stop when we are at the window size
                if (!params.unsorted_output && p->first > maxOutputPos) {
                    break; // no more to do
                } else {
                    for (vector<BamAlignment>::iterator a = p->second.begin(); a != p->second.end(); ++a) {
                        writer.SaveAlignment(*a);
                    }
                }
            }
            if (p != alignmentSortQueue.begin()) {
                alignmentSortQueue.erase(alignmentSortQueue.begin(), p);
            }
        }
    } // end GetNextAlignment loop

    if (!params.dry_run) {
        map<long int, vector<BamAlignment> >::iterator p = alignmentSortQueue.begin();
        for ( ; p != alignmentSortQueue.end(); ++p) {
            for (vector<BamAlignment>::iterator a = p->second.begin(); a != p->second.end(); ++a)
                writer.SaveAlignment(*a);
        }
    }

    reader.Close();
    writer.Close();

    cerr << "total reads:\t" << total_reads << endl;
    cerr << "realigned:\t" << total_realigned << endl;
    cerr << "improved:\t" << total_improved << endl;

}

int main (int argc, char** argv) {
    
    Parameters params(argc, argv);

    if (!params.read_input.empty()) {
	// one-off read alignment
	// assemble local DAG and align read, report output
        construct_dag_and_align_single_sequence(params);
    } else if (params.realign_bam) {
        realign_bam(params);
    }

    return 0;

}
