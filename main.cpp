#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <google/sparse_hash_map>
//#include <json/json.h>
//#include <tr1/functional>

#include "gliamodels.h"
#include "nodealign.h"
#include "traceback.h"
#include "gsw.h"
#include "examples.h"
#include "fastq.h"
#include "show.h"
#include "ghash.h"
#include "getSeeds.h"
#include "jsreader.h"
#include "seqtools.h"
#include "parameters.h"
#include "construct.h"
#include "vcflib/Variant.h"
#include "fastahack/Fasta.h"

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"

using namespace std;
using namespace vcf;
using namespace BamTools;


short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

char shortInt2QualityChar(short i) {
    return static_cast<char>(i + 33);
}

// TODO: Move to GHASH   <-done?
int hashfasta(string fasta_file_name, int hashsize, vector<fasta_entry> &ref_genome) {	

	load_fasta_file(fasta_file_name, ref_genome);
	
	vector<fasta_entry>::iterator t;
       
	for (t = ref_genome.begin(); t != ref_genome.end(); ++t) {
		cout << t->name << endl;
		//cout << t->sequence << endl;
		hashcontig(*t, t->kmer_hash, hashsize);
		cout << t->name << " ...hash complete" << endl;
	
		sortContigHash(t->kmer_hash);
		cout << t->name << " ...sort complete" << endl;
				
	}
	
	return 0;
}

void countMismatchesAndGaps(
    BamAlignment& alignment,
    //vector<CigarOp>& cigarData,
    Cigar& cigar,
    string referenceSequence,
    int& mismatches,
    int& gaps,
    int& gapslen,
    int& softclips,
    int& mismatchQsum,
    int& softclipQsum
    ) {

    int sp = 0;
    int rp = 0;
    for (Cigar::const_iterator c = cigar.begin();
         c != cigar.end(); ++c) {
        int l = c->length;
        char t = c->type;
        if (t == 'M') { // match or mismatch
            for (int i = 0; i < l; ++i) {
                if (alignment.QueryBases.at(rp) != referenceSequence.at(sp)) {
                    ++mismatches;
                    mismatchQsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
                }
                ++sp;
                ++rp;
            }
        } else if (t == 'D') { // deletion
            ++gaps;
            gapslen += l;
            sp += l;  // update reference sequence position
        } else if (t == 'I') { // insertion
            ++gaps;
            gapslen += l;
            rp += l;  // update read position
        } else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
            softclips += l;
            for (int i = 0; i < l; ++i) {
                softclipQsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
                ++rp;
            }
        } else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
        } else if (t == 'N') { // skipped region in the reference not present in read, aka splice
            sp += l;
        }
    }

}

void gswalign(vector<sn*>& nlist,
              string& read,
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

    result_F = gsw(read, nlist,
                   params.match, params.mism, params.gap);
    score_F = result_F->top_score.score;
    backtrace_F = master_backtrack(result_F, trace_report_F);
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
        string readrc = reverseComplement(read);
        result_R = gsw(readrc, nlist,
                       params.match, params.mism, params.gap);
        score_R = result_R->top_score.score;
        backtrace_R = master_backtrack(result_R, trace_report_R);
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
    vector<Variant> variants;
    VariantCallFile vcffile;

    if (!params.vcf_file.empty()) {
        vcffile.open(params.vcf_file);
        Variant var(vcffile);
    
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
    bt backtrace;
    mbt trace_report;
    int score;
    string strand;
    gswalign(nlist,
             read,
             params,
             backtrace,
             trace_report,
             score,
             strand);

    cout << score << " " << strand
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

bool shouldRealign(BamAlignment& alignment, string& ref, long int offset, Parameters& params) {
    return true;
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

    suppress_output = true;
    BamWriter writer;
    if (!suppress_output && !writer.Open("stdout", reader.GetHeaderText(), reader.GetReferenceData())) {
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

    VariantCallFile vcffile;
    if (!params.vcf_file.empty()) {
        vcffile.open(params.vcf_file);
    } else {
        cerr << "realignment requires VCF file" << endl;
        exit(1);
    }
    Variant var(vcffile);

    BamAlignment alignment;
    int flanking_window = dag_window_size/2; // streaming sort half the dag window size
    map<long unsigned int, vector<BamAlignment> > alignmentSortQueue;

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

    while (reader.GetNextAlignment(alignment)) {

        long unsigned int initialAlignmentPosition = alignment.Position;
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
            }

            // TODO get sequence length and use to bound noted window size (edge case)
            ref = reference.getSubSequence(seqname, dag_start_position - 1, dag_window_size); // 0/1 conversion

            // get variants for new DAG
            vector<Variant> variants;
            if (!vcffile.setRegion(currentSeqname,
                                   dag_start_position,
                                   dag_start_position + ref.size())) {
                cerr << "could not set region on VCF file to " << currentSeqname << ":"
                     << dag_start_position << "-" << dag_start_position + ref.size()
                     << endl;
                exit(1);
            } else {
                while (vcffile.getNextVariant(var)) {
                    if (var.position + var.ref.length() <= dag_start_position + ref.size()) {
                        variants.push_back(var);
                    }
                }
            }

            // clear nlist
            for (vector<sn*>::iterator s = nlist.begin(); s != nlist.end(); ++s) {
                delete *s;
            }
            nlist.clear();

            // and build the DAG
            constructDAG(nlist,
                         ref,
                         currentSeqname,
                         variants,
                         dag_start_position);

            if (params.display_dag) {
                cout << "DAG generated from input variants over "
                     << seqname << ":" << dag_start_position << "-" << dag_window_size
                     << endl;
                //displayDAG(nlist.back());
                for (vector<sn*>::iterator n = nlist.begin(); n != nlist.end(); ++n) {
                    cout << *n << endl;
                }
                cout << endl;
            }

        }

        if (!alignment.IsMapped() || shouldRealign(alignment, ref, dag_start_position, params)) {

            try {

                Cigar cigar;
                string read = alignment.QueryBases;
                bt backtrace;
                mbt trace_report;
                int score;
                string strand;
                gswalign(nlist,
                         read,
                         params,
                         backtrace,
                         trace_report,
                         score,
                         strand);

                if (params.dry_run) {
                    cout << read << endl;
                    cout << score << " " << strand
                         << " seq:" << trace_report.x << " read:" << trace_report.y
                         << " " << trace_report.gcigar << " " << trace_report.fcigar << endl;
                }

            } catch (...) {
                cerr << "exception when realigning " << alignment.Name
                     << " at position " << referenceIDToName[alignment.RefID]
                     << ":" << alignment.Position
                     << " " << alignment.QueryBases << endl;
            }
        }

        if (!suppress_output) {
            alignmentSortQueue[alignment.Position].push_back(alignment);
            // ensure correct order if alignments move
            if (initialAlignmentPosition > (unsigned int) flanking_window) {
                long unsigned int maxOutputPos = initialAlignmentPosition - flanking_window;
                map<long unsigned int, vector<BamAlignment> >::iterator p = alignmentSortQueue.begin();
                for ( ; p != alignmentSortQueue.end(); ++p) {
                    if (p->first > maxOutputPos) {
                        break; // no more to do
                    } else {
                        //for (vector<BamAlignment>::iterator a = p->second.begin(); a != p->second.end(); ++a)
                            //writer.SaveAlignment(*a);
                    }
                }
                if (p != alignmentSortQueue.begin())
                    alignmentSortQueue.erase(alignmentSortQueue.begin(), p);
            }
        }

    }

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
