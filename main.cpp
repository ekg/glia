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

#include "gssw.h"

using namespace std;
//using namespace vcf;
using namespace BamTools;


gssw_graph_mapping*
gswalign(gssw_graph* graph,
         //vector<Cigar>& cigars, // by node id
         //vector<long int>& refpositions,  // by node id
         ReferenceMappings& ref_map,
         string& read,
         string& qualities,
         Parameters& params,
         long int& position,
         int& score,
         Cigar& flat_cigar,
         string& strand,
         int8_t* nt_table,
         int8_t* score_matrix) {

    sn* result_F;
    sn* result_R;
    int score_F=0;
    int score_R=0;
    mbt trace_report_F;
    mbt trace_report_R;
    bt backtrace_F;
    bt backtrace_R;
    gssw_graph_mapping* gmf = NULL;
    gssw_graph_mapping* gmr = NULL;

    if (params.debug) cerr << "aligning forward" << endl;

    const char* cread = read.c_str();
    gssw_graph_fill(graph, cread, nt_table, score_matrix, params.gap_open, params.gap_extend, 15, 2);
    gmf = gssw_graph_trace_back (graph,
                                 cread,
                                 read.length(),
                                 params.match,
                                 params.mism,
                                 params.gap_open,
                                 params.gap_extend);

    if (params.display_backtrace) {
        cerr << "==== forward alignment ====" << endl;
        cerr << graph_mapping_to_string(gmf) << endl;
    }
    if (params.display_all_nodes) {
        cerr << "==== forward alignment ====" << endl;
        gssw_graph_print_score_matrices(graph, read.c_str(), read.length());
    }

    score_F = gmf->score;

    // check if the reverse complement provides a better alignment
    if (params.alignReverse) {
        if (params.debug) cerr << "aligning reverse" << endl;
        string readrc = reverseComplement(read);
        string qualitiesrc = qualities;
        reverse(qualitiesrc.begin(), qualitiesrc.end());
        const char* creadrc = readrc.c_str();
        gssw_graph_fill(graph, creadrc, nt_table, score_matrix, params.gap_open, params.gap_extend, 15, 2);
        gmr = gssw_graph_trace_back (graph,
                                     creadrc,
                                     readrc.length(),
                                     params.match,
                                     params.mism,
                                     params.gap_open,
                                     params.gap_extend);

        if (params.display_backtrace) {
            cerr << "==== reverse alignment ====" << endl;
            cerr << graph_mapping_to_string(gmr) << endl;
            //gssw_print_graph_mapping(gmr);
        }
        if (params.display_all_nodes) {
            cerr << "==== reverse alignment ====" << endl;
            gssw_graph_print_score_matrices(graph, read.c_str(), read.length());
            //gssw_print_graph_mapping(gmr);
        }
        score_R = gmr->score;

    }

    gssw_graph_mapping* gm = NULL;
    if (score_F > score_R) {
        score = score_F;
        gm = gmf;
        strand = "+";
        if (params.alignReverse) {
            gssw_graph_mapping_destroy(gmr);
        }
    } else {
        score = score_R;
        gm = gmr;
        strand = "-";
        gssw_graph_mapping_destroy(gmf);
        // reverse complement things
        read = reverseComplement(read);
        reverse(qualities.begin(), qualities.end());
    }

    // print graph mapping
    if (params.debug) {
        cerr << graph_mapping_to_string(gm) << endl;
    }

    // determine the reference-relative position
    ReferenceMapping& rm = ref_map.get_node(gm->cigar.elements[0].node);
    if (rm.is_ref()) {
        position = gm->position + rm.ref_position;
    } else { // flatten to previous position
        position = rm.ref_position;
    }

    int read_pos = 0;
    // flatten the winning alignment
    for (int i = 0; i < gm->cigar.length; ++i) {
        gssw_node* n = gm->cigar.elements[i].node;
        gssw_cigar* c = gm->cigar.elements[i].cigar;
        Cigar graph_relative_cigar = Cigar(c);

        //Cigar& ref_relative_cigar = cigars.at(n->id);

        ReferenceMapping& rm = ref_map.get_node(gm->cigar.elements[i].node);
        Cigar& ref_relative_cigar = rm.cigar;

        // if we are not in the reference coordinate space, flatten
        // to enable downstream detection algorithms to work on the BAM stream
        if (ref_relative_cigar.refLen() != ref_relative_cigar.readLen()) {

            flat_cigar.append(ref_relative_cigar);
            // first check that we don't simply match the reference
            // although map to a node in the graph which is not in the reference path

            // flatten things back into the reference space
            //if (params.flatten) {
            if (true) {
                string s = string(n->seq);
                //cerr << read.substr(read_pos, graph_relative_cigar.readLen()) << endl << s << endl;
                read.replace(read_pos, graph_relative_cigar.readLen(), s);
                short average_qual = (short) averageQuality(qualities.substr(read_pos, graph_relative_cigar.readLen()));
                qualities.replace(read_pos, graph_relative_cigar.readLen(),
                                  string(s.size(), shortInt2QualityChar(average_qual)));
                //cerr << read << endl << qualities << endl;
            }
        } else {
            flat_cigar.append(graph_relative_cigar);
        }
        read_pos += graph_relative_cigar.readLen();

        // do the edges
        if (i < gm->cigar.length - 1) {
            // check for edge mapping, e.g. deletion
            //cerr << "checking for next edge mapping (e.g. deletion)" << endl;
            gssw_node* next = gm->cigar.elements[i+1].node;
            //cerr << "getting edge " << n << " to " << next << endl;
            ReferenceMapping& rm = ref_map.get_edge(n, next);
            Cigar& edge_ref_relative_cigar = rm.cigar;
            //cerr << "new cigar from " << n->id << " to " << next->id  <<" is " << edge_ref_relative_cigar << endl;
            if (edge_ref_relative_cigar.refLen() != edge_ref_relative_cigar.readLen()) {
                // NB these should only be deletions, as edges won't represent inserted sequences
                flat_cigar.append(edge_ref_relative_cigar);
                //cerr << "flattening! " << graph_relative_cigar.readLen()  << endl;
            }
        }

        /*
        cerr << "building cigar: " << read_pos << endl
             << flat_cigar << endl
             << read << endl;
        */
    }
    //cerr << flat_cigar << endl;

/*
        if (node->isref) { // if we're in the reference coordinate space
            backtrace.cigar = cigar;
        } else {

            backtrace.x = 0; // matches to start of variant now
            trace_report.read.replace(y, cigar.readLen(), node->sequence);
            backtrace.cigar = node->cigar;
            reverse(backtrace.cigar.begin(), backtrace.cigar.end()); // re-reverse
            trace_report.qualities.replace(y, cigar.readLen(),
                                           string(node->sequence.size(), shortInt2QualityChar(30)));
            //trace_report.read.replace(y, backtrace.cigar.readLen(), node->sequence.substr(0, x));
            // and 1bp of reference before the implied divergent sequence
        }
*/

    //cout << "x: " << trace_report.x << " y: " << trace_report.y << endl;
    return gm;

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
    //vector<Cigar> cigars;
    //vector<long int> refpositions;
    ReferenceMappings ref_map;
    gssw_graph* graph = gssw_graph_create(0);
    int8_t* nt_table = gssw_create_nt_table();
	int8_t* mat = gssw_create_score_matrix(params.match, params.mism);
    constructDAGProgressive(graph,
                            ref_map,
                            targetSequence,
                            target.startSeq,
                            variants,
                            offset,
                            nt_table,
                            mat,
                            params.flat_input_vcf);

    if (params.display_dag) {
        cout << "DAG generated from input variants:" << endl;
    }


    // run the alignment

    string read = params.read_input;
    string qualities(read.size(), shortInt2QualityChar(30));
    int score;
    long int position;
    string strand;
    Cigar flat_cigar;
    gssw_graph_mapping* gm = gswalign(graph,
                                      ref_map,
                                      read,
                                      qualities,
                                      params,
                                      position,
                                      score,
                                      flat_cigar,
                                      strand,
                                      nt_table,
                                      mat);
    cerr << graph_mapping_to_string(gm) << endl;
    gssw_graph_mapping_destroy(gm);

    /*
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
    */
}

bool shouldRealign(BamAlignment& alignment,
                   string& ref,
                   long int offset,
                   Parameters& params,
                   AlignmentStats& stats) {

    if (allN(alignment.QueryBases)) {
        if (params.debug) {
            cerr << "not realigning because query is all Ns! " << alignment.Name << endl;
        }
        return false;
    }
    if (!alignment.IsMapped()) {
        if (params.debug) {
            cerr << "realigning because read " << alignment.Name << " is not mapped " << endl;
        }
        return true;
    }
    
    if (alignment.CigarData.empty()) {
        cerr << "realigning because alignment " << alignment.Name << " @ " << alignment.Position
             << " has empty (or corrupted?) CIGAR" << endl;
        return true;
    }

    Cigar cigar(alignment.CigarData);
    countMismatchesAndGaps(alignment, cigar, ref, offset, stats, params.debug);

    if (stats.mismatch_qsum >= params.mismatch_qsum_threshold
        || stats.softclip_qsum >= params.softclip_qsum_threshold
        || stats.gaps >= params.gap_count_threshold
        || stats.gapslen >= params.gap_length_threshold) {
        if (params.debug) {
            cerr << "realigning because read " << alignment.Name
                 << " meets mismatch (q" << stats.mismatch_qsum << " in " << stats.mismatches << ")" //<< " vs. " << params.mismatch_qsum_threshold << "),"
                 << " softclip (q" << stats.softclip_qsum << " in " << stats.softclips << ")" //<< " vs. " << params.softclip_qsum_threshold << "),"
                 << " gap count (" << stats.gaps << ")" //" vs. " << params.gap_count_threshold << "),"
                 << " or gap length (" << stats.gapslen << ")" //<< " vs. " << params.gap_length_threshold << ") "
                 << " thresholds" << endl;
            cerr << cigar << endl;
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

    //Cigar cigar(alignment.CigarData);
    //countMismatchesAndGaps(alignment, cigar, ref, offset, stats);
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
        if (!vcffile.open(params.vcf_file)) {
            cerr << "could not open VCF file " << params.vcf_file << endl;
            exit(1);
        }
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

    long int dag_start_position = 0;
    string currentSeqname;
    string ref;
    //vector<Cigar> cigars; // contains the Cigar strings of nodes in the graph
    //vector<long int> refpositions; // contains the reference start coords of nodes in the graph
    ReferenceMappings ref_map;
    gssw_graph* graph = gssw_graph_create(0);
    int8_t* nt_table = gssw_create_nt_table();
    int8_t* mat = gssw_create_score_matrix(params.match, params.mism);

    int total_reads = 0;
    int total_realigned = 0;
    int total_improved = 0;
    bool emptyDAG = false; // if the dag is constructed over empty sequence
                           // such as when realigning reads mapped to all-N sequence
    if (params.debug) {
        cerr << "about to start processing alignments" << endl;
    }

    while (reader.GetNextAlignment(alignment)) {

        string& seqname = referenceIDToName[alignment.RefID];

        if (params.debug) {
            cerr << "--------------------------------------------" << endl
                 << "processing alignment " << alignment.Name << " at "
                 << seqname << ":" << alignment.Position << endl;
        }

        /*
        if (!alignment.IsMapped() && graph->size == 0) {
            if (params.debug) {
                cerr << "unable to build DAG using unmapped read "
                     << alignment.Name << " @ "
                     << seqname << ":" << alignment.Position
                     << " no previous mapped read found and DAG currently empty" << endl;
            }
            alignmentSortQueue[dag_start_position+dag_window_size].push_back(alignment);
            continue;
        }
        */

        ++total_reads;

        BamAlignment originalAlignment = alignment;
        long unsigned int initialAlignmentPosition = alignment.Position;
        //if (dag_start_position == 1) {
        //    dag_start_position = max(1, (int)initialAlignmentPosition - dag_window_size/2);
        //}

        // should we construct a new DAG?  do so when 3/4 of the way through the current one
        // center on current position + 1/2 dag window
        // TODO check this scheme using some scribbles on paper
        // alignment.IsMapped()
        if ((seqname != currentSeqname
             || ((alignment.Position + (alignment.QueryBases.size()/2)
                  > (3*dag_window_size/4) + dag_start_position)))
            && alignment.Position < reference.sequenceLength(seqname)) {

            if (seqname != currentSeqname) {
                if (params.debug) {
                    cerr << "switched ref seqs" << endl;
                }
                dag_start_position = max((long int) 0,
                                         (long int) (alignment.GetEndPosition() - dag_window_size/2));
            // recenter DAG
            } else if (!ref_map.empty()) {
                dag_start_position = dag_start_position + dag_window_size/2;
                dag_start_position = max(dag_start_position,
                                         (long int) (alignment.GetEndPosition() - dag_window_size/2));
            } else {
                dag_start_position = alignment.Position - dag_window_size/2;
            }
            dag_start_position = max((long int)0, dag_start_position);

            // TODO get sequence length and use to bound noted window size (edge case)
            //cerr << "getting ref " << seqname << " " << max((long int) 0, dag_start_position) << " " << dag_window_size << endl;

            // get variants for new DAG
            vector<vcf::Variant> variants;
            if (!vcffile.setRegion(seqname,
                                   dag_start_position + 1,
                                   dag_start_position + dag_window_size)) {
                // this is not necessarily an error; there should be a better way to check for VCF file validity
                /*
                cerr << "could not set region on VCF file to " << currentSeqname << ":"
                     << dag_start_position << "-" << dag_start_position + ref.size()
                     << endl;
                */
                //exit(1);
            } else {

                // check first variant
                if (vcffile.getNextVariant(var)) {
                    while (var.position <= dag_start_position + 1) {
                        //cerr << "var position == dag_start_position " << endl;
                        dag_start_position -= 1;
                        vcffile.setRegion(seqname,
                                          dag_start_position + 1,
                                          dag_start_position + dag_window_size);
                        if (!vcffile.getNextVariant(var)) { break; }
                    }
                }

                vcffile.setRegion(seqname,
                                  dag_start_position + 1,
                                  dag_start_position + dag_window_size);

                while (vcffile.getNextVariant(var)) {
                    if (params.debug) cerr << "getting variant at " << var.sequenceName << ":" << var.position << endl;
                    //cerr << var.position << " + " << var.ref.length() << " <= " << dag_start_position << " + " << dag_window_size << endl;
                    //cerr << var.position << " >= " << dag_start_position << endl;
                    if (var.position + var.ref.length() <= dag_start_position + dag_window_size
                        && var.position >= dag_start_position) {
                        variants.push_back(var);
                    }
                }

            }

            //cerr << "dag_start_position " << dag_start_position << endl;
            ref = reference.getSubSequence(seqname,
                                           max((long int) 0, dag_start_position),
                                           dag_window_size); // 0/1 conversion

            // clear graph and metadata
            ref_map.clear();
            //cigars.clear();
            //refpositions.clear();
            gssw_graph_destroy(graph);

            if (params.debug) { cerr << "constructing DAG" << endl; }
            // and build the DAG
            graph = gssw_graph_create(0);
            constructDAGProgressive(graph,
                                    ref_map,
                                    ref,
                                    seqname,
                                    variants,
                                    dag_start_position,
                                    nt_table,
                                    mat,
                                    params.flat_input_vcf);

            if (params.debug) {
                cerr << "graph has " << graph->size << " nodes" << endl;
                cerr << "DAG generated from input variants over "
                     << seqname << ":" << dag_start_position << "-" << dag_start_position + dag_window_size
                     << endl;
            }
            if (params.display_dag) {
                gssw_graph_print(graph);
                /*
                for (Backbone::iterator b = backbone.begin(); b != backbone.end(); ++b) {
                    cout << b->first << " "
                         << b->first->id << " "
                         << b->second.ref_position << " "
                         << b->second.cigar << endl
                         << b->first->seq << endl;
                }
                */
            }

            if (graph->size == 1 && allN(ref) || graph->size == 0) {
                if (params.debug) {
                    cerr << "DAG is empty (1 node, all N).  Alignment is irrelevant." << endl;
                }
                emptyDAG = true;
            } else {
                emptyDAG = false;
            }

        }

        AlignmentStats stats_before;
        bool was_mapped = alignment.IsMapped();
        bool has_realigned = false;
        if (was_mapped) {
            if (dag_start_position + dag_window_size < alignment.GetEndPosition()) {
                ref = reference.getSubSequence(seqname,
                                               max((long int) 0, dag_start_position),
                                               alignment.GetEndPosition() - dag_start_position); // 0/1 conversion
            }
        }

        if (params.debug) {
            if (emptyDAG) {
                cerr << "cannot realign against empty (all-N single node) graph" << endl;
            }
        }

        if (!emptyDAG && shouldRealign(alignment, ref, dag_start_position, params, stats_before)) {

            ++total_realigned;

            if (params.debug) {
                cerr << "realigning: " << alignment.Name
                     << " " << alignment.QueryBases << endl
                     << " aligned @ " << alignment.Position
                     << " to variant graph over "
                     << seqname
                     << ":" << dag_start_position
                     << "-" << dag_start_position + dag_window_size << endl;
            }

            try {

                Cigar flat_cigar;
                string read = alignment.QueryBases;
                string qualities = alignment.Qualities;
                int score;
                long int position;
                string strand;
                gssw_graph_mapping* gm =
                    gswalign(graph,
                             ref_map,
                             read,
                             qualities,
                             params,
                             position,
                             score,
                             flat_cigar,
                             strand,
                             nt_table,
                             mat);
                //
                gssw_graph_mapping_destroy(gm);

                if (params.dry_run) {

                    if (strand == "-" && !alignment.IsMapped()) {
                        read = reverseComplement(read);
                    }
                    cout << read << endl;
                    cout << graph_mapping_to_string(gm) << endl;
                    cout << score << " " << strand << " "
                         << position << " "
                         << flat_cigar << endl;

                } else {

                    /*
                    if (strand == "-") {
                        read = reverseComplement(trace_report.read);
                    }
                   */
 
                    // TODO the qualities are not on the right side of the read
                    if (strand == "-" && alignment.IsMapped()) {
                        // if we're realigning, this is always true unless we swapped strands
                        alignment.SetIsReverseStrand(true);
                        //reverse(alignment.Qualities.begin(), alignment.Qualities.end()); // reverse qualities
                    }
                    //alignment.QueryBases = reverseComplement(trace_report.read);
                    alignment.QueryBases = read;
                    alignment.Qualities = qualities;

                    alignment.Position = position;// + 1;// + 1;//(trace_report.node->position - 1) + trace_report.x;
                    alignment.SetIsMapped(true);
                    if (!alignment.MapQuality) {
                        alignment.MapQuality = 20; // horrible hack...  at least approximate with alignment mismatches against graph
                    }

                    // check if somehow we've ended up with an indel at the ends
                    // if so, grab the reference sequence right beyond it and add
                    // a single match to the cigar, allowing variant detection methods
                    // to run on the results without internal modification
                    Cigar& cigar = flat_cigar;
                    //cerr << flat_cigar << " " << flat_cigar.readLen() << " " << flat_cigar.refLen() << endl;
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
                        newCigar.append(flat_cigar);
                        flat_cigar = newCigar;
                    }
                    if (cigar.back().isIndel() ||
                        (cigar.back().isSoftclip() && cigar.at(cigar.size()-2).isIndel())) {
                        string refBase = reference.getSubSequence(seqname,
                                                                  alignment.Position
                                                                  + flat_cigar.refLen(),
                                                                  flankSize);
                        if (cigar.back().isSoftclip()) {
                            alignment.QueryBases.erase(alignment.QueryBases.end()-cigar.back().length,
                                                       alignment.QueryBases.end());
                            alignment.Qualities.erase(alignment.Qualities.end()-cigar.back().length,
                                                      alignment.Qualities.end());
                            cigar.pop_back();
                        }
                        Cigar newCigar; newCigar.push_back(CigarElement(flankSize, 'M'));
                        flat_cigar.append(newCigar);
                        //flat_cigar.append(newCigar);
                        alignment.QueryBases.append(refBase);
                        alignment.Qualities.append(string(flankSize, shortInt2QualityChar(30)));
                    }

                    flat_cigar.toCigarData(alignment.CigarData);
                    //cerr << flat_cigar << " " << flat_cigar.readLen() << " " << flat_cigar.refLen() << endl;

                    if (dag_start_position + dag_window_size < alignment.GetEndPosition()) {
                        ref = reference.getSubSequence(seqname,
                                                       max((long int) 0, dag_start_position),
                                                       alignment.GetEndPosition() - dag_start_position); // 0/1 conversion
                    }

                    AlignmentStats stats_after;
                    countMismatchesAndGaps(alignment, flat_cigar, ref, dag_start_position, stats_after, params.debug);
                    /*
                    if ((!was_mapped || (stats_before.softclip_qsum >= stats_after.softclip_qsum
                                         && stats_before.mismatch_qsum >= stats_after.mismatch_qsum))
                         && acceptRealignment(alignment, ref, dag_start_position, params, stats_after)) {
                    */
                    /*
                    if ((!was_mapped || (stats_before.softclip_qsum + stats_before.mismatch_qsum
                                         >= stats_after.softclip_qsum + stats_after.mismatch_qsum))
                         && acceptRealignment(alignment, ref, dag_start_position, params, stats_after)) {
                    */

                    // we accept the new alignment if...
                    if (!was_mapped  // it wasn't mapped previously
                        // or if we have removed soft clips or mismatches (per quality) from the alignment
                        //|| ((stats_before.softclip_qsum >= stats_after.softclip_qsum
                        //     && stats_before.mismatch_qsum >= stats_after.mismatch_qsum)
                        || ((stats_before.softclip_qsum + stats_before.mismatch_qsum
                             >= stats_after.softclip_qsum + stats_after.mismatch_qsum)
                            // and if we have added gaps, we have added them to remove mismatches or softclips
                            && (stats_before.gaps >= stats_after.gaps // accept any time we reduce gaps while not increasing softclips/mismatches
                                || (stats_before.gaps < stats_after.gaps // and allow gap increases when they improve the alignment
                                    && (stats_before.softclip_qsum 
                                        + stats_before.mismatch_qsum
                                        >
                                        stats_after.softclip_qsum
                                        + stats_after.mismatch_qsum))))
                            // and the alignment must not have more than the acceptable number of gaps, softclips, or mismatches
                            // as provided in input parameters
                        && acceptRealignment(alignment, ref, dag_start_position, params, stats_after)) {

                        // keep the alignment
                        // TODO require threshold of softclips to keep alignment (or count of gaps, mismatches,...)
                        if (params.debug) {
                            cerr << "realigned " << alignment.Name << " to graph, which it maps to with "
                                 << stats_after.mismatch_qsum << "q in mismatches and "
                                 << stats_after.softclip_qsum << "q in soft clips" << endl;
                        }
                        ++total_improved;
                        has_realigned = true;
                    } else {
                        // reset to old version of alignment
                        if (params.debug) {
                            cerr << "failed realignment of " << alignment.Name << " to graph, which it maps to with: " 
                                 << stats_after.mismatch_qsum << "q in mismatches " << "(vs " << stats_before.mismatch_qsum << "q before), and "
                                 << stats_after.softclip_qsum << "q in soft clips " << "(vs " << stats_before.softclip_qsum << "q before) " << endl;
                        }
                        has_realigned = false;
                        alignment = originalAlignment;
                    }
                }
                //} // try block
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

        // ensure correct order if alignments move
        long int maxOutputPos = initialAlignmentPosition - dag_window_size;
        // if we switched sequences we need to flush out all the reads from the previous one
        string lastSeqname = currentSeqname;
        if (seqname != currentSeqname) {
            // so the max output position is set past the end of the last chromosome
            if (!currentSeqname.empty()) {
                maxOutputPos = reference.sequenceLength(currentSeqname) + dag_window_size;
            }
            currentSeqname = seqname;
        }

        if (!params.dry_run) {
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
            if (!params.only_realigned || has_realigned) {
                alignmentSortQueue[alignment.Position].push_back(alignment);
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

    gssw_graph_destroy(graph);
    free(nt_table);
	free(mat);

    reader.Close();
    writer.Close();

    if (params.debug) {
        cerr << "total reads:\t" << total_reads << endl;
        cerr << "realigned:\t" << total_realigned << endl;
        cerr << "improved:\t" << total_improved << endl;
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
