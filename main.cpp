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

using namespace std;
using namespace vcf;

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
	

int main (int argc, const char * argv[])		// For the Stand Alone version
{

    
    Parameters params(argc, (char **)argv);
    
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

    vcffile.open(params.vcf_file);
    Variant var(vcffile);
    
    vcffile.setRegion(params.target);
    while (vcffile.getNextVariant(var)) {
	if (var.position + var.ref.length() <= target.stopPos) {
	    variants.push_back(var);
	}
    }

    if (variants.empty()) {
	cerr << "no variants" << endl;
    }

    long offset = target.startPos;

    // Declare the target DAG to align against.
    vector<sn*> nlist;

    //origIndel(nlist);

    constructDAG(nlist, targetSequence, variants, offset);

    if (params.display_dag) {
	cout << "DAG generated from input variants:" << endl;
	//displayDAG(nlist.back());
	for (vector<sn*>::iterator n = nlist.begin(); n != nlist.end(); ++n) {
	    cout << *n << endl;
	}
	cout << endl;
    }

    if (params.useFile == false) {
    
        string read = params.read_input;

	sn* result_F;
	sn* result_R;
	int score_F=0;
	int score_R=0;
	mbt trace_report_F;
	mbt trace_report_R;
	bt backtrace_F;
	bt backtrace_R;
	string strand = "+";
	
        result_F = gsw(read, nlist,
		       params.match, params.mism, params.gap);
	score_F = result_F->top_score.score;
	backtrace_F = master_backtrack(result_F, trace_report_F);
	vector<sn*>& nodes_F = trace_report_F.node_list;

	if (params.display_backtrace) {
	    cout << "==== forward alignment ====" << endl;
	    for (vector<sn*>::iterator n = nodes_F.begin();
		 n != nodes_F.end(); ++n) {
		displayAlignment(*n);
	    }
	} else if (params.display_all_nodes) {
	    cout << "==== forward alignment ====" << endl;
	    for (vector<sn*>::iterator n = nlist.begin();
		 n != nlist.end(); ++n) {
		displayAlignment(*n);
	    }
	}

	// check if the reverse complement provides a better alignment
	if (params.alignReverse) {
	    result_R = gsw(reverseComplement(read), nlist,
			   params.match, params.mism, params.gap);
	    score_R = result_R->top_score.score;
	    backtrace_R = master_backtrack(result_R, trace_report_R);
	    vector<sn*>& nodes_R = trace_report_R.node_list;
	    if (params.display_backtrace) {
		cout << "==== reverse alignment ====" << endl;
		for (vector<sn*>::iterator n = nodes_R.begin();
		     n != nodes_R.end(); ++n) {
		    displayAlignment(*n);
		}
	    } else if (params.display_all_nodes) {
		cout << "==== reverse alignment ====" << endl;
		for (vector<sn*>::iterator n = nlist.begin();
		     n != nlist.end(); ++n) {
		    displayAlignment(*n);
		}
	    }
	}

        //displayNode(result);
        //displayAlignment(result);
        //displayAlignment(nlist[0]);

	bt* backtrace; // best backtrace
	mbt* trace_report;
	int score;
	if (score_F > score_R) {
	    backtrace = &backtrace_F;
	    trace_report = &trace_report_F;
	    score = score_F;
	    strand = "+";
	} else {
	    backtrace = &backtrace_R;
	    trace_report = &trace_report_R;
	    score = score_R;
	    strand = "-";
	}

        //cout << "x: " << trace_report.x << " y: " << trace_report.y << endl;
        cout << score << " " << strand
	     << " seq:" << trace_report->x << " read:" << trace_report->y
	     << " " << trace_report->cigar << endl;
        
    } else {
        
        // Declere a read object from a fastq file   -> WHAT's the Best Place to do this  <-
        fastq_entry read_q;
        
        // Open FastQ File
        // TODO: use the parameter
        ifstream filehandle (params.fastq_file.c_str());		// open fastq file

        // Main input file loop
        if (filehandle.is_open()) {
            while (filehandle.good()) {
                
                // Load up the read with the four lines from the FastQ file
                read_q = getNextRead(filehandle);
                cout << read_q.readname << endl;
                
                // Reverse Complement the sequence information of the read.
                read_q.rc_sequence = reverseComplement(read_q.sequence);
                
                
                /* --- alignment -- */
                
                cout << "Root node of the DAG:\t"<< nlist.front()->name << endl
                << "Sequence:\t" << read << endl;
                
                // Check if the nodes get cleaned up.
                sn* result = gsw(read_q.sequence, nlist,
				 params.match, params.mism, params.gap);
                
                cout << "End of Alignment:\t" << result->name << "\t"
                << "top score:\t"<<result->top_score.score<<endl;
                
                
                displayNode(result);
                displayAlignment(result);
                //displayAlignment(nlist[0]);
                
                mbt trace_report;
                master_backtrack(result, trace_report);
                
                cout << "x: " << trace_report.x << " y: " << trace_report.y << endl;
                cout << trace_report.cigar << endl;
                
            }
            filehandle.close();
        } else {
            cout << "Unable to open read file.\n";
           }
    
		
	}
	


	//write_file();
	return 0;	
}
