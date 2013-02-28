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
    
    cout << "read: " << params.read_input << endl;
    //cout << "fastq file:" << params.fastq_file << endl;
    cout << "fasta reference:" << params.fasta_reference << endl;

    cout << "vcf file " << params.vcf_file << endl;
    cout << "target " << params.target << endl;


    // get variants in target
    vector<Variant> variants;
    VariantCallFile vcffile;

    vcffile.open(params.vcf_file);
    Variant var(vcffile);
    
    vcffile.setRegion(params.target);
    while (vcffile.getNextVariant(var)) {
	variants.push_back(var);
    }

    // get sequence of target
    FastaReference reference;
    reference.open(params.fasta_reference);
    FastaRegion target(params.target);

    string targetSequence = reference.getSubSequence(target.startSeq,
						     target.startPos - 1,
						     target.length());    


    long offset = target.startPos;

    // Declare the target DAG to align against.
    vector<sn*> nlist;

    constructDAG(nlist, targetSequence, variants, offset);
    //origIndel(nlist);
    //json_example(nlist);

    
    if (params.useFile == false) {
    
        string read = params.read_input;
    
        sn* result = gsw(read, nlist);
        
        cout << "End of Alignment:\t" << result->name << "\t"
        << "top score:\t"<<result->top_score.score<<endl;
        
        
        displayNode(result);
        displayAlignment(result);
        //displayAlignment(nlist[0]);
        
        mbt trace_report;
        master_backtrack(result, trace_report);
        
        cout << "x: " << trace_report.x << " y: " << trace_report.y << endl;
        cout << trace_report.cigar << endl;
        
        
    } else {
        
        
        /* Ancient 
         //string read = "CTTCTTCTTCTTCTTCTTCTTCTTCCTTCTTCTTCTTCTTCTTCTTCTTC";
         //string read = "ATCGAA";
         */
        
        // Declere a read object from a fastq file   -> WHAT's the Best Place to do this  <-
        fastq_entry read_q;
        
        // Open FastQ File
        // TODO: use the parameter
        ifstream filehandle ("test.fastq");		// open fastq file
        
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
                sn* result = gsw(read_q.sequence, nlist);
                
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






int realmain (int argc, char * const argv[])
{
	
	/* ___ PART 1:  Hash the Genome "Build Step" __________________  */
	
    
    // The File Path of the Reference Genome (JSON?)
	string fasta_file_name = "/Users/kural/Downloads/chr20.fa";
	
	// The Reference Genome (a vector of fasta entries)
	vector<fasta_entry> ref_genome;	
	
	// kmer hash build size.. very important
	int hashsize = 25;

	/* Loads up the genome & hashes it. Might carve out in the future as 
	   a build step, after a binary serialization is done 
	   (i.e. can write out a "processed genome build"  */
	hashfasta(fasta_file_name, hashsize, ref_genome);


	
	/* ___ PART 2:  Read FASTQ and Align Reads _____________________ */

	
	// testing
	//string read = "CTTCTTCTTCTTCTTCTTCTTCTTCCTTCTTCTTCTTCTTCTTCTTCTTC";
	string read = "ATCGAA";
	
	// Declere a read object from a fastq file
	fastq_entry read_q;					
	
	// Open FastQ File
	ifstream filehandle ("test.fastq");		// open fastq file
	
	// Main input file loop
	if (filehandle.is_open()) 
	{
		while (filehandle.good()) 
		{

			// Load up the read with the four lines from the FastQ file
			read_q = getNextRead(filehandle);
			cout << read_q.readname << endl;
			
			// Reverse Complement the sequence information of the read.
			read_q.rc_sequence = reverseComplement(read_q.sequence);

			
			
			// For each chromosome,  get clusters ...
			vector<fasta_entry>::iterator t;
			for (t = ref_genome.begin(); t != ref_genome.end(); ++t) {
				
				// vector<pair<int, int> > clusters;

				lookupRead(read, t->kmer_hash, hashsize);
			}
			
			
			/* --- alignment -- */
			
			//string read = "ATCGAA";
			
			
			// Declare the target DAG to align against. 
			vector<sn*> nlist;
			origIndel(nlist);

			
			
			cout << "in the main:\t"<< nlist.front()->name << "\tread:\t" << read << endl;
			sn* result = gsw(read, nlist);
			
			cout<<result->name<<", top score:\t"<<result->top_score.score<<endl;
			
			displayNode(result);
			displayAlignment(result);
			displayAlignment(nlist[0]);
			
			mbt trace_report;
			master_backtrack(result, trace_report);
			
			cout << "x: " << trace_report.x << " y: " << trace_report.y << endl;
			cout << trace_report.cigar << endl;
			
		}
		filehandle.close();
		
	}
	else cout << "Unable to open read file.\n"; 
	
	//write_file();
	
	return 0;
}




