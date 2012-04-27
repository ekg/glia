#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include <google/sparse_hash_map>
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


using namespace std;


string reverseComplement(string read) {
	return "hello!";
}

// TODO: Move to GHASH
int hashfasta(string fasta_file_name, int hashsize, vector<fasta_entry> &ref_genome) {	

		
	// map<string, vector<int> > contig_hash;
	// google::sparse_hash_map<string, vector <int>, tr1::hash<string> > contig_hash;
	/* Disabled after moving the hash to the fasta_entry struct
	google::sparse_hash_map<string, vector <int> > contig_hash;
	contig_hash.set_deleted_key("");
	*/	

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
	

int main (int argc, char * const argv[])
{
	
	/* ___ PART 1:  Hash the Genome ________________________________  */
    
    // The File Path of the Reference Genome (JSON?)
	string fasta_file_name = "/Users/kural/Downloads/chr20.fa";
	
	// The Reference Genome (a vector of fasta entries)
	vector<fasta_entry> ref_genome;	
	
	// kmer hash build size.. very important (JSON?)
	int hashsize = 25;

	/* Loads up the genome & hashes it. Might carve out in the future as 
	   a build step, after a binary serialization is done */
	hashfasta(fasta_file_name, hashsize, ref_genome);


	
	/* ___ PART 2:  Read FASTQ and Align Reads _____________________ */

	
	//load_fastq_file();
	//write_file();
	string read = "CTTCTTCTTCTTCTTCTTCTTCTTCCTTCTTCTTCTTCTTCTTCTTCTTC";
	
	/* todo: 
	   - simulate some reads (Single End)
	   - create fastq read objects to hold a) the reverse complement data b) the seeds.
	   - get seeds
	   - once rest of it is done, align the read
	 */
	
	
	vector<fasta_entry>::iterator t;
	for (t = ref_genome.begin(); t != ref_genome.end(); ++t) {

		lookupRead(read, t->kmer_hash, hashsize);
	}
		

	/* --- alignment -- */
	
	// string read = "ATCGAA";

	
	vector<sn*> nlist;
	origIndel(nlist);

	
	cout << "in the main:"<< endl<<nlist.front()->name << endl;
	sn* result = gsw(read, nlist);
	
	cout<<result->name<<", top score: "<<result->top_score.score<<endl;
	
	displayNode(result);
	displayAlignment(result);
	displayAlignment(nlist[0]);
	
	mbt trace_report;
	master_backtrack(result, trace_report);
	
	cout << "x: " << trace_report.x << " y: " << trace_report.y << endl;
	cout << trace_report.cigar << endl;
	
	return 0;
	
}

