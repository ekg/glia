/*
 *  fastq.cpp
 *  glia
 *
 *  Created by Deniz Kural on 12/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "fastq.h"

using namespace std;

/* File writer prototype 
   Haven't checked in a while */
int write_bam_file () 
{
	ofstream outfile;
	outfile.open("test.out");
	outfile << "Writing to the file.\n";
	outfile.close();
	return 0;
}



// Return a Read Quad from a fastq file
fastq_entry getNextRead(ifstream& filehandle)
{
	fastq_entry readquad;
	getline (filehandle,readquad.readname);
	// cout << readquad.readname << endl;
	getline (filehandle,readquad.sequence);
	getline (filehandle,readquad.description);
	getline (filehandle,readquad.quality);
	return readquad;
}


/* Function to Load a FASTQ File
 See prototype loader on dump.cpp */

int load_fastq_file () {
	// Declare Variables
	fastq_entry readquad;					// declare fastq read
	ifstream filehandle ("test.fastq");		// open fastq file
	
	// Main input file loop
	if (filehandle.is_open()) 
	{
		while (filehandle.good()) 
		{
			readquad = getNextRead(filehandle);
			cout << readquad.readname << endl;
		}
		filehandle.close();
		
	}
	else cout << "Unable to open read file.\n"; 
	return 0;
}


/* Desc: Load a Fasta File to memory
 * Inputs: FASTA file name, Vector of fasta entries
 * Output: Fill the vector with entries contig name, sequence
 */

int load_fasta_file (string fasta_file_name, vector<fasta_entry> &refs) {
	ifstream filehandle (fasta_file_name.c_str());						// open fasta file	
    // ifstream filehandle ("/Users/kural/Documents/short1.fasta");		// test
	// vector<fasta_entry> refs;
	fasta_entry fastarecord;										
	string readline;

	// Main input file loop
	if (filehandle.is_open()) {
		while (filehandle.good()) {
			
			getline(filehandle, readline);
			/* finish current fastarecord, start new one & give name
			 * build exception for beginning
			 */
			if (readline[0] == '>') {							// if beginning of new entry
				if (fastarecord.name.empty()) {				
					fastarecord.name = readline;
					fastarecord.kmer_hash.set_deleted_key("");
				} else {
					refs.push_back(fastarecord);
					fastarecord.name = readline;
					fastarecord.sequence = "";
				}
			} else {
				fastarecord.sequence.append(readline);
			}
		}
		
		refs.push_back(fastarecord);
		filehandle.close();
		
		/*
		vector<fasta_entry>::iterator t;									
		
		for (t = refs.begin(); t != refs.end(); ++t) {
			cout << t->name << endl;
			cout << t->sequence << endl;
		}
		*/
	}
	else cout << "Unable to open read file.\n"; 
	return 0;
}


