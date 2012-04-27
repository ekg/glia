/*
 *  fastq.h
 *  glia
 *
 *  Created by Deniz Kural on 12/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef FASTQ_H
#define FASTQ_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <google/sparse_hash_map>

// data structure for a fastq entry
struct fastq_entry {
	std::string readname;
	std::string sequence;
	std::string description;
	std::string quality;
	
	std::string complement;
	// add clusters ..?
	// add alignment .. ?
	
};

// data structure for a fastq entry
struct fasta_entry {
	std::string name;
	std::string sequence;
	// string seq_length;
	google::sparse_hash_map<std::string, std::vector <int> > kmer_hash;
};

// flesh out into bam writer
int write_bam_file ();

// reads fastq entries into a fastq structure
fastq_entry getNextRead(std::ifstream& filehandle);

// fastq loader function
int load_fastq_file (); 

// fasta loader function
int load_fasta_file (std::string fasta_file_name, std::vector<fasta_entry> &refs);

#endif

