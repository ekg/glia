/*
 *  ghash.h
 *  glia
 *
 *  Created by Deniz Kural on 1/4/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GHASH_H
#define GHASH_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>

#include <google/sparse_hash_map>
//#include <tr1/functional>

#include "fastq.h"

struct sh_eqstr {
	bool operator()(const char* s1, const char* s2) const {
		return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
	}
};


// int hashcontig(fasta_entry &fasta_record, std::map<std::string, std::vector<int> > &ghash);
// int hashcontig(fasta_entry &fasta_record, google::sparse_hash_map<std::string, std::vector<int>, std::tr1::hash<std::string> > &ghash);
int hashcontig(fasta_entry &fasta_record, google::sparse_hash_map<std::string, std::vector<int> > &ghash, int hashsize);


int sortContigHash(google::sparse_hash_map<std::string, std::vector <int> > &ghash);



#endif