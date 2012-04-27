/*
 *  ghash.cpp
 *  glia
 *
 *  Created by Deniz Kural.
 *  Copyright 2010-2012 Deniz Kural. All rights reserved.
 *
 */

#include "ghash.h"
using namespace std;

/* 
 Desc: Accept a fasta record,  output a hash of the record.
 Input: Fasta Entry,  Hash Structure, Hash Size
 Output: Modified Hash Structure
 */

int hashcontig(fasta_entry &fasta_record, google::sparse_hash_map<string, vector <int> >  &ghash, int hashsize) {
	
	bool test_mode = true;
	
	// todo: Change this later / make it a variable
	// todo: Make sure that fasta string length > hash size
	string kmer;

	// map<string, vector<int> > ghash;
	
	/* Transform the Reference to Upper Case */
	cout << "transforming reference into upper case..." << endl;
	transform(fasta_record.sequence.begin(), fasta_record.sequence.end(), fasta_record.sequence.begin(), ::toupper);
	cout << "transform Complete, size: "<< fasta_record.sequence.size() << endl;
	
	/* For Test Purposes We limit hash map size */
	int i_limit; 
	
	if (test_mode == true) {
		i_limit = 1000000;
	} else {
		i_limit = fasta_record.sequence.size() - hashsize + 1;
	}
	/* Test Logic End */

			
	for (int i = 0; i < i_limit; ++i) {
		kmer = fasta_record.sequence.substr(i, hashsize);		// optimization: kmer not necessary [beware]
		//cout << kmer << " " << i << endl;
		
		ghash[kmer].push_back(i);
		//cout << "+ " << ghash[kmer].front() << " " << ghash[kmer].back() << endl;
	}
	
	
	cout << "hashcontig ghash size :" << ghash.size() << endl;

	return 0;
}


int sortContigHash(google::sparse_hash_map<string, vector <int> > &ghash) {

	// iterate over hash keys (kmers)
	google::sparse_hash_map<string, vector <int> >::iterator t;
	
	cout << "sortcontig ghash size :" << ghash.size() << endl;
	
	for (t = ghash.begin(); t != ghash.end(); t++) {
		// sort each vector in place
		sort((*t).second.begin(), (*t).second.end());
		
		// print kmer : # occurences : 1st occurance
		// cout << (*t).first << " : " << (*t).second.size() << " : " << (*t).second.front() << endl;
	}
	
	
	// ...
	// it should then be written/serialized somewhere
	return 0;
}
	



// deprecated

/* 
string::iterator t;									
for (t = fasta_record.sequence.begin(); t != fasta_record.sequence.end() - hashsize; ++t) {		
    cout<<*t<<"-";  
}
*/
