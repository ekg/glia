/*
 *  examples.cpp
 *  glia
 *
 *  Created by Deniz Kural on 12/27/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "examples.h"

using namespace std;


// sn := string-node data structure
/*
struct sn {
	std::string sequence;
	std::string name;
	std::vector<sn*> p5;			// now this will be sn pointers
	std::vector<sn*> p3;			// now this will be sn pointers
	int seq_len;
	ts top_score;
	std::vector<std::vector<int> > score_matrix;
	std::vector<std::vector<char> > arrow_matrix;
	std::vector<std::vector<sn*> > parent_matrix;
	int depth;						// leave uninitialized?
}; 
 */


// 	string read = "ATCGAA";

/* Create an Indel Example */
int origIndel(vector<sn*> &nlist) {
	
	// Initialize Nodes
	sn* n1;
	sn* n2;
	sn* n3;
	n1 = new sn;
	n2 = new sn;
	n3 = new sn;
	
	n1->name = "n1";
	n1->sequence = "TTGGAT";
	n1->seq_len = n1->sequence.length();
	n1->depth = -1;
	
	
	//This is the inserted sequence...
	n2->name  = "n2";
	n2->sequence = "CGAATT";
	n2->seq_len = n2->sequence.length();
	n2->depth = -1;

	n3->name = "n3";
	n3->sequence = "ATGGG";
	n3->seq_len = n3->sequence.length();
	n3->depth = -1;
	
	//Connect Nodes
	n1->p3.push_back(n2);
	n1->p3.push_back(n3);

	n2->p3.push_back(n3);

	n2->p5.push_back(n1);
	n3->p5.push_back(n1);
	n3->p5.push_back(n2);

	//vector<sn*> nlist; 
	
	cout << "the instance in the function: " << n1->name << n1->sequence << endl;
	
	nlist.push_back(n1);
	nlist.push_back(n2);
	nlist.push_back(n3);
	
	cout << "the vector in the function: " << nlist[1]->name << endl;
	
	return 0;
}

