//
//  construct.cpp
//  
//
//  Created by Deniz Kural on 2/19/13.
//
//


#include "construct.h"

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
	int depth;				// leave uninitialized?
}; 
 */





int constructDAG(vector<sn*> &nlist, 
