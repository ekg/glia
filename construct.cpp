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

int constructDAG(vector<sn*> &nlist, string &targetSequence, 
		 vector<Variants> &variants, long offset) {


  long  current_pos;
  long prevs_pos = targetSequence.size();
  string right_seq = "";
  string left_seq = "":
  
for(vector<Variants>::reverse_iterator rit = variants.rbegin(); 
      rit != variants.rend(); ++rit) {

    // Construct Right-Node    
    current_position = rit->position - offset;
    
    // Var Type changes this
    right_seq = targetSequence.substr(current_pos, (prev_pos - current_pos));
    left_seq = targetSequence.substr(0,current_pos);

    
    // Construct Right Node
    sn* right_node;
    
    sn* na1;
    sn* na2;

    sn* left_node;

    right_node = new sn;

    na1 = new sn;
    na2 = new sn;

    left_node = new sn; 

    // Fill out Right Node
    right_node->name.append("ref");
    right_node->name.append(to_string(current_pos));

    right_node->sequence = right_sequence;
    right_node->seq_len = right_node->sequence.length();
    right_node->depth = -1;
    
    // Fill out Allele Nodes

    // Fill out Left Nodes


    // Connect Nodes
    
    right_node->p5.push_back(na1);
    right_node->p5.push_back(na1);
    
    na1->p3.push_back(right_node);
    na2->p3.push_back(right_node);
    
    na1->p5.push_back(left_node);
    na1->p5.push_back(left_node);

    left_node->p3.push_back(na1);
    left_node->p3.push_back(na2);
    
    nlist.push_back(na1);
    nlist.push_back(na2);
    nlist.push_back(right_node);
    

    }
  
  last-node

    
    
    
 

}


