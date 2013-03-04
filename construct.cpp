//
//  construct.cpp
//  
//
//  Created by Deniz Kural on 2/19/13.
//
//


#include "construct.h"
#include "convert.h"

using namespace std;
using namespace vcf;

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


int constructDAG(vector<sn*> &nlist, string &targetSequence, string& sequenceName,
		 vector<Variant> &variants, long offset) {


    long prev_pos = targetSequence.size();
    string p3_ref_seq;

    vector<sn*> pp3_var_nodes; // previous p3 nodes

    for(vector<Variant>::reverse_iterator rit = variants.rbegin(); 
	rit != variants.rend(); ++rit) {

	Variant& var = *rit;
	int current_pos = var.position - offset + var.ref.size();

	// Construct Right-Node    
	p3_ref_seq = targetSequence.substr(current_pos, prev_pos - current_pos);

	// new previous position is at the start of the variant
	prev_pos = var.position - offset;

	// Construct Right Node
	sn* p3_ref_node = new sn(
	    p3_ref_seq
	    ,
	    sequenceName
	    + ":"
	    + convert(var.position + var.ref.size())
	    + "-"
	    + convert(var.position
		      + var.ref.size()
		      + p3_ref_seq.size())
	    + ".ref.r");

	// connect to old p3 nodes
	for (vector<sn*>::iterator n = pp3_var_nodes.begin(); n != pp3_var_nodes.end(); ++n) {
	    p3_ref_node->p3.push_back(*n);
	    (*n)->p5.push_back(p3_ref_node);
	}
	pp3_var_nodes.clear();

	// construct the ref node
	sn* ref_node = new sn(
	    var.ref
	    ,
	    sequenceName
	    + ":"
	    + convert(var.position)
	    + "-"
	    + convert(var.position
		      + var.ref.size())
	    + ".ref.0");

	// stash p3_ and current ref nodes
	nlist.push_back(p3_ref_node);
	nlist.push_back(ref_node);

	// connect p3_ref <-> ref
	p3_ref_node->p5.push_back(ref_node);
	ref_node->p3.push_back(p3_ref_node);

	// store the current ref in the pp3 nodes for connection on next iteration
	pp3_var_nodes.push_back(ref_node);

	// Fill and connect Allele Nodes to p3_ref_node

	int i = 1;
	for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++i) {
	    sn* alt_node = new sn(
		*a
		,
		sequenceName
		+ ":"
		+ convert(var.position)
		+ "-"
		+ convert(var.position
			  + var.ref.size())
		+ ".alt." + convert(i));
	    // save in nlist
	    nlist.push_back(alt_node);
	    // retain for connection to ref p3_ref_node of next variant
	    pp3_var_nodes.push_back(alt_node);
	    // connect to current p3_ref_node
	    alt_node->p3.push_back(p3_ref_node);
	    // and connect the p5 of the p3_ref_node to the alt node
	    p3_ref_node->p5.push_back(alt_node);
	}

    }

    // last node construction and connection

    // Construct Right-Node    
    p3_ref_seq = targetSequence.substr(0, prev_pos);

    // Construct Right Node
    sn* p3_ref_node = new sn(
	p3_ref_seq
	,
	sequenceName
	+ ":"
	+ convert(offset)
	+ "-"
	+ convert(prev_pos + offset)
	+ ".ref.r");

    // connect to old p3 nodes
    for (vector<sn*>::iterator n = pp3_var_nodes.begin(); n != pp3_var_nodes.end(); ++n) {
	p3_ref_node->p3.push_back(*n);
	(*n)->p5.push_back(p3_ref_node);
    }
    pp3_var_nodes.clear();

    // save reference
    nlist.push_back(p3_ref_node);

    reverse(nlist.begin(), nlist.end());

}
