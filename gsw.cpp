/*
 *  gsw.cpp
 *  glia
 *
 *  Created by Deniz Kural.
 *  Copyright 2011 Deniz Kural. All rights reserved.
 *
 */


#include "gsw.h"

using namespace std;


// compare parents of a node by depth
struct cmp_gd_depths {
    bool operator () (sn* node_x, sn* node_y) {				// warning: got rid of consts!
		return getDepth(node_x) > getDepth(node_y);        
    }
} my_gd_compare;


// constant version of compare parents of a node by depth,
// assumes that depth already exists correctly
struct cmp_parent_depths {
    bool operator () (const sn* node_x, const sn* node_y) {
		return node_x->depth > node_y->depth;        
    }
} my_depth_compare;


int getDepth(sn* node) {
    /* Return the topological order of a node in a DAG
       - Check if it has depth, if so don't repeat the backtrack
       - Else, check if it has a parent node, if so recurse
       - Otherwise declare the node to be a head node 
    */
    
    int d;										// optimization:  get rid of this
    if (node->depth == -1) {
        if (node->p5.empty() == false) {
            //d = 1 + max( [getDepth(p) for p in node.p5] )
            for (vector<sn*>::iterator n = node->p5.begin() ;n != node->p5.end(); ++n) {
                getDepth(*n);
            }
            sort(node->p5.begin(), node->p5.end(), my_gd_compare);
            d = node->p5[0]->depth + 1;			// check to see if this is same as .front()
            node->depth = d;
            return d;
        } else {
            return 1;
        }
    } else {
        return node->depth;
    }
}	


sn* sequenceDagAlign(string sequence, vector<sn*> nlist, int maxdepth,
		     const int match, const int mism, const int gap) {

    int top_score = 0;
    sn* top_node = NULL;

    for (int i = 1; i < maxdepth + 1; ++i) {						// why 1, +1? 
        vector<sn*>::iterator t;									// good place for vector iteration?
        for (t=nlist.begin(); t!=nlist.end(); ++t) {
            if ((*t)->depth == i) {
                StringNodeAlign(sequence, sequence.length(), **t,
                                match, mism, gap);
                if ((*t)->top_score.score > top_score) {
                    top_node = (*t);
                    top_score = (*t)->top_score.score;
                }
            }
        }
    }
    // Note that top score is not returned, but just used internally. 
    return top_node;
}


sn* gsw(string read, vector<sn*> nlist,
	const int match, const int mism, const int gap) {
    /* Create a local POSET and record the poset layer ordering on the Node object.
     * This might not be necessary in the full framework - for example when a reference DAG is already pre-sorted
     * and comes with appropriate depth. This assumes a local unsorted framework. 
     */
    int maxdepth;

    vector<sn*>::iterator t;
    for (t=nlist.begin(); t!=nlist.end(); ++t) {
        (*t)->depth = getDepth(*t);
        // as long as this isn't an empty node (anchor)
        // initialize the scoring matrix for this specific read
        (*t)->initScore(read.size());
    }
		
    sort(nlist.begin(), nlist.end(), my_depth_compare);
    maxdepth = -1;
    for (t=nlist.begin(); t!=nlist.end(); ++t) {
        int d = (*nlist.begin())->depth;
        if (d > maxdepth) maxdepth = d;
    }

    /* Submits the read & each node in the correct sequence to the block-GSW algorithm:
       - Picks the root node(s),  proceeds according to topological sort so child nodes inherit parents
       - Performs Alignment in random order for all nodes of equal topological order
       - Repeat
	 
       The JSON component for reporting is a bit of a hack and should be fixed
    */


    //Calls the above function to actually perform alignment
    sn* result = sequenceDagAlign(read, nlist, maxdepth,
                                  match, mism, gap);
	
    //alignment = nodealign.master_backtrack(result);
    

    return result;
}


