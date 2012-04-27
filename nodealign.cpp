/*  
 * Python Notice:
 * Created on Apr 15, 2011
 * @author: kural
 *
 * Functions here contain the graph alignment algorithm,
 * implemented in the block-matrix approach.
 *
 * nodealign.cpp
 * glia
 *
 * Created by Deniz Kural on 6/28/11.
 * Copyright 2011 Deniz Kural. All rights reserved.
 *
 */

#include "nodealign.h"

using namespace std;


/* 
 * Declare Data Structures
 * TODO: Some to graduate into classes?
 */


// compare parent nodes by score.
struct cmp_parent_nodes {
	int y_i;
    bool operator () (const sn* node_x, const sn* node_y) {
		return node_x->score_matrix[y_i][node_x->seq_len] > node_y->score_matrix[y_i][node_x->seq_len];        
    }
};



/* Set Alignment Scores */
// todo: add gap-extend, consider floats, consider user input */

const char match = 10;
const char mism = -10;
const char gap = -10;

// const char gap_open = -10
// const char gap_extend = -2


/*
 * INPUTS:  read sequence, read length, node
 * OUTPUTS: alignment
 *
 * TODO:    Accept node pointer instead?
 * TODO:    Pass outputs as inputs-by-reference
 * TODO:    Check consistency of char vs int
 *
 */
int StringNodeAlign(string read, int read_length, sn &node) {
	/* Main Alignment Algorithm for String and Node */
	
	// Initialize variables
	char dia_score;
	char dia;
	
	sn* top_parent;					// pointer or not)?

	
	// Initialize Top Score before the alignment starts
	ts top_score;
	top_score.score = 0;
	top_score.x = 0;
	top_score.y = 0;

	
	// declare the score matrix
	vector<vector<int> > score_matrix (read_length+1, vector<int>(node.seq_len+1, 0));

	// declare the backtrace matrix
	vector<vector<char> > arrow_matrix (read_length+1, vector<char>(node.seq_len+1, 0));

	// declare the parent matrix
	vector<vector<sn*> > parent_matrix (read_length+1, vector<sn*>(node.seq_len+1, &node));
	

	/* Inherit The Top values for the leftmost utility column */
	
	// Seed the Service Left Column to 0 or parents
	if (node.p5.empty()) {
		for (int i = 0; i < read_length +1; i++) {
			parent_matrix[i][0] = &node; }                // redundant, (replace with what?)
	} else {
		for (int i = 1; i < read_length + 1; i++) {
			// i = 1 because coordinate [0,0] gets seeded with the top row above
		
			// THIS LINE BELOW IS NOT PORTED!!!
			//top_parent = sorted([ [p.mscore[i][-1], p.mscore[i][-1][0] ] for p in node.p5] key=itemgetter(1), reverse=True)[0][0]
			
			//cmp_parent_nodes.y_i = i;
			cmp_parent_nodes my_compare;
			my_compare.y_i = i;
			sort(node.p5.begin(), node.p5.end(), my_compare);
			
			top_parent = node.p5[0];
			
			// DEBUG PRINT
			//displayAlignment(top_parent);
			//cout << "top parent: " << top_parent->name;
			//cout << "score at y:" << i << "x:" << top_parent->seq_len << " -- "; 
			//cout << top_parent->score_matrix[i][top_parent->seq_len]<<endl;
			
			score_matrix[i][0] = top_parent->score_matrix[i][top_parent->seq_len];  // convert to pointers?
			parent_matrix[i][0] = top_parent;
			arrow_matrix[i][0] = top_parent->arrow_matrix[i][top_parent->seq_len];  // convert to pointers?
			
			
			//Update Top Score
			if (score_matrix[i][0] > top_score.score) {
				top_score.score = score_matrix[i][0];
				top_score.y = i;
				top_score.x = 0;
			}
		}
	}

	
	/* Do the Alignment */
	for (int y = 1; y < read_length + 1; y++) {
		for (int x = 1; x < node.seq_len + 1; x++) {
				
		/* Figure out if the base is a match or a mismatch
		 Note that the indexes x-1 and y-1 refer to the sequences and not to the matrix
		 The first letters indexed at 0 correspond to  (1,1)  on the matrix */
			
			// debug print
			// if (node.id == "n1") {
			// cout <<node.sequence[x-1]<<x-1<<read[y-1]<<y-1;
			// }
			
			if (node.sequence[x-1] == read[y-1]) {
				dia_score = match;
				dia = 'm';
			} else {
				dia_score = mism;
				dia = 'x';
			}

			/* Calculate the potential scores */
				
			// debug prints
			// cout node.mscore[y][x];
			// cout node.mscore[y-1][x-1];
			// cout node.mscore[y][x-1];

			char dia_total = score_matrix[y-1][x-1] + dia_score;
			char up_total = score_matrix[y-1][x] + gap;
			char side_total = score_matrix[y][x-1] + gap;
			
			
			// Compare scores & update matrices (score + arrow)
			// todo: Add gap_extend
			// todo: Separate out the update function to a small function..
			
			if (dia_total >= up_total) {
				if (dia_total >= side_total) {
					if (dia_total < 0)
						dia_total = 0;
					arrow_matrix[y][x] = dia;
					score_matrix[y][x] = dia_total;
				} else {
					if (side_total < 0)
						side_total = 0;
					arrow_matrix[y][x] = 's';
					score_matrix[y][x] = side_total;
				}
			} else if (up_total >= side_total) {
				if (up_total < 0)
					up_total = 0;
				arrow_matrix[y][x] = 'u';
				score_matrix[y][x] = up_total;		
			} else {
				if (side_total < 0)
					side_total = 0;
				arrow_matrix[y][x] = 's';
				score_matrix[y][x] = side_total;
			}

			/* Assign the matrix
			 the format:  score, arrow, node.id
			 py: node.mscore[y][x] = [score[1], score[0], node.id, node]  */

			/* Update Top Score */
			if (score_matrix[y][x] > top_score.score) {
				top_score.score = score_matrix[y][x];
				top_score.y = y;
				top_score.x = x;
			}
		
		}	// close the row
	}	// close the column
			
	node.top_score = top_score;
	node.score_matrix = score_matrix;
	node.parent_matrix = parent_matrix;
	node.arrow_matrix = arrow_matrix;
	return 0;
}



	
	
	
