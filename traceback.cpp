/*
 *  traceback.cpp
 *  glia
 *
 *  Created by Deniz Kural.
 *  Copyright 2011 Deniz Kural. All rights reserved.
 *
 */


#include "traceback.h"

using namespace std;


/* Room for optimization after debugging:
 * (1) Get rid of int x, int y
 * (2) Get rid of intermediary cigar
 * (3) Get rid of intermediary x,y re-assignment
 */

int master_backtrack(sn* node, mbt &trace_report) {
	// Declare CIGAR string
	string backstr = "";
	vector<string> node_list;
	string cigar;
	
	
	// Recover starting matrix coordinates from the highest scoring entry
	int x = node->top_score.x;
	int y = node->top_score.y;
	
	//Go into the recursion
	bt result = backtrack(node, x, y, backstr, node_list);
	
	// cout << "cigar: " << result.backstr << endl;
	
	x = result.x;
	y = result.y;
	cigar = result.backstr;
	
	// py: return {'node_id':node_id, 'x':x, 'y':y, 'cigar':cigar, 'nodes':nodes}
	trace_report.x = x; trace_report.y = y; trace_report.cigar = cigar;
	trace_report.node_list = node_list;  trace_report.node_name = node_list.back();
	
	return 0;
}


/* 
 * The below does more than get a simple backtrace:  
 * It tries to recover all the aligned nodes - there is some rationale
 * to do this;  but it might not be very strong.  Might add / improve later
 * 
 * Will do final reversal in the end for C++ reasons
 * Might be a good idea to add a BAM-flattener
 */


bt backtrack(sn* node, int x, int y, string backstr, vector<string> &node_list) {
	/* (1) Recover score from the node if possible, check for index error
	 * (2) If score is 0, return existing cigar, node, position
	 * (3) Otherwise,  add the new string, and call the function again!
	 * (4) Do it anyways, but if the node changes, append node name, and change x coordinate.
	 *
	 * TODO: Update to Reflect Gaps
	 */
	
	// Declare and initialize backtrack data Structure
	bt backtrace;
	backtrace.x = x;
	backtrace.y = y;
	
	int score;
	char arrow;
	
	// Look for the score in the given coordinates. Try&Catch included for debugging. 
	try {
		score = node->score_matrix[y][x];
	}
	catch (exception& e) {
		cout<<"Standard Exception: "<<e.what()<<endl;
		cout<<"INDEX ERROR: "<<node->name<<x<<y<<score<<backstr<<endl;
	}
	
	
	/* Read the Arrow Matrix & Iterate over the Backtrace */
	if (score == 0) {
		// py: backstr = node.name + ':' + backstr;
		backstr.append(":");
		backstr.append(node->name);
		node_list.push_back(node->name);
		backtrace.backstr = backstr;
		
		return backtrace;
		
	} else {
		arrow = node->arrow_matrix[y][x];					// TODO: check about pointers
		if (arrow == 'm') {
			// py: backstr = "M" + backstr
			backstr.append("M");
			x = x - 1;
			y = y - 1;
		} else if (arrow == 'x') {
			// py: backstr = "X" + backstr
			backstr.append("X");
			x = x - 1;
			y = y - 1;
		} else if (arrow == 'u') {
			// py: backstr = 'D' + backstr
			backstr.append("D");
			y = y - 1;
		} else if (arrow == 's') {
			// py: backstr = 'I' + backstr
			backstr.append("I");
			x = x - 1;
		} else {
			cout<<"BackTrace Error: Unknown Type";           // add proper error checking
		}
		
		
		sn* new_node = node->parent_matrix[y][x];
		
		if (new_node->name != node->name) {							// might be better way to cmp
			x = new_node->seq_len;
			
			// py: backstr = '|' + node.name + ':' + backstr
			backstr.append(":");
			backstr.append(node->name);
			node_list.push_back(node->name);
			backstr.append("|");
		}
		
		return backtrack(new_node, x, y, backstr, node_list);
	}
}

