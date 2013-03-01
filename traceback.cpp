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

bt master_backtrack(sn* node, mbt &trace_report) {
    // Declare CIGAR string

    string cigar;
    string backstr = "";
    vector<string> node_list;
	
    // Recover starting matrix coordinates from the highest scoring entry
    int x = node->top_score.x;
    int y = node->top_score.y;
	
    //Go into the recursion with max score coordinates, empty cigar & empty node list.
    vector<bt> trace;
    bt result = backtrack(node, x, y, trace, backstr, node_list);
	
    // cout << "cigar: " << result.backstr << endl;
	
    x = result.x;
    y = result.y;
    cigar = result.backstr;

    reverse(node_list.begin(), node_list.end());
    reverse(trace.begin(), trace.end());

    vector<string>::iterator n = node_list.begin();
    vector<bt>::iterator t = trace.begin();

    stringstream cigarss;
    for ( ; t != trace.end() && n != node_list.end(); ++t, ++n) {
	reverse(t->backstr.begin(), t->backstr.end());
	cigarss << t->x << "," << t->y << ":" << *n << ":" << t->backstr;
	if ((t+1) != trace.end()) cigarss << "|";
    }
	
    // py: return {'node_id':node_id, 'x':x, 'y':y, 'cigar':cigar, 'nodes':nodes}
    trace_report.x = x;
    trace_report.y = y;
    trace_report.cigar = cigarss.str();
	
    trace_report.node_list = node_list;
    trace_report.node_name = node_list.back();
	
    return result;
}


/* 
 * The below does more than get a simple backtrace:  
 * It tries to recover all the aligned nodes - there is some rationale
 * to do this;  but it might not be very strong.  Might add / improve later
 * 
 * Will do final reversal in the end for C++ reasons
 * Might be a good idea to add a BAM-flattener (can also do in Python)
 */


bt backtrack(sn* node, int x, int y, vector<bt>& trace, string& backstr, vector<string> &node_list) {
    /* (1) Recover score from the node if possible, check for index error
     * (2) If score is 0, return existing cigar, node, position
     * (3) Otherwise,  add the new string, and call the function again!
     * (4) Do it anyways, but if the node changes, append node name, and change x coordinate.
     *
     * TODO: Update to Reflect Gaps
     * Long Term Todo: See how many recursions break it!
     */
	
    // Declare and initialize backtrack data Structure
    // When inside the bottom of recursing stack, this also reassigns final results
    bt backtrace;

    backtrace.x = x;
    backtrace.y = y;
	
    int score;
    char arrow;
	
    // Look for the score in the given coordinates. Try&Catch included for debugging. 
    try {
	score = node->matrix[y][x].score;
    }
    catch (exception& e) {
	cout<<"Standard Exception: "<<e.what()<<endl;
	cout<<"INDEX ERROR: "<<node->name<<x<<y<<score<<backstr<<endl;
    }
	
	
    /* Read the Arrow Matrix & Iterate over the Backtrace */
    if (score == 0) {
	// py: backstr = node.name + ':' + backstr;
	node_list.push_back(node->name);
	backtrace.backstr = backstr;
	trace.push_back(backtrace);
	return backtrace;

    } else {
	arrow = node->matrix[y][x].arrow;	  // TODO: check about pointers  ??
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
	    cout<<"BackTrace Error: Unknown Type";     // add proper error checking
	}
		
	sn* new_node = node->matrix[y][x].parent;

	if (new_node->name != node->name) {   // might be better way to cmp
	    x = new_node->seq_len;
	    node_list.push_back(node->name);
	    backtrace.backstr = backstr;
	    trace.push_back(backtrace);
	    backstr.clear();
	}
		
	return backtrack(new_node, x, y, trace, backstr, node_list);
    }
}

