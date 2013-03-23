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

bt master_backtrack(sn* node, mbt &trace_report, string& read, string& qualities) {
    // Declare CIGAR string

    string gcigar;
    Cigar fcigar;
    string backstr = "";
    vector<sn*> node_list;

    trace_report.read = read; // copy for possible modification
    trace_report.qualities = qualities;
	
    // Recover starting matrix coordinates from the highest scoring entry
    int x = node->top_score.x;
    int y = node->top_score.y;
	
    //Go into the recursion with max score coordinates, empty cigar & empty node list.
    vector<bt> trace;

    bt graphresult = graphbacktrack(node, x, y, trace, gcigar, node_list);
    reverse(node_list.begin(), node_list.end());
    reverse(trace.begin(), trace.end());
    vector<bt>::iterator t = trace.begin();
    stringstream gcigarss;
    for ( ; t != trace.end() && t != trace.end(); ++t) {
        reverse(t->backstr.begin(), t->backstr.end());
        gcigarss << t->node->name << ":"
                << t->x << "," << t->y << ";"
                << t->backstr;
        if ((t+1) != trace.end()) gcigarss << "|";
    }
	
    trace.clear();
    node_list.clear();

    /*
    if (!node->isref) {
        cerr << 3 << endl;
        cerr << node->sequence << " size " << node->sequence.size() << " " << x << endl;
        if (node->sequence.size() > x) {
            trace_report.implied_p3  = node->sequence.substr(x-1);
        }
        cerr << 4 << endl;
    }
    */
    bt flatresult = flatbacktrack(node, x, y, trace, fcigar, node_list, trace_report);
    reverse(node_list.begin(), node_list.end());
    reverse(trace.begin(), trace.end());
    int softclipsAfter = read.length() - y;;
    int softclipsBefore = trace.front().y;
    t = trace.begin();
    Cigar fmcigar;
    fmcigar.append(Cigar(softclipsBefore, 'S'));
    for ( ; t != trace.end() && t != trace.end(); ++t) {
        reverse(t->cigar.begin(), t->cigar.end());
        fmcigar.append(t->cigar);
    }
    fmcigar.append(Cigar(softclipsAfter, 'S'));
	
    // cout << "gcigar: " << result.backstr << endl;
	
    x = flatresult.x;
    y = flatresult.y;
    //gcigar = result.backstr;
    // todo fcigar

    // py: return {'node_id':node_id, 'x':x, 'y':y, 'cigar':cigar, 'nodes':nodes}
    trace_report.x = x;
    trace_report.y = y;

    trace_report.fcigar = fmcigar;
    trace_report.gcigar = gcigarss.str();
	
    trace_report.node_list = node_list;
    trace_report.node_name = node_list.front()->name;
    trace_report.node = node_list.front();

    return flatresult;
}


/* 
 * The below does more than get a simple backtrace:  
 * It tries to recover all the aligned nodes - there is some rationale
 * to do this;  but it might not be very strong.  Might add / improve later
 * 
 * Will do final reversal in the end for C++ reasons
 * Might be a good idea to add a BAM-flattener (can also do in Python)
 */


bt graphbacktrack(sn* node, int x, int y, vector<bt>& trace, string& backstr, vector<sn*> &node_list) {
    /* (1) Recover score from the node if possible, check for index error
     * (2) If score is 0, return existing cigar, node, position
     * (3) Otherwise,  add the new string, and call the function again!
     * (4) Do it anyways, but if the node changes, append node name, and change x coordinate.
     *
     * TODO: Update to Reflect Gaps
     * Long Term Todo: See how many recursions break it! -- unlimited as long as we have tail recursion
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
        node_list.push_back(node);
        backtrace.node = node;
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
            backstr.append("I");
            y = y - 1;
        } else if (arrow == 's') {
            // py: backstr = 'I' + backstr
            backstr.append("D");
            x = x - 1;
        } else {
            cout<<"BackTrace Error: Unknown Type";     // add proper error checking
        }
        sn* new_node = node->matrix[y][x].parent;

        if (new_node->name != node->name) {   // might be better way to cmp
            x = new_node->seq_len;
            node_list.push_back(node);
            backtrace.node = node;
            backtrace.backstr = backstr;
            trace.push_back(backtrace);
            backstr.clear();
        }
		
        return graphbacktrack(new_node, x, y, trace, backstr, node_list);
    }
}


bt flatbacktrack(sn* node, int x, int y, vector<bt>& trace, Cigar& cigar,
                 vector<sn*> &node_list, mbt& trace_report) {
	
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
        cout<<"INDEX ERROR: "<<node->name<<x<<y<<score<<endl;
    }

    /* Read the Arrow Matrix & Iterate over the Backtrace */
    if (score == 0) {

        node_list.push_back(node);
        backtrace.node = node;
        if (node->isref) { // if we're in the reference coordinate space
            backtrace.cigar = cigar;
        } else {
            /*
            cout << "end using cigar " << cigar << " would insert sequence from "
                 << 0 << " of " << node->sequence.size() << "bp in "
                 << endl << node << endl << " to " << y << " of "
                 << cigar.readLen() << "bp in read " << trace_report.read << endl;
            */
            backtrace.x = 0; // matches to start of variant now
            trace_report.read.replace(y, cigar.readLen(), node->sequence);
            backtrace.cigar = node->cigar;
            reverse(backtrace.cigar.begin(), backtrace.cigar.end()); // re-reverse
            trace_report.qualities.replace(y, cigar.readLen(),
                                           string(node->sequence.size(), shortInt2QualityChar(30)));
            //trace_report.read.replace(y, backtrace.cigar.readLen(), node->sequence.substr(0, x));
            // and 1bp of reference before the implied divergent sequence
        }
        trace.push_back(backtrace);
        cigar.clear();
        return backtrace;

    } else {

        arrow = node->matrix[y][x].arrow;	  // TODO: check about pointers  ??
        if (arrow == 'm') {
            cigar.append(Cigar(1,'M'));
            x = x - 1;
            y = y - 1;
        } else if (arrow == 'x') {
            cigar.append(Cigar(1,'M'));
            x = x - 1;
            y = y - 1;
        } else if (arrow == 'u') {
            cigar.append(Cigar(1,'I'));
            y = y - 1;
        } else if (arrow == 's') {
            cigar.append(Cigar(1,'D'));
            x = x - 1;
        } else {
            cout<<"BackTrace Error: Unknown Type";     // add proper error checking
        }

        sn* new_node = node->matrix[y][x].parent;

        if (new_node->name != node->name) {   // might be better way to cmp
            if (node->isref) { // if we're in the reference coordinate space
                backtrace.cigar = cigar;
            } else {
                // YO  TODO... FLATTEN INSERTIONS into read space
                // TOOD you need the start of the traceback...  the coordinates of the score max in x and y
                if (node_list.empty()) {
                    /*
                    cout << "begin using cigar " << cigar << " would insert sequence from "
                         << 0 << " of " << node->sequence.size() << "bp in "
                         << endl << node << endl << " to " << y << " of "
                         << cigar.readLen() << "bp in read " << trace_report.read << endl;
                    */
                    trace_report.read.replace(y, cigar.readLen(), node->sequence);
                    trace_report.qualities.replace(y, cigar.readLen(),
                                                   string(node->sequence.size(), shortInt2QualityChar(30)));
                    backtrace.cigar = node->cigar;
                } else {
                    /*
                    cout << "using cigar " << cigar << " would insert sequence from "
                         << x << " of " << cigar.refLen() << "bp in "
                         << endl << node << endl << " to " << y << " of "
                         << cigar.readLen() << "bp in read " << trace_report.read << endl;
                    */
                    trace_report.read.replace(y, cigar.readLen(), node->sequence.substr(x, cigar.refLen()));
                    trace_report.qualities.replace(y, cigar.readLen(),
                                                   string(cigar.refLen(), shortInt2QualityChar(30)));
                    backtrace.cigar = node->cigar;
                }
                //trace_report.read.replace(y, backtrace.cigar.readLen(), node->sequence.substr(x));
                // and 1bp of reference before the implied divergent sequence
                //trace_report.read.append(node->sequence.substr(0, x));
                reverse(backtrace.cigar.begin(), backtrace.cigar.end()); // re-reverse
            }
            x = new_node->seq_len;
            node_list.push_back(node);
            backtrace.node = node;
            trace.push_back(backtrace);
            cigar.clear();
        }

        return flatbacktrack(new_node, x, y, trace, cigar, node_list, trace_report);
    }
}

