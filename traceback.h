/*
 *  traceback.h
 *  glia
 *
 *  Created by Deniz Kural on 12/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TRACEBACK_H  
#define TRACEBACK_H

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>

#include "gliamodels.h"
#include "nodealign.h"
#include "cigar.h"
#include "utility.h"


// bt := backtrace return object
struct bt {
    int x;
    int y;
    std::string backstr;
    Cigar cigar;
    sn* node;
bt() : node(NULL) { }
};


// mbt := data structure for trace report
struct mbt {
    int x;
    int y;
    std::string gcigar; // graph cigar
    Cigar fcigar; // flattened to reference
    std::vector<sn*> node_list;
    std::string node_name;
    sn* node; // start node
    // when flattening, we might want to pull in sequence from
    // nodes which are only partly overlapped, and then append it
    // to our alignments.  these hold such sequence when it exists.
    std::string read; // modified, possibly flattened read
    std::string qualities; // modified, possibly flattened quality bases
};


// recursive backtrack
bt graphbacktrack(sn* node,
	     int x, int y,
	     std::vector<bt>& trace,
	     std::string& backstr,
	     std::vector<sn*> &node_list);

bt flatbacktrack(sn* node,
                 int x, int y,
                 std::vector<bt>& trace,
                 Cigar& cigar,
                 std::vector<sn*> &node_list,
                 mbt& trace_report);

// mbt := data structure for trace report
struct mbt;

// master backtrack
bt master_backtrack(sn* node, mbt &trace_report, std::string& read, std::string& qualities);


#endif
