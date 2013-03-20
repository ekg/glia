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
    std::string node_name;				// why use outside of context?
    sn* node; // start node
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
                 std::vector<sn*> &node_list);

// mbt := data structure for trace report
struct mbt;

// master backtrack
bt master_backtrack(sn* node, mbt &trace_report);


#endif
