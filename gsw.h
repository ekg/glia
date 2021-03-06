/*
 *  gsw.h
 *  glia
 *
 *  Created by Deniz Kural.
 *  Copyright 2011 Deniz Kural. All rights reserved.
 *
 */

#ifndef GSW_H
#define GSW_H

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>

#include "gliamodels.h"
#include "nodealign.h"
#include "traceback.h"


// struct cmp_parent_depths;

/* Return the topological order of a node in a DAG
Check if it has depth, if so don't repeat the backtrack
Else, check if it has a parent node, if so recurse
Otherwise declare the node to be a head node */
int getDepth(sn* node);


sn* sequenceDagAlign(std::string sequence, std::vector<sn*> nlist, int maxdepth,
		     const int match, const int mism, const int gap);

sn* gsw(std::string read, std::vector<sn*> nlist,
	const int match, const int mism, const int gap);

#endif
