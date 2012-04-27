/*  
 * Python Notice:
 * Created on Apr 15, 2011
 * @author: kural
 *
 * Functions here contain the graph alignment algorithm,
 * implemented in the block-matrix approach.
 *
 * nodealign.h
 * glia
 *
 * Created by Deniz Kural on 6/28/11.
 * Copyright 2011 Deniz Kural. All rights reserved.
 *
 */

#ifndef NODEALIGN_H  
#define NODEALIGN_H


#include <vector>
#include <string>
#include <iostream>
#include "gliamodels.h"
#include "show.h"


// compare parent nodes by score.
struct cmp_parent_nodes;

/* Main Alignment Algorithm for String and Node */
int StringNodeAlign(std::string read, int read_length, sn &node);


#endif
