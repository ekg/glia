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

#include "gliamodels.h"
#include "nodealign.h"
#include "traceback.h"


// struct cmp_parent_depths;

int getDepth(sn* node);

sn* sequenceDagAlign(std::string sequence, std::vector<sn*> nlist, int maxdepth);

sn* gsw(std::string read, std::vector<sn*> nlist);

#endif