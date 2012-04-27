/*
 *  getSeeds.h
 *  glia
 *
 *  Created by Deniz Kural on 1/24/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef GETSEEDS_H
#define GETSEEDS_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <queue>

#include <google/sparse_hash_map>

#include "fastq.h"

int lookupRead(std::string read, google::sparse_hash_map<std::string, std::vector<int> > &ghash, int hash_size);


#endif