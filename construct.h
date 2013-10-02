//
//  construct.h
//  
//
//  Created by Deniz Kural on 2/19/13.
//
//

#ifndef CONSTRUCT_H
#define CONSTRUCT_H

#include <string>
#include <vector>
#include "Variant.h"
#include "gliamodels.h"
#include "cigar.h"


using namespace std;
using namespace vcf;

int constructDAG(std::vector<sn*> &nlist, std::string& targetSequence, std::string& sequenceName,
                 std::vector<vcf::Variant> &variants, long offset);

#endif /* defined(CONSTRUCT_H) */
