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
#include "vcflib/Variant.h"
#include "gliamodels.h"


using namespace std;


int constructDAG(std::vector<sn*> &nlist, std::string &targetSequence, 
		 std::vector<vcf::Variant> &variants, long offset);

#endif /* defined(CONSTRUCT_H) */
