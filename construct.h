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
#include "gssw.h"


using namespace std;
//using namespace vcf;

int constructDAG(gssw_graph* graph,
                 std::vector<Cigar> &cigars,
                 std::string& targetSequence,
                 std::string& sequenceName,
                 std::vector<vcf::Variant> &variants,
                 long offset,
                 int8_t* nt_table,
                 int8_t* score_matrix);

#endif /* defined(CONSTRUCT_H) */
