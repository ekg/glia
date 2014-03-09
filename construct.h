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

struct BackboneElement {
    long int ref_position;
    Cigar cigar;
BackboneElement() : ref_position(0), cigar(Cigar()) { }
BackboneElement(long int r,
                Cigar c) : ref_position(r), cigar(c) { }
};

struct Backbone : map<gssw_node*, BackboneElement> {
    void add(gssw_node* n, long int r, Cigar c) {
        (*this)[n] = BackboneElement(r, c);
    }
};
// std::vector<Cigar> &cigars,
// std::vector<long int> &refpositions,


int constructDAG(gssw_graph* graph,
                 Backbone& backbone,
                 std::string& targetSequence,
                 std::string& sequenceName,
                 std::vector<vcf::Variant> &variants,
                 long offset,
                 int8_t* nt_table,
                 int8_t* score_matrix);

#endif /* defined(CONSTRUCT_H) */
