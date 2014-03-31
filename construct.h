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

struct ReferenceMapping {
    long int ref_position;
    Cigar cigar;
ReferenceMapping() : ref_position(0), cigar(Cigar()) { }
ReferenceMapping(long int r,
                 Cigar c) : ref_position(r), cigar(c) { }
    bool is_ref(void) {
        return cigar.size() == 1 && cigar.front().type == 'M';
    }
};

struct ReferenceMappings {
    map<gssw_node*, ReferenceMapping> nodes;
    map<pair<gssw_node*, gssw_node*>, ReferenceMapping> edges;
    void add_node(gssw_node* n, long int r, Cigar c) {
        nodes[n] = ReferenceMapping(r, c);
    }
    void del_node(gssw_node* n) {
        nodes.erase(n);
    }
    void add_edge(gssw_node* n, gssw_node* m, long int r, Cigar c) {
        edges[make_pair(n, m)] = ReferenceMapping(r, c);
    }
    void del_edge(gssw_node* n, gssw_node* m) {
        edges.erase(make_pair(n, m));
    }
    ReferenceMapping& get_node(gssw_node* n) {
        return nodes[n];
    }
    ReferenceMapping& get_edge(gssw_node* n, gssw_node* m) {
        return edges[make_pair(n, m)];
    }
    bool empty(void) {
        return nodes.size() == 0 && edges.size() == 0;
    }
    void clear(void) {
        nodes.clear();
        edges.clear();
    }
};


int constructDAG(gssw_graph* graph,
                 ReferenceMappings& refMappings,
                 std::string& targetSequence,
                 std::string& sequenceName,
                 std::vector<vcf::Variant> &variants,
                 long offset,
                 int8_t* nt_table,
                 int8_t* score_matrix);

#endif /* defined(CONSTRUCT_H) */
