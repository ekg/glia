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
    map<gssw_node*, ReferenceMapping>::iterator get_node(gssw_node* n) {
        map<gssw_node*, ReferenceMapping>::iterator p = nodes.find(n);
        if (p == nodes.end()) {
            cerr << "ERROR: could not find reference mapping for node " << n << endl;
            exit(1);
        } else {
            return p;
        }
    }
    map<pair<gssw_node*, gssw_node*>, ReferenceMapping>::iterator get_edge(gssw_node* n, gssw_node* m) {
        map<pair<gssw_node*, gssw_node*>, ReferenceMapping>::iterator p = edges.find(make_pair(n, m));
        if (p == edges.end()) {
            cerr << "ERROR: could not find reference mapping for edge from node " << n << " to " << m << endl;
            exit(1);
        } else {
            return p;
        }
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

int constructDAGProgressive(gssw_graph* graph,
                            ReferenceMappings& ref_map,
                            string &targetSequence,
                            string& sequenceName,
                            vector<vcf::Variant> &variants,
                            long offset,
                            int8_t* nt_table,
                            int8_t* score_matrix);


void topological_sort(map<long, set<gssw_node*> >& nodes,
                      list<gssw_node*>& sorted_nodes);

void visit_node(gssw_node* node,
                list<gssw_node*>& sorted_nodes,
                set<gssw_node*>& unmarked_nodes,
                set<gssw_node*>& temporary_marks);

#endif /* defined(CONSTRUCT_H) */
