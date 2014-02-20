//
//  construct.cpp
//  
//
//  Created by Deniz Kural and Erik Garrison on 2/19/13.
//
//


#include "construct.h"
#include "convert.h"
#include "gssw.h"

using namespace std;
//using namespace vcf;


int constructDAG(vector<sn*> &nlist,
                 gssw_graph* graph,
                 string &targetSequence,
                 string& sequenceName,
                 vector<vcf::Variant> &variants,
                 long offset,
                 int8_t* nt_table,
                 int8_t* score_matrix) {

    /*
    cerr << "constructing DAG over "
         << targetSequence.size()
         << " and " << variants.size()
         << " variants with offset " << offset << endl;

    cerr << "target:\t" << targetSequence << endl;
    */

    //graph = gssw_graph_create(1);

    long prev_pos = targetSequence.size();
    string p3_ref_seq;

    vector<sn*> pp3_var_nodes; // previous p3 nodes

    for(vector<vcf::Variant>::reverse_iterator rit = variants.rbegin(); 
        rit != variants.rend(); ++rit) {

        vcf::Variant& var = *rit;
        int current_pos = (long int) var.position - (long int) offset + var.ref.size();

        // Construct Right-Node
        p3_ref_seq = targetSequence.substr(current_pos, prev_pos - current_pos);

        // new previous position is at the start of the variant
        prev_pos = var.position - offset;

        // Construct Right Node
        sn* p3_ref_node = new sn(
            p3_ref_seq
            ,
            sequenceName
            + ":"
            + convert(var.position + var.ref.size())
            + "-"
            + convert(var.position
                      + var.ref.size()
                      + p3_ref_seq.size())
            + ".ref.r"
            ,
            var.position + var.ref.size() // 1-based?
            ,
            Cigar(convert(p3_ref_seq.size()) + "M") // cigar
            );
        gssw_node* gssw_p3_ref_node = (gssw_node*)gssw_node_create((void*)p3_ref_node, graph->size, p3_ref_seq.c_str(), nt_table, score_matrix);
        p3_ref_node->node = gssw_p3_ref_node;

        // construct the ref node of variant
        sn* ref_node = new sn(
            var.ref
            ,
            sequenceName
            + ":"
            + convert(var.position)
            + "-"
            + convert(var.position
                      + var.ref.size())
            + ".ref.0"
            ,
            var.position // 1-based?
            ,
            Cigar(convert(var.ref.size()) + "M") // cigar
            );
        gssw_node* gssw_ref_node = (gssw_node*)gssw_node_create((void*)ref_node, graph->size+1, var.ref.c_str(), nt_table, score_matrix);
        ref_node->node = gssw_ref_node;

        // TODO cleanup 0M reference hacks

        // stash p3_ and current ref nodes

        // if we have a 0M reference node, we need to attach the pp3 nodes
        // to the current ref and alt nodes and ignore this node.
        // this boolean tells us when.
        bool zero_length_p3_ref = p3_ref_node->cigar.refLen() == 0;

        if (!zero_length_p3_ref) {
            nlist.push_back(p3_ref_node);
            gssw_graph_add_node(graph, gssw_p3_ref_node);
        }
        nlist.push_back(ref_node);
        gssw_graph_add_node(graph, gssw_ref_node);

        // connect to old p3 nodes
        cerr << "connect to old p3 nodes" << endl;
        for (vector<sn*>::iterator n = pp3_var_nodes.begin(); n != pp3_var_nodes.end(); ++n) {
            p3_ref_node->p3.push_back(*n);
            if (!zero_length_p3_ref) {
                (*n)->p5.push_back(p3_ref_node);
                gssw_nodes_add_edge(gssw_p3_ref_node, (*n)->node);
            }
        }
        pp3_var_nodes.clear();

        // connect p3_ref <-> ref
        cerr << "connect p3_ref <-> ref" << endl;
        if (zero_length_p3_ref) {
            // if 0-length, transfer connections across
            for (vector<sn*>::iterator p = p3_ref_node->p3.begin(); p != p3_ref_node->p3.end(); ++p) {
                (*p)->p5.push_back(ref_node);
                ref_node->p3.push_back(*p);
                gssw_nodes_add_edge(gssw_ref_node, (*p)->node);
            }
        } else {
            p3_ref_node->p5.push_back(ref_node);
            ref_node->p3.push_back(p3_ref_node);
            gssw_nodes_add_edge(gssw_ref_node, gssw_p3_ref_node);
        }

        // store the current ref in the pp3 nodes for connection on next iteration
        pp3_var_nodes.push_back(ref_node);

        // Fill and connect Allele Nodes to p3_ref_node

        // get the cigars for the alt alleles
        map<string, vector<vcf::VariantAllele> > vavs = var.parsedAlternates();

        int i = 1;
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++i) {
            Cigar cigar(vavs[*a]);
            sn* alt_node = new sn(
                *a
                ,
                sequenceName
                + ":"
                + convert(var.position)
                + "-"
                + convert(var.position
                          + var.ref.size())
                + ".alt." + convert(i)
                ,
                var.position
                ,
                cigar
                );
            gssw_node* gssw_alt_node = (gssw_node*)gssw_node_create((void*)alt_node, graph->size+(i-1), a->c_str(), nt_table, score_matrix);
            alt_node->node = gssw_alt_node;
            // save in nlist
            nlist.push_back(alt_node);
            // save in graph
            gssw_graph_add_node(graph, gssw_alt_node);
            // retain for connection to ref p3_ref_node of next variant
            pp3_var_nodes.push_back(alt_node);
            if (!zero_length_p3_ref) {
                // connect to current p3_ref_node
                alt_node->p3.push_back(p3_ref_node);
                // and connect the p5 of the p3_ref_node to the alt node
                p3_ref_node->p5.push_back(alt_node);
                // connect to p3_ref_node
                cerr << "connect to p3_ref_node" << endl;
                gssw_nodes_add_edge(gssw_alt_node, gssw_p3_ref_node);
            } else {
                for (vector<sn*>::iterator p = p3_ref_node->p3.begin();
                     p != p3_ref_node->p3.end(); ++p) {
                    (*p)->p5.push_back(alt_node);
                    alt_node->p3.push_back(*p);
                    cerr << "connect alt to p3_ref p3" << endl;
                    gssw_nodes_add_edge(gssw_alt_node, (*p)->node);
                }
            }
        }

    }

    // last node construction and connection

    p3_ref_seq = targetSequence.substr(0, prev_pos);

    sn* p3_ref_node = new sn(
        p3_ref_seq
        ,
        sequenceName
        + ":"
        + convert(offset)
        + "-"
        + convert(prev_pos + offset)
        + ".ref.r"
        ,
        offset
        ,
        Cigar(convert(prev_pos) + "M")
        );
    gssw_node* gssw_p3_ref_node = (gssw_node*)gssw_node_create((void*)p3_ref_node, graph->size, p3_ref_seq.c_str(), nt_table, score_matrix);
    p3_ref_node->node = gssw_p3_ref_node;

    bool zero_length_p3_ref = p3_ref_node->cigar.refLen() == 0;
    if (!zero_length_p3_ref) { // save only if there is sequence data
        // connect to old p3 nodes
        for (vector<sn*>::iterator n = pp3_var_nodes.begin(); n != pp3_var_nodes.end(); ++n) {
            p3_ref_node->p3.push_back(*n);
            (*n)->p5.push_back(p3_ref_node);
            gssw_nodes_add_edge(p3_ref_node->node, (*n)->node);
        }
        pp3_var_nodes.clear();
        
        // save reference
        nlist.push_back(p3_ref_node);
        gssw_graph_add_node(graph, p3_ref_node->node);
    }

    reverse(nlist.begin(), nlist.end());
    // now... to reverse the graph
    // which is partially-ordered now btw (but reversed).  great.
    int j = graph->size-1;
    for (int i = 0; i < graph->size/2; ++i) {
        gssw_node* tmp = graph->nodes[i];
        graph->nodes[i] = graph->nodes[j - i];
        graph->nodes[j - i] = tmp;
    }
    

}
