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


int constructDAG(gssw_graph* graph,
                 vector<Cigar> &cigars,
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

    long prev_pos = 0;
    string p5_ref_seq;

    vector<gssw_node*> p5_var_nodes; // previous p5 nodes

    vector<gssw_node*> ref_backbone;

    // in construction, we assume nodes are already in DAG form
    // really, this may be a problem; imagine overlapping variants
    // so can the construction algorithm be adjusted to allow this
    // to do so we should maintain the reference backbone and build new nodes off of this

    for(vector<vcf::Variant>::iterator it = variants.begin(); it != variants.end(); ++it) {

        vcf::Variant& var = *it;
        vector<gssw_node*> var_nodes;

        // two possibilities assuming complete positional normalization and no overlapping variants (TODO implement check!)
        // 1) we have a bubble
        //    - make a reference node from here to the last position, link it to the last alleles
        //    - make a node for each of the ref and alt alleles in this variant
        //    - link the last ref node to the new nodes representing the alleles at this variant,
        //    - and remember the ref and alts of this variant for next time around, so they can be linked downstream
        // 2) we have successive variants
        //    - don't bother making a reference node--- it's embedded in the last variant record
        //    - make a node for each of the ref and alt alleles in this variant
        //    - link them the ref and alt alleles in the last variant
        //    - remember the ref and alts from this variant for the next time around

        int current_pos = (long int) var.position - (long int) offset;
        //cerr << var << endl;

        // Construct Right-Node
        p5_ref_seq = targetSequence.substr(prev_pos, current_pos - prev_pos);

        // new previous position is at the start of the variant
        prev_pos = var.position + var.ref.size() - offset;

        Cigar* p5_ref_cigar = NULL;
        gssw_node* p5_ref_node = NULL;
        if (!p5_ref_seq.empty()) {
            cigars.push_back(Cigar(convert(p5_ref_seq.size()) + "M"));
            p5_ref_cigar = &cigars.back();
            p5_ref_node = (gssw_node*)gssw_node_create((void*)p5_ref_cigar, graph->size, p5_ref_seq.c_str(), nt_table, score_matrix);
            gssw_graph_add_node(graph, p5_ref_node);
            ref_backbone.push_back(p5_ref_node);
            //cerr << "connect to old p5 nodes" << endl;
            for (vector<gssw_node*>::iterator n = p5_var_nodes.begin(); n != p5_var_nodes.end(); ++n) {
                gssw_nodes_add_edge(*n, p5_ref_node);
            }
            p5_var_nodes.clear();
        }

        // construct the ref node of variant
        cigars.push_back(Cigar(convert(var.ref.size()) + "M"));
        Cigar* ref_cigar = &cigars.back();
        gssw_node* ref_node = (gssw_node*)gssw_node_create((void*)ref_cigar, graph->size, var.ref.c_str(), nt_table, score_matrix);
        gssw_graph_add_node(graph, ref_node);
        ref_backbone.push_back(ref_node);
        // store the current ref in the pp3 nodes for connection on next iteration
        var_nodes.push_back(ref_node);

        // if we have variants in succession, we need to attach the p5 nodes
        // to the current ref and alt nodes

        if (!p5_ref_node) {
            // transfer connections across to this ref
            for (vector<gssw_node*>::iterator n = p5_var_nodes.begin(); n != p5_var_nodes.end(); ++n) {
                gssw_nodes_add_edge(*n, ref_node);
            }
        } else {
            // otherwise we naturally connect ref to ref
            gssw_nodes_add_edge(p5_ref_node, ref_node);
        }

        // Fill and connect Allele Nodes to p3_ref_node

        // get the cigars for the alt alleles
        map<string, vector<vcf::VariantAllele> > vavs = var.parsedAlternates();

        int i = 1;
        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a, ++i) {
            // create node for this alt, first establishing its ref-relative cigar
            cigars.push_back(Cigar(vavs[*a]));
            Cigar* cigar = &cigars.back();
            gssw_node* alt_node = (gssw_node*)gssw_node_create((void*)cigar, graph->size, a->c_str(), nt_table, score_matrix);
            // save in graph
            gssw_graph_add_node(graph, alt_node);
            // retain for connection to ref p3_ref_node of next variant
            var_nodes.push_back(alt_node);
            if (!p5_ref_node) {
                // carry across previous alternates and ref, if we are at successive variants
                for (vector<gssw_node*>::iterator n = p5_var_nodes.begin(); n != p5_var_nodes.end(); ++n) {
                    gssw_nodes_add_edge(*n, alt_node);
                }
            } else {
                gssw_nodes_add_edge(p5_ref_node, alt_node);
            }
        }

        p5_var_nodes.clear();
        p5_var_nodes.insert(p5_var_nodes.begin(), var_nodes.begin(), var_nodes.end());

    }

    int current_pos = targetSequence.size() - offset;//(long int) var.position - (long int) offset;
    p5_ref_seq = targetSequence.substr(prev_pos, current_pos - prev_pos);

    Cigar* p5_ref_cigar = NULL;
    gssw_node* p5_ref_node = NULL;
    if (!p5_ref_seq.empty()) {
        cigars.push_back(Cigar(convert(p5_ref_seq.size()) + "M"));
        p5_ref_cigar = &cigars.back();
        p5_ref_node = (gssw_node*)gssw_node_create((void*)p5_ref_cigar, graph->size, p5_ref_seq.c_str(), nt_table, score_matrix);
        gssw_graph_add_node(graph, p5_ref_node);
        ref_backbone.push_back(p5_ref_node);
        //cerr << "connect to old p5 nodes" << endl;
        for (vector<gssw_node*>::iterator n = p5_var_nodes.begin(); n != p5_var_nodes.end(); ++n) {
            gssw_nodes_add_edge(*n, p5_ref_node);
        }
        p5_var_nodes.clear();
    }
    
    /*
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
    */
    

}
