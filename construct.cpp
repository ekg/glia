//
//  construct.cpp
//  
//
//  Created by Deniz Kural and Erik Garrison on 2/19/13.
//
//


#include "construct.h"
#include "convert.h"

using namespace std;
//using namespace vcf;


int constructDAG(vector<sn*> &nlist, string &targetSequence, string& sequenceName,
                 vector<vcf::Variant> &variants, long offset) {

    /*
    cerr << "constructing DAG over "
         << targetSequence.size()
         << " and " << variants.size()
         << " variants with offset " << offset << endl;

    cerr << "target:\t" << targetSequence << endl;
    */

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

        // TODO cleanup 0M reference hacks

        // stash p3_ and current ref nodes

        // if we have a 0M reference node, we need to attach the pp3 nodes
        // to the current ref and alt nodes and ignore this node.
        // this boolean tells us when.
        bool zero_length_p3_ref = p3_ref_node->cigar.refLen() == 0;

        if (!zero_length_p3_ref) {
            nlist.push_back(p3_ref_node);
        }
        nlist.push_back(ref_node);

        // connect to old p3 nodes
        for (vector<sn*>::iterator n = pp3_var_nodes.begin(); n != pp3_var_nodes.end(); ++n) {
            p3_ref_node->p3.push_back(*n);
            if (!zero_length_p3_ref) {
                (*n)->p5.push_back(p3_ref_node);
            }
        }
        pp3_var_nodes.clear();

        // connect p3_ref <-> ref
        if (zero_length_p3_ref) {
            // if 0-length, transfer connections across
            for (vector<sn*>::iterator p = p3_ref_node->p3.begin(); p != p3_ref_node->p3.end(); ++p) {
                (*p)->p5.push_back(ref_node);
                ref_node->p3.push_back(*p);
            }
        } else {
            p3_ref_node->p5.push_back(ref_node);
            ref_node->p3.push_back(p3_ref_node);
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
            // save in nlist
            nlist.push_back(alt_node);
            // retain for connection to ref p3_ref_node of next variant
            pp3_var_nodes.push_back(alt_node);
            if (!zero_length_p3_ref) {
                // connect to current p3_ref_node
                alt_node->p3.push_back(p3_ref_node);
                // and connect the p5 of the p3_ref_node to the alt node
                p3_ref_node->p5.push_back(alt_node);
            } else {
                for (vector<sn*>::iterator p = p3_ref_node->p3.begin();
                     p != p3_ref_node->p3.end(); ++p) {
                    (*p)->p5.push_back(alt_node);
                    alt_node->p3.push_back(*p);
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

    bool zero_length_p3_ref = p3_ref_node->cigar.refLen() == 0;
    if (!zero_length_p3_ref) { // save only if there is sequence data
        // connect to old p3 nodes
        for (vector<sn*>::iterator n = pp3_var_nodes.begin(); n != pp3_var_nodes.end(); ++n) {
            p3_ref_node->p3.push_back(*n);
            (*n)->p5.push_back(p3_ref_node);
        }
        pp3_var_nodes.clear();
        
        // save reference
        nlist.push_back(p3_ref_node);
    }

    reverse(nlist.begin(), nlist.end());

}
