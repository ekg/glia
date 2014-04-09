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
                 ReferenceMappings& ref_map,
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

    long prev_pos = offset;
    string p5_ref_seq;

    vector<gssw_node*> p5_var_nodes; // previous p5 nodes

    //vector<gssw_node*> ref_backbone;

    long int id = 0;

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

        int current_pos = (long int) var.position - 1;
        //cerr << var << endl;

        // Construct Right-Node
        p5_ref_seq = targetSequence.substr(prev_pos - offset, current_pos - prev_pos);

        Cigar* p5_ref_cigar = NULL;
        gssw_node* p5_ref_node = NULL;
        if (!p5_ref_seq.empty()) {
            //cigars.push_back(Cigar(convert(p5_ref_seq.size()) + "M"));
            //refpositions.push_back(prev_pos + offset);
            //p5_ref_cigar = &cigars.back();
            p5_ref_node = (gssw_node*)gssw_node_create((void*)NULL, graph->size, p5_ref_seq.c_str(), nt_table, score_matrix);
            gssw_graph_add_node(graph, p5_ref_node);
            ref_map.add_node(p5_ref_node, prev_pos, Cigar(convert(p5_ref_seq.size()) + "M"));
            //cerr << "connect to old p5 nodes" << endl;
            for (vector<gssw_node*>::iterator n = p5_var_nodes.begin(); n != p5_var_nodes.end(); ++n) {
                gssw_nodes_add_edge(*n, p5_ref_node);
            }
            p5_var_nodes.clear();
        }

        // new previous position is at the start of the variant
        prev_pos = current_pos + var.ref.size();

        // construct the ref node of variant
        //cigars.push_back(Cigar(convert(var.ref.size()) + "M"));
        //refpositions.push_back(current_pos + offset);
        //cerr << "cpo ref pos .back() " << refpositions.back() << endl;
        //Cigar* ref_cigar = &cigars.back();
        gssw_node* ref_node = (gssw_node*)gssw_node_create((void*)NULL, graph->size, var.ref.c_str(), nt_table, score_matrix);
        gssw_graph_add_node(graph, ref_node);
        ref_map.add_node(ref_node, current_pos, Cigar(convert(var.ref.size()) + "M"));
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

        for (vector<string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
            // create node for this alt, first establishing its ref-relative cigar
            //cigars.push_back(Cigar(vavs[*a]));
            //refpositions.push_back(current_pos + offset);
            //Cigar* cigar = &cigars.back();
            gssw_node* alt_node = (gssw_node*)gssw_node_create((void*)NULL, graph->size, a->c_str(), nt_table, score_matrix);
            // save in graph
            gssw_graph_add_node(graph, alt_node);
            ref_map.add_node(alt_node, current_pos, Cigar(vavs[*a]));
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

    int current_pos = targetSequence.size() + offset;//(long int) var.position - (long int) offset;
    p5_ref_seq = targetSequence.substr(prev_pos - offset, current_pos - prev_pos);
    //cerr << "getting " << prev_pos << " - " << offset << " for " << current_pos << " - " << prev_pos << " bp" << endl;

    Cigar* p5_ref_cigar = NULL;
    gssw_node* p5_ref_node = NULL;
    if (!p5_ref_seq.empty()) {
        //cigars.push_back(Cigar(convert(p5_ref_seq.size()) + "M"));
        //refpositions.push_back(prev_pos + offset);
        //p5_ref_cigar = &cigars.back();
        p5_ref_node = (gssw_node*)gssw_node_create((void*)NULL, graph->size, p5_ref_seq.c_str(), nt_table, score_matrix);
        gssw_graph_add_node(graph, p5_ref_node);
        ref_map.add_node(p5_ref_node, prev_pos, Cigar(convert(p5_ref_seq.size()) + "M"));
        //ref_backbone.push_back(p5_ref_node);
        //cerr << "connect to old p5 nodes" << endl;
        for (vector<gssw_node*>::iterator n = p5_var_nodes.begin(); n != p5_var_nodes.end(); ++n) {
            gssw_nodes_add_edge(*n, p5_ref_node);
        }
        p5_var_nodes.clear();
    }
    

}

pair<gssw_node*, gssw_node*>
divide_ref_path(map<long, gssw_node*>& ref_path,
                ReferenceMappings& ref_map,
                map<long, set<gssw_node*> >& nodes,
                long pos,
                int8_t* nt_table,
                int8_t* score_matrix,
                int& id) {

    cerr << "dividing path at " << pos << endl;

    // now do the same thing for our ending node
    map<long, gssw_node*>::iterator ref = ref_path.upper_bound(pos);
    if (ref == ref_path.end()) {
        cerr << "at last node" << endl;
    }
    for (map<long, gssw_node*>::iterator x = ref_path.begin(); x != ref_path.end(); ++x) {
        cerr << x->first << ":" << x->second << endl;
    }
    --ref;
    cerr << "ref len " << ref->second->len << endl;
    cerr << ref->second << endl;
    //--ref; // step to previous
    // we should now be pointing to node past where we should insert
    if (ref == ref_path.begin()) {
        cerr << "splitting first ref node" << endl;
        //cerr << "variant is out of bounds!!??" << endl;
        //exit(1);
    }
    cerr << "Getting positions and stuff" << endl;
    long ref_node_pos = ref->first;
    gssw_node* old_node = ref->second;
    cerr << "which by the way are " << ref_node_pos << " " << old_node << endl;
    
    if (ref_node_pos == pos) {
        cerr << "ref node position (" << ref_node_pos << ") == pos (" << pos << ")" << endl;
        map<long, gssw_node*>::iterator n = ref; --n;
        cerr << ref->second->id << ":" << ref->second->seq << endl
             << n->second->id << ":" << n->second->seq << endl;
        return make_pair(n->second, ref->second);
    } else {
        // divide the ref node at our alt starting position
        int diff = pos - ref_node_pos;
        cerr << "pos = " << pos << " " << " ref node pos = " << ref_node_pos << endl;
        cerr << "ref len " << ref->second->len << endl;
        cerr << "ref = " << ref->second->seq << endl;
        // make our right node
        string left_ref_node_seq = string(ref->second->seq);//, diff);
        left_ref_node_seq = left_ref_node_seq.substr(0, diff);
        cerr << "left_ref_node_seq = " << left_ref_node_seq << endl;
        gssw_node* left_ref_node = (gssw_node*)gssw_node_create((void*)NULL, id++, left_ref_node_seq.c_str(), nt_table, score_matrix);
        // replace node connections
        gssw_node** p = old_node->prev;
        for (int i = 0; i < old_node->count_prev; ++i, ++p) {
            ref_map.del_edge(*p, old_node);
            ref_map.add_edge(*p, left_ref_node, ref_node_pos, Cigar(0, 'M'));
            gssw_node_replace_next(*p, old_node, left_ref_node);
            gssw_node_add_prev(left_ref_node, *p);
        }

        // make our right node
        //diff = diff;
        cerr << ref->second->len << endl;
        cerr << ref->second->seq << endl;
        string right_ref_node_seq = string(ref->second->seq);
        cerr << "reee... " << right_ref_node_seq << endl;
        cerr << diff << endl;
        right_ref_node_seq = right_ref_node_seq.substr(diff);
        cerr << "right_ref_node_seq = " << right_ref_node_seq << endl;
        gssw_node* right_ref_node = (gssw_node*)gssw_node_create((void*)NULL, id++, right_ref_node_seq.c_str(), nt_table, score_matrix);
        // ahem and connect current to previous
        gssw_node** n = old_node->next;
        for (int i = 0; i < old_node->count_next; ++i, ++n) {
            ref_map.del_edge(old_node, *n);
            ref_map.add_edge(right_ref_node, *n, pos + diff, Cigar(0, 'M'));
            gssw_node_replace_prev(*n, old_node, right_ref_node);
            gssw_node_add_next(right_ref_node, *n);
        }

        cerr << "hemmememmem" << endl;
        // connect left to right
        gssw_nodes_add_edge(left_ref_node, right_ref_node);
        ref_map.add_edge(left_ref_node, right_ref_node, pos, Cigar(0, 'M'));

        cerr << "destruction" << endl;
        // destroy old node
        gssw_node_destroy(old_node);
        nodes[ref_node_pos].erase(old_node);
        ref_map.del_node(old_node);

        // replace with new ones
        // left
        cerr << "and renewal" << endl;
        cerr << ref_node_pos << " " << left_ref_node << endl;
        cerr << "nodes @ " << ref_node_pos << " = " << nodes[ref_node_pos].size() << endl;
        //nodes.erase(ref_node_pos);
        if (nodes[ref_node_pos].size() == 0) {
            nodes.erase(ref_node_pos);
        }
        nodes[ref_node_pos].insert(left_ref_node);
        cerr << "wheetet" << endl;
        cerr << 1 << endl;
        ref_path[ref_node_pos] = left_ref_node;
        cerr << 2 << endl;
        ref_map.add_node(left_ref_node, ref_node_pos, Cigar(left_ref_node_seq.size(), 'M'));
        cerr << 3 << endl;
        // right
        cerr << "also renewal for the left.  there ya go buddy." << endl;
        nodes[pos].insert(right_ref_node);
        cerr << 4 << endl;
        ref_path[pos] = right_ref_node;
        cerr << 5 << endl;
        ref_map.add_node(right_ref_node, pos, Cigar(right_ref_node_seq.size(), 'M'));
        cerr << " and hey!" << endl;

        return make_pair(left_ref_node, right_ref_node);
    }

}


int constructDAGProgressive(gssw_graph* graph,
                            ReferenceMappings& ref_map,
                            string &targetSequence,
                            string& sequenceName,
                            vector<vcf::Variant> &variants,
                            long offset,
                            int8_t* nt_table,
                            int8_t* score_matrix) {


// algorithm
// maintain a core reference path upon which we add new variants as they come
// addition procedure is the following
// find reference node overlapping our start position
// if it is already the end of a node, add the new node
// if it is not the end of a node, break it, insert edges from old->new
// go to end position of alt allele (could be the same position)
// if it already has a break, just point to the next node in line
// if it is not broken, break it and point to the next node
// add new node for alt alleles, connect to start and end node in reference path
// store the ref mapping as a property of the edges and nodes (this allows deletion edges and insertion subpaths)

// probably the name "Backbone" should be changed....?  confusing
// maybe to "reference mapping" or something?

    cerr << "here......!!" << endl;
    int id = 0;
    map<long, gssw_node*> reference_path;
    map<long, set<gssw_node*> > nodes; // for maintaining a reference-sorted graph
    cerr << "target sequence " << targetSequence << endl;
    cerr << "starts at " << offset << endl;

    gssw_node* ref_node = (gssw_node*)gssw_node_create((void*)NULL, id++, targetSequence.c_str(), nt_table, score_matrix);
    reference_path[offset] = ref_node;
    nodes[offset].insert(ref_node);
    ref_map.add_node(ref_node, offset, Cigar(convert(targetSequence.size()) + "M"));
    cerr << "added ref node " << ref_node->seq << endl;

    for(vector<vcf::Variant>::iterator it = variants.begin(); it != variants.end(); ++it) {

        vcf::Variant& var = *it;
        cerr << var << endl;
        int current_pos = (long int) var.position - 1;
        // decompose the alt
        map<string, vector<vcf::VariantAllele> > alternates = var.parsedAlternates();
        for (map<string, vector<vcf::VariantAllele> >::iterator va = alternates.begin(); va !=alternates.end(); ++va) {
            vector<vcf::VariantAllele>& alleles = va->second;

            for (vector<vcf::VariantAllele>::iterator a = alleles.begin(); a != alleles.end(); ++a) {

                vcf::VariantAllele& allele = *a;
                // reference alleles are provided naturally by the reference itself
                if (allele.ref == allele.alt) {
                    continue;
                }

                cerr << "variant allele " << allele << endl;
                long allele_start_pos = allele.position - 1;  // 0/1 based conversion... thanks vcflib!
                long allele_end_pos = allele_start_pos + allele.ref.size();
                cerr << "allele start, end = " << allele_start_pos << " " << allele_end_pos << endl;

                gssw_node* left_ref_node = NULL;
                gssw_node* middle_ref_node = NULL;
                gssw_node* right_ref_node = NULL;
                cerr << "here ... " << endl;
                pair<gssw_node*, gssw_node*> ref_nodes = divide_ref_path(reference_path,
                                                                         ref_map,
                                                                         nodes,
                                                                         allele_start_pos,
                                                                         nt_table,
                                                                         score_matrix,
                                                                         id);
                cerr << "returned from divide ref path " << ref_nodes.first << " and " << ref_nodes.second << endl;
                left_ref_node = ref_nodes.first;
                // if the ref portion of the allele is not empty, then we need to make another cut
                if (!allele.ref.empty()) {
                    ref_nodes = divide_ref_path(reference_path,
                                                ref_map,
                                                nodes,
                                                allele_end_pos,
                                                nt_table,
                                                score_matrix,
                                                id);
                    middle_ref_node = ref_nodes.first;
                    right_ref_node = ref_nodes.second;
                } else {
                    right_ref_node = ref_nodes.second;
                }

                // create a new alt node and connect the pieces from before
                if (!allele.alt.empty()) {
                    gssw_node* alt_node = (gssw_node*)gssw_node_create((void*)NULL, id++, allele.alt.c_str(), nt_table, score_matrix);
                    nodes[allele_start_pos].insert(alt_node);
                    //ref_map.add_node(alt_node, allele_start_pos, );
                    ref_map.add_node(alt_node, allele_start_pos, Cigar(allele));
                    gssw_nodes_add_edge(left_ref_node, alt_node);
                    ref_map.add_edge(left_ref_node, alt_node, allele_start_pos, Cigar(0, 'M'));
                    gssw_nodes_add_edge(alt_node, right_ref_node);
                    // NB, should check that the allele_start_pos + allele.ref.size() == allele ref end, aka start of next ref
                    ref_map.add_edge(alt_node, right_ref_node, allele_start_pos + allele.ref.size(), Cigar(0, 'M'));
                    // the overlapping reference sequence is already connected
                } else {// otherwise, we have a deletion
                    gssw_nodes_add_edge(left_ref_node, right_ref_node);
                    ref_map.add_edge(left_ref_node, right_ref_node, allele_start_pos, Cigar(allele.ref.size(), 'D'));
                }
            }
        }
    }

    /*
    // add "N" reference node to start if we began with a variant
    // ensures proper fill in graph alignment
    if (nodes.begin()->second.size() > 1) {
        cerr << "adding null ref node" << endl;
        gssw_node* null_ref_node = (gssw_node*)gssw_node_create((void*)NULL, id++, "N", nt_table, score_matrix);
        long p = nodes.begin()->first;
        reference_path[p] = null_ref_node;
        set<gssw_node*>& nexts = nodes.begin()->second;
        nodes[p].insert(null_ref_node);
        ref_map.add_node(null_ref_node, p, Cigar("1M"));
        ref_map.add_edge(null_ref_node, reference_path[p+1], p+1, Cigar(0, 'M'));
        // add connections
        for (set<gssw_node*>::iterator n = nexts.begin(); n != nexts.end(); ++n) {
            gssw_nodes_add_edge(null_ref_node, *n);
        }
    }
    */

    // topologically sort nodes before inserting into gssw_graph
    // this should be done for us due to ordering in nodes[]

    list<gssw_node*> sorted_nodes;
    topological_sort(nodes, sorted_nodes);

    id = 0;
    for (list<gssw_node*>::iterator n = sorted_nodes.begin(); n != sorted_nodes.end(); ++n) {
        (*n)->id = id++;
        cout << "node: " << *n << " id " << (*n)->id << endl;
        gssw_graph_add_node(graph, *n);
    }

    /*
    //set<gssw_node*> added_nodes;
    // execute a breadth-first search to generate ordered nodes
    int depth = 0;
    gssw_node* n = *(nodes.begin()->second.begin());
    while (true) {
        for (int i = 0; i < n->count_next; ++i) {
            
        }
    }
    */
    

}


    /*
Tarjan's topological sort

L <- Empty list that will contain the sorted nodes
while there are unmarked nodes do
    select an unmarked node n
    visit(n) 
function visit(node n)
    if n has a temporary mark then stop (not a DAG)
    if n is not marked (i.e. has not been visited yet) then
        mark n temporarily
        for each node m with an edge from n to m do
            visit(m)
        mark n permanently
        add n to head of L
    */


void topological_sort(map<long, set<gssw_node*> >& nodes,
                      list<gssw_node*>& sorted_nodes) {
    set<gssw_node*> unmarked_nodes;
    set<gssw_node*> temporary_marks;
    for (map<long, set<gssw_node*> >::iterator s = nodes.begin(); s != nodes.end(); ++s) {
        set<gssw_node*>& ns = s->second;
        for (set<gssw_node*>::iterator n = ns.begin(); n != ns.end(); ++n) {
            unmarked_nodes.insert(*n);
        }
    }
    while (!unmarked_nodes.empty()) {
        gssw_node* node = *(unmarked_nodes.begin());
        visit_node(node,
                   sorted_nodes,
                   unmarked_nodes,
                   temporary_marks);
    }
}

void visit_node(gssw_node* node,
                list<gssw_node*>& sorted_nodes,
                set<gssw_node*>& unmarked_nodes,
                set<gssw_node*>& temporary_marks) {
    cerr << node << endl;
    /*
    if (temporary_marks.find(node) != temporary_marks.end()) {
        cerr << "cannot sort graph because it is not a DAG!" << endl;
        exit(1);
    }
    */
    if (unmarked_nodes.find(node) != unmarked_nodes.end()) {
        temporary_marks.insert(node);
        for (int i = 0; i < node->count_next; ++i) {
            visit_node(node->next[i],
                       sorted_nodes,
                       unmarked_nodes,
                       temporary_marks);
        }
        unmarked_nodes.erase(node);
        sorted_nodes.push_front(node);
    }
}
