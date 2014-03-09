#include "utility.h"

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

char shortInt2QualityChar(short i) {
    return static_cast<char>(i + 33);
}

bool allATGC(string& s) {
    for (string::iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
            return false;
        }
    }
    return true;
}

bool allN(string& s) {
    for (string::iterator c = s.begin(); c != s.end(); ++c) {
        if (*c != 'N') {
            return false;
        }
    }
    return true;
}

string graph_mapping_to_string(gssw_graph_mapping* gm) {
    gssw_graph_cigar* g = &gm->cigar;
    int32_t i;
    int32_t c = 0;
    gssw_node_cigar* nc = g->elements;
    stringstream gmss;
    gmss << gm->score << "@" << gm->position << ":";
    for (i = 0; i < g->length; ++i, ++nc) {
        gmss << nc->node->id << "[";
        gssw_cigar* c = nc->cigar;
        int j;
        int l = c->length;
        gssw_cigar_element* e = c->elements;
        for (j=0; j < l; ++j, ++e) {
            gmss << e->length << e->type;
        }
        gmss << "]";
    }
    return gmss.str();
}
