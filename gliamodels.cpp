/*
 *  gliamodels.cpp
 *  glia
 *
 *  Created by Deniz Kural on 12/26/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "gliamodels.h"

std::ostream& operator<<(std::ostream& o, const sn* s) {
    if (s) {
	o << (void*)s << ";"
	  << s->name << ":"
	  << s->sequence << ";"
	  << s->p5.front()->name << "-" << s->p3.front()->name << ";"
	  << s->top_score << ";"
	  << s->depth;
    }
    return o;
}

std::ostream& operator<<(std::ostream& o, const ts& t) {
    o << t.score << ";" << t.x << "," << t.y;
    return o;
}

std::ostream& operator<<(std::ostream& o, const ms& m) {
    o << m.score << ";" << m.arrow << "," << m.parent;
    return o;
}

void sn::initScore(size_t read_length) {
    ms o;
    o.parent = this;
    matrix.resize(read_length+1, std::vector<ms>(seq_len+1, o));
}
