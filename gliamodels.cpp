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
	  << s->top_score << ";"
	  << s->depth << ";";
	o << "p5:";
	for (std::vector<sn*>::const_iterator n = s->p5.begin(); n != s->p5.end(); ++n) {
	    o << (*n)->name << (n+1 == s->p5.end() ? "" : ",");
	}
	o << ";";
	o << "p3:";
	for (std::vector<sn*>::const_iterator n = s->p3.begin(); n != s->p3.end(); ++n) {
	    o << (*n)->name << (n+1 == s->p3.end() ? "" : ",");
	}
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
    matrix.clear();
    matrix.resize(read_length+1, std::vector<ms>(seq_len+1, o));
}

int displayDAG(const sn* s) {
    std::cout << s << std::endl;
    for (std::vector<sn*>::const_iterator n = s->p3.begin(); n != s->p3.end(); ++n) {
	displayDAG(*n);
    }
}
