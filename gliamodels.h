/*
 *  gliamodels.h
 *  glia
 *
 *  Created by Deniz Kural.
 *  Copyright 2010-2012 Deniz Kural. All rights reserved.
 *
 */


#ifndef GLIAMODELS_H  
#define GLIAMODELS_H

#include <iostream>
#include <vector>
#include <string>
#include <iostream>

// forward declaration of structures
struct ms;
struct ts;
struct sn;


// ts := top score data structure,  3-tuplet of int
// keeps track of matrix node with highest score
struct ts {
    int score;
    int x;
    int y;
ts() :  score(0),
	x(0),
	y(0) { }
};


// ms := matrix score structure
struct ms {
    int score;
    char arrow;
    sn* parent;
ms() :  score(0),
	arrow(0),
	parent(NULL) { }
};


// sn := string-node data structure
struct sn {
    std::string sequence;
    std::string name;
    std::vector<sn*> p5;			// now this will be sn pointers
    std::vector<sn*> p3;			// now this will be sn pointers
    int seq_len;
    ts top_score;
    std::vector<std::vector<ms> > matrix;
    int depth;			        // leave uninitialized? NO!!!!!!!!

// initialization
sn() :  seq_len(0),
	depth(-1) { }

    void initScore(size_t read_length);

};

std::ostream& operator<<(std::ostream& o, const sn* s);
std::ostream& operator<<(std::ostream& o, const ts& t);


#endif
