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


// ts := top score data structure,  3-tuplet of int
// keeps track of matrix node with highest score
struct ts {
	int score;
	int x;
	int y;
};


// sn := string-node data structure
struct sn {
	std::string sequence;
	std::string name;
	std::vector<sn*> p5;			// now this will be sn pointers
	std::vector<sn*> p3;			// now this will be sn pointers
	int seq_len;
	ts top_score;
	std::vector<std::vector<int> > score_matrix;
	std::vector<std::vector<char> > arrow_matrix;
	std::vector<std::vector<sn*> > parent_matrix;
	int depth;						// leave uninitialized?
};


#endif
