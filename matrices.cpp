/*
 *  matrices.cpp
 *  glia
 *
 *  Created by Deniz Kural on 12/21/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "matrices.h"
#include <vector>
#include <iostream>
#include <string>

using namespace std;

int matrix_trial()
{
	vector<vector<int> > score_matrix (3, vector<int>(2,0));
	
	score_matrix[0][3] = 1;
	cout<<score_matrix[0][0];
	
	return 62;
}

int string_trial() {
	string st1 = "apple";
	string st2 = "apple";

	if(st1 == st2) {
		cout << "wow they really are equal!" << endl;
	}
	
	return 0;
}
	