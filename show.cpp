/*
 *  show.cpp
 *  glia
 *
 *  Created by Deniz Kural on 1/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "show.h"

using namespace std;


/*
vector<string> parentNames(sn* node) {
	
	vector<string> parentNames;
	
	vector<sn*>::iterator t;									// good place for vector iteration?
	for (t=node.p5.begin(); t!=node.p5.end(); ++t) {
		parent_names.push_back((*t)->name);
	}
	
	return parentNames;
}

int displayStrings(vector<string> values) {
	
	vector<sn*>::iterator t;									// good place for vector iteration?
	for (t=node.p5.begin(); t!=node.p5.end(); ++t) {
		parent_names.push_back((*t)->name);
	}
	
*/

int displayNode(sn* node) {

	cout << "Node Summary" << endl << "----------" << endl;
	cout << "NodeName: " << node->name << endl;
	cout << "Sequence: " << node->sequence << endl;
	cout << "Depth:    " << node->depth << endl;
	cout << "Parents:  ";
	
	vector<sn*>::iterator t;									
	
	for (t = node->p5.begin(); t != node->p5.end(); ++t) {
		cout<<(*t)->name<<", ";
	}
	cout<<endl;
	
	cout << "Children: ";
		
	for (t = node->p3.begin(); t != node->p3.end(); ++t) {
		cout<<(*t)->name<<", ";
	}
	cout<<endl<<endl;
	
	return 0;
}




int displayAlignment(sn* node) {
	
	
	vector<vector<int> >::iterator y;
	vector<int>::iterator x;
	
	cout << "Score:"<< endl << "----------" << endl;
	
	for (y = node->score_matrix.begin(); y != node->score_matrix.end(); ++y) {
		cout << "\t";
		for (x = (*y).begin(); x != (*y).end();  ++x) {
			
				 cout << *x << "\t";
		}
		cout << endl;
	}
	
	cout << endl;
	cout << "Back Trace:" << endl << "----------" << endl;
	
	vector<vector<char> >::iterator p;
	vector<char>::iterator q;
	
	for (p = node->arrow_matrix.begin(); p != node->arrow_matrix.end(); ++p) {
		cout << "\t";
		for (q = (*p).begin(); q != (*p).end();  ++q) {
			
			cout << *q << "\t";
		}
		cout << endl;
	}
	
	cout << endl;
	cout << "Parent:" << endl << "----------" << endl;
	
	vector<vector<sn*> >::iterator a;
	vector<sn*>::iterator b;
	
	for (a = node->parent_matrix.begin(); a != node->parent_matrix.end(); ++a) {
		cout << "\t";
		for (b = (*a).begin(); b != (*a).end();  ++b) {
			
			cout << (*b)->name << "\t";
		}
		cout << endl;
	}
	
	cout << endl;
	return 0;
}
	
	


