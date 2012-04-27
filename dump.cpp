/*
 *  dump.cpp
 *  glia
 *
 *  Created by Deniz Kural on 6/20/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include <fstream>
#include <iostream>
#include "dump.h"

using namespace std;

int load_file () {
	string line;
	ifstream myfile ("test.fastq");
	
	if (myfile.is_open()) {
		while (myfile.good()) {
			getline (myfile,line);
			cout << line << endl;
		}
		myfile.close();
		
	}
	else cout << "Unable to open read file"; 
	return 0;
}



int chplay() {
	unsigned char x, y;
	x = 3;
	y = 5;
	cout << x+y;
	
	return 0;
}



/* File writer prototype */
int write_file () {
	ofstream outfile;
	outfile.open("test.out");
	outfile << "Writing to the file.\n";
	outfile.close();
	return 0;
}