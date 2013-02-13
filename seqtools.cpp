/*
 *  seqtools.cpp
 *  glia
 *
 *  Created by Deniz Kural on 8/2/12.
 *  Copyright 2012. All rights reserved.
 *
 */

#include "seqtools.h"

using namespace std;

/* Returns the Reverse Complement of a DNA Sequence, from the alphabet {A,T,C,G,N} */
string reverseComplement(string read) {
	
	// Declare the (empty) reverse complement read as a string
	string rc_read;	
	
	// Reverse Read
	rc_read.assign(read.rbegin(), read.rend());
	
	// Complement.  Note that not IUPAC compliant. Uses the alphabet {A,T,C,G,N}
	string::iterator t;
	for (t = rc_read.begin(); t != rc_read.end(); ++t) {
		switch (*t) {
			case 'A':
				*t = 'T';
				break;
			case 'T':
				*t = 'A';
				break;
			case 'C':
				*t = 'G';
				break;
			case 'G':
				*t = 'C';
				break;
			case 'N':
				*t = 'N';
				break;
			default:
				cout << "Unknown Nucleotide!";
				break;
		}
	}
	
	// Return the Read  (faster if done through pointers?)
	return rc_read;
}
