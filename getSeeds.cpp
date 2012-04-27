/*
 *  getSeeds.cpp
 *  glia
 *
 *  Created by Deniz Kural on 1/24/12.
 *  Copyright 2012 Deniz Kural. All rights reserved.
 *
 */

#include "getSeeds.h"
using namespace std;


// compare vector pointers by integers pointed
struct cmp_seeds {
	bool operator () (const pair< vector<int>::iterator, vector<int>::iterator> p1, const pair< vector<int>::iterator, vector<int>::iterator> p2) const {
		return *(p1.first) >  *(p2.first);      
	}
};


/* short description: look up the read and determine GSW candidates.
   long description: look up each kmker in read from the hashtable, which returns a sorted vector, and cluster online through feeding the vectors through a priority queue.
   inputs: read, the hash, hash size
   outputs: candidates for SW
 */

int lookupRead(string read, google::sparse_hash_map<string, vector<int> > &ghash, int hash_size) {
	
	
	/* clustering logic */

	int last_position;						// todo: check to see if 0-based or 1-based genome
	int current_position;
	int cluster_begin;
	
	int cluster_distance = 5;				//  super important
	
	vector<pair<int, int> > clusters;
	
	
	/* clustering logic end */

	
	cmp_seeds mycompare;
	
	pair< vector<int>::iterator, vector<int>::iterator> c_lead;
	priority_queue< pair< vector<int>::iterator, vector<int>::iterator>, vector<pair< vector<int>::iterator, vector<int>::iterator> >, cmp_seeds> cluster_feeder(mycompare);
	
		
	// look up kmers in hash
	string read_kmer;
	int i_limit = read.size() - hash_size + 1;
	

	for (int i = 0; i < i_limit; ++i) {
		read_kmer = read.substr(i, hash_size);
		if (ghash[read_kmer].empty() != true) {				// check if key exists in hash
			cluster_feeder.push( make_pair( ghash[read_kmer].begin(), ghash[read_kmer].end()));
		}
	}
	
	/* clustering logic */
	
	if (cluster_feeder.empty() != true) {
		c_lead = cluster_feeder.top();
		last_position = *(c_lead.first);
		cluster_begin = *(c_lead.first);
	}
	
	/* clustering logic end */

	
	// sort all positions
		
	while ( cluster_feeder.empty() != true) {					// this while belongs under the if  [huh?]
		c_lead = cluster_feeder.top();
		cluster_feeder.pop();
	
		/* clustering logic */
		
		current_position = *(c_lead.first); 
		
		if ( current_position - last_position < cluster_distance ) {
			// extend cluster
			last_position = current_position;
		} else {
			// finish old cluster, add to cluster library, 
			clusters.push_back(make_pair(cluster_begin, last_position));

			// start new cluster
			last_position = current_position;
			cluster_begin = current_position;
		}
		
		//cout << *(c_lead.first) << " " << *(c_lead.second) << endl;
        //this could be a "yield" to a separate clustering component

		/* end clustering logic */
		
		c_lead.first++;
		
		if (c_lead.first != c_lead.second) {
			cluster_feeder.push(c_lead);
		}
	}
	
	// todo: pretty print clusters
	// todo: return clusters
		
	return 0;
}



/*
 *
 * - iterate read & recover kmers
 *
 * for each contig:
 *  - look up kmers in hash
 *  - sort all positions
 *  - cluster
 *  - return clusters
 */

