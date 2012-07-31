/*
 *  getSeeds.cpp
 *  glia
 *
 *  Created by Deniz Kural on 1/24/12.
 *  Copyright 2010-2012 Deniz Kural. All rights reserved.
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
   inputs: read, the hash (of a single string), hash size
   outputs: candidates for SW
 
   notes: Probably should give it an object to return
 */

int lookupRead(string read, google::sparse_hash_map<string, vector<int> > &ghash, int hash_size) {
	
	
	/* clustering logic */

	int last_position;						// todo: check to see if 0-based or 1-based genome
	int current_position;					// line of last extension! 
	int cluster_begin;
	
	/* Issue to test: Currently, if cluster_distance does NOT assume int + hash_size when calculating end of cluster! */
	int cluster_distance = 5;				//  super important. essentially defines max indel size.
	
	
	// Declare the container holding the clusters: defined as a collection of intervals on a linear genome.
	vector<pair<int, int> > clusters;
	
	
	/* clustering logic end */

	
	cmp_seeds mycompare;
	
	// A pair of iterators
	pair< vector<int>::iterator, vector<int>::iterator> c_lead;
	
	// Priority Queue called cluster_feeder
	priority_queue< pair< vector<int>::iterator, vector<int>::iterator>, vector<pair< vector<int>::iterator, vector<int>::iterator> >, cmp_seeds> cluster_feeder(mycompare);
	
		
	/* look up kmers in hash in the following manner:
	   - Iterate over each kmer in a read in a windowed fashion, 
	   - see if it exists in the hashmap of the chromosome,  
	   - if so make a pair of iterators marking the beginning and end of the vector of integers returned by the hashmap (i.e. the matches for the kmer!)
	   - add it to the priority queue, that gets compared by the first values that the begin pointer points to.
	*/
	
	string read_kmer;
	int i_limit = read.size() - hash_size + 1;
	

	for (int i = 0; i < i_limit; ++i) {
		read_kmer = read.substr(i, hash_size);
		if (ghash[read_kmer].empty() != true) {				// check if key exists in hash
			cluster_feeder.push( make_pair( ghash[read_kmer].begin(), ghash[read_kmer].end()));
		}
	}
	
	
	/* clustering logic */
	
	// We look up the top of the queue only to seed last position and cluster begin. We do not pop to do it again.
	if (cluster_feeder.empty() != true) {
		c_lead = cluster_feeder.top();
		last_position = *(c_lead.first);
		cluster_begin = *(c_lead.first);
	
	/* clustering logic end */

	
		// sort all positions
		
		while ( cluster_feeder.empty() != true) {
			
		
			// We now start going through the queue - we lookup the top and this time actually pop
			c_lead = cluster_feeder.top();
			cluster_feeder.pop();
			
			
			/* clustering logic */
			
			current_position = *(c_lead.first); 
			
			if ( current_position - last_position < cluster_distance ) {
				// extend cluster
				last_position = current_position;
			} 
			else {
				// finish old cluster, add to cluster library, 
				clusters.push_back(make_pair(cluster_begin, last_position));				// NOTE: Consider adding hashsize -  last_position+hashsize
				
				//this could be a "yield" to a separate component
				cout << cluster_begin << " " << last_position << endl;
				
				//there could be a filtering mechanism here to get rid of unpromising clusters et cetera.
				
				// start new cluster
				last_position = current_position;
				cluster_begin = current_position;
			}
		
			//cout << *(c_lead.first) << " " << *(c_lead.second) << endl;					// legacy. note that c_lead.second is somewhat meaningless.
			/* end clustering logic */
		
			c_lead.first++;
			
			if (c_lead.first != c_lead.second) {
				cluster_feeder.push(c_lead);
			}
		}
		
		clusters.push_back(make_pair(cluster_begin, last_position));
		// if you add filtering / alignment make sure it is also called one last time here.
		
	}

	
	// print clusters
	/*
	vector<pair<int, int> >::iterator t;
	
	for (t = clusters.begin(); t != clusters.end(); ++t) {
		cout << t->first << " ... " << t->second << endl;
	}
	*/
	 
	
	// todo: return clusters
		
	return 0;
}



/*
 *
 * - iterate read & recover kmers
 *
 * for each contig:
 *  - look up kmers in hash
 *  - sort all positions .. happens while clustering through priority queue
 *  - cluster
 *  - return clusters
 */

