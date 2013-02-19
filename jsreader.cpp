/*
 *  jsreader.cpp
 *  glia
 *
 *  Created by Deniz Kural on 8/1/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "jsreader.h"

using namespace std;

int json_example(vector<sn*> &nlist) {
	
	map<string, sn*> node_lookup;
	//node_lookup["n1"] = *sn;
	
	
	//cout << "Starting JSON Example" << endl;
	Json::Value root;   // will contains the root value after parsing.
	Json::Reader dagreader;
	Json::StyledWriter dagwriter;

	cout << "Loading JSON file" << endl;
	ifstream input_json_file("/Users/kural/indel.json");

	bool parsingSuccessful = dagreader.parse( input_json_file, root );
	

	if ( !parsingSuccessful ) {
		
		// report to the user the failure and their locations in the document.
		cout << "Failed to parse configuration\n"
			 << dagreader.getFormattedErrorMessages();
    
		return 0;
	}
	
	// example: retrieve string:string pair
	//string encoding = dagroot.get("encoding", "UTF-8" ).asString();
	
	const Json::Value dag = root["dag"];
	
			
	for ( int index = 0; index < dag.size(); ++index ) {
		
		Json::Value node = dag[index];
		
		sn* cnode;
		cnode = new sn;
		
		cnode->name = node.get("id","NA").asString();
		cnode->sequence = node.get("sequence","").asString();
		cnode->seq_len = cnode->sequence.length();
		cnode->depth = node.get("depth",-1).asInt();
		
		// add to the lookup table
		node_lookup[cnode->name] = cnode;
		nlist.push_back(cnode); 
		
	}
	
	for ( int index = 0; index < dag.size(); ++index ) {
		
		Json::Value node = dag[index];

		Json::Value p3 = node["p3"];
		for (int j = 0; j < p3.size(); ++j ) {
			//cout << p3.get(j,"NA").asString() << endl;
			nlist[index]->p3.push_back(node_lookup[p3.get(j,"NA").asString()]);
		}
	
		Json::Value p5 = node["p5"];
		for (int j = 0; j < p5.size(); ++j ) {
			//cout << p5.get(j,"NA").asString() << endl;
			nlist[index]->p5.push_back(node_lookup[p5.get(j,"NA").asString()]);
		}		

	}
	
	cout << dag;
	return 0;
}
	
		
	
	
		
    //loadPlugIn( plugins[index].asString() );
	
	//setIndentLength( root["indent"].get("length", 3).asInt() );
	//setIndentUseSpace( root["indent"].get("use_space", true).asBool() );
	
	// ...
	// At application shutdown to make the new configuration document:
	// Since Json::Value has implicit constructor for all value types, it is not
	// necessary to explicitly construct the Json::Value object:

	//root["encoding"] = getCurrentEncoding();
	//root["indent"]["length"] = getCurrentIndentLength();
	//root["indent"]["use_space"] = getCurrentIndentUseSpace();
	
	// Make a new JSON document for the configuration. Preserve original comments.
	//string outputConfig = dagwriter.write( root );

	// You can also use streams.  This will put the contents of any JSON
	// stream at a particular sub-value, if you'd like.
	//cin >> dagroot["subtree"];
	
	// And you can write to a stream, using the StyledWriter automatically.
	//cout << root;
	
    //return 0;

