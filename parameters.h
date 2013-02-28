//
//  parameters.h
//  Based on marthlab/ekg param use of getopt
//  Created by Deniz Kural on 2/19/13.
//


#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <getopt.h>

using namespace std;

// Command line parameters Class
class Parameters {
    
    friend ostream &operator<<(ostream &out, const Parameters &p);
    
public:
    
    // i/o parameters:
    string read_input;           // -s --sequence
    string fastq_file;          // -q --fastq-file
    string fasta;               // -f --fasta-reference
    string target;             // -t --target

    string outputFile;          // -o --bam-output

    // operation parameters
    bool useFile;               // -x --fastq-file
    bool alignReverse;          // -r --reverse-complement

    
    // functions
    Parameters(int argc, char * const argv);
    void usage(char **argv);
    void simpleUsage(char **argv);
    
    //reporting
    string commandline;
    
};


#endif /* defined(PARAMETERS_H) */
