//
//  parameters.cpp
//  
//
//  Created by Deniz Kural on 2/19/13.
//
//

#include "parameters.h"

using namespace std;


void Parameters::simpleUsage(char ** argv) {
    cout
    << "usage: " << argv[0] << " -f [REFERENCE] [OPTIONS] >[OUTPUT]" << endl
    << "Please see README at http://www.github.com/denizkural/clia ." << endl
    << "author:   Deniz Kural <denizkural@gmail.com>" << endl;
}

void Parameters::usage(char** argv) {
    cout
    << "usage: " << argv[0] << " -f [REFERENCE] [OPTIONS] >[OUTPUT]" << endl
    << "Please see README at http://www.github.com/denizkural/clia ." << endl
    << "author:   Deniz Kural <denizkural@gmail.com>" << endl;
}


Parameters::Parameters(int argc, char* const argv) {
    
    if (argc == 1) {
        simpleUsage(argv);
        exit(1);
    }
    
    // record command line parameters
    commandline = argv[0];
    for (int i = 1; i < argc; ++i) {
        commandline += " ";
        commandline += argv[i];
    }
    
    // set defaults
    
    // i/o parameters:
    read_input = "";            // -s --sequence
    fastq_file = ""             // -q --fastq-file
    fasta = "";                 // -f --fasta-reference
    target = "";                // -t --target
    
    outputFile = "";            // -o --bam-output
    
    // operation parameters
    useFile = FALSE;            // -x --use-file
    alignReverse = TRUE;        // -r --reverse-complement
    
    int c; // counter for getopt
    
    static struct option long_options[] =
    {
        {"help", no_argument, 0, 'h'},
        {"sequence", required_argument, 0, 's'}
        {"fastq-file", required_argument, 0, 'q'}
        {"fasta-reference", required_argument, 0, 'f'}
        {"target", required_argument, 0, 't'}
        {"output-file", required_argument, 0, 'o'}
        {"use-file", no_argument, 0, 'x'}
        {"reverse-complement", no_argument, 0, 'r'}
        
        {0, 0, 0, 0}
        
    };
    
    while (true) {
        int option_index = 0;
        c = getopt_long(argc, argv,
                        "hcO4ZKjH0diN5aI_k=wluVXJY:b:G:M:x:@:A:f:t:r:s:v:n:B:p:m:q:R:Q:U:$:e:T:P:D:^:S:W:F:C:L:8:z:1:3:E:7:2:9:%:",
                        long_options, &option_index);
        
        if (c == -1) // end of options
            break;

        switch (c) {
            // -s --sequence
            case 's':
                read_input = optarg;
                break;
            // -q --fastq-file
            case 'q':
                fastq_file = optarg;
                break;
            // -f --fasta-reference
            case 'f':
                fasta_reference = optarg;
                break;
            // -t --target
            case 't':
                target = optarg;
                break;
            // -o --output-file
            case 'o':
                outputFile = optarg;
            // -x --use-file
            case 'x':
                useFile = TRUE;
            // -r --reverse-complement
                alignReverse = FALSE;

            // - --help
            case 'h':
                usage(argv);
                exit(0);
                break;
                
            case '?': // print a suggestion about the most-likely long option which the argument matches
            {
                string bad_arg(argv[optind - 1]);
                option* opt = &long_options[0];
                option* closest_opt = opt;
                int shortest_distance = levenshteinDistance(opt->name, bad_arg);
                ++opt;
                while (opt->name != 0) {
                    int distance = levenshteinDistance(opt->name, bad_arg);
                    if (distance < shortest_distance) {
                        shortest_distance = distance;
                        closest_opt = opt;
                    }
                    ++opt;
                }
                cerr << "did you mean --" << closest_opt->name << " ?" << endl;
                exit(1);
            }
                break;

            default:
                abort ();
        }
    }
    
    
    // TODO:  Add a lot more of this stuff.
    if (fasta_reference == "") {
        cerr << "Please specify a fasta reference file." << endl;
        exit(1);
    }
    
    if (useFile == TRUE && fastq_file == "") {
        cerr << "Please specify a fastq input file." << endl;
    }
    
}

