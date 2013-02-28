//
//  parameters.cpp
//  
//
//  Created by Deniz Kural on 2/19/13.
//
//

#include "parameters.h"

using namespace std;


// Levenshtein Distance Algorithm: C++ Implementation
// by Anders Sewerin Johansen
// http://www.merriampark.com/ldcpp.htm

#include <string>
#include <vector>

int levenshteinDistance(const std::string source, const std::string target) {
    // Step 1
    const int n = source.length();
    const int m = target.length();
    if (n == 0) {
        return m;
    }
    if (m == 0) {
        return n;
    }
    
    // Good form to declare a TYPEDEF
    typedef std::vector< std::vector<int> > Tmatrix;
    
    Tmatrix matrix(n+1);
    
    // Size the vectors in the 2.nd dimension. Unfortunately C++ doesn't
    // allow for allocation on declaration of 2.nd dimension of vec of vec
    
    for (int i = 0; i <= n; i++) {
        matrix[i].resize(m+1);
    }
    
    // Step 2
    for (int i = 0; i <= n; i++) {
        matrix[i][0]=i;
    }
    for (int j = 0; j <= m; j++) {
        matrix[0][j]=j;
    }
    
    // Step 3
    for (int i = 1; i <= n; i++) {
        const char s_i = source[i-1];
        
        // Step 4
        for (int j = 1; j <= m; j++) {
            const char t_j = target[j-1];
            
            // Step 5
            int cost;
            if (s_i == t_j) {
                cost = 0;
            }
            else {
                cost = 1;
            }
            
            // Step 6
            const int above = matrix[i-1][j];
            const int left = matrix[i][j-1];
            const int diag = matrix[i-1][j-1];
            int cell = min( above + 1, min(left + 1, diag + cost));
            
            // Step 6A: Cover transposition, in addition to deletion,
            // insertion and substitution. This step is taken from:
            // Berghel, Hal ; Roach, David : "An Extension of Ukkonen's
            // Enhanced Dynamic Programming ASM Algorithm"
            // (http://www.acm.org/~hlb/publications/asm/asm.html)
            
            if (i>2 && j>2) {
                int trans=matrix[i-2][j-2]+1;
                if (source[i-2]!=t_j) trans++;
                if (s_i!=target[j-2]) trans++;
                if (cell>trans) cell = trans;
            }
            
            matrix[i][j]=cell;
        }
    }
    
    // Step 7
    return matrix[n][m];
}


void Parameters::simpleUsage(char ** argv) {
    cout
    << "usage: " << argv[0] << " -s [SEQUENCE] -f [REFERENCE] -t [TARGET] -v [VCF-FILE] > [OUTPUT]" << endl
    << "Please see README at http://www.github.com/denizkural/clia ." << endl
    << "author:   Deniz Kural <denizkural@gmail.com>" << endl;
}

void Parameters::usage(char** argv) {
    cout
    << "usage: " << argv[0] << " -s [SEQUENCE] -f [REFERENCE] -t [TARGET] -v [VCF-FILE] > [OUTPUT]" << endl
    << "Please see README at http://www.github.com/denizkural/clia ." << endl
    << "author:   Deniz Kural <denizkural@gmail.com>" << endl;
}


Parameters::Parameters(int argc, char** argv) {
    
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
    fastq_file = "";            // -q --fastq-file
    fasta_reference = "";                 // -f --fasta-reference
    target = "";                // -t --target
    vcf_file = "";              // -v --vcf-file
    
    outputFile = "";            // -o --bam-output
    
    // operation parameters
    useFile = false;            // -x --use-file
    alignReverse = true;        // -r --reverse-complement
    
    int c; // counter for getopt
    
    static struct option long_options[] =
    {
        {"help", no_argument, 0, 'h'},
        {"sequence", required_argument, 0, 's'},
        {"fastq-file", required_argument, 0, 'q'},
        {"fasta-reference", required_argument, 0, 'f'},
        {"target", required_argument, 0, 't'},
        {"vcf-file", required_argument, 0, 'v'},
        {"output-file", required_argument, 0, 'o'},
        {"use-file", no_argument, 0, 'x'},
        {"reverse-complement", no_argument, 0, 'r'},
        
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
                
            // -v --vcf-file
            case 'v':
                target = optarg;
                break;
                
            // -o --output-file
            case 'o':
                outputFile = optarg;
                break;
                
            // -x --use-file
            case 'x':
                useFile = true;
                break;
                
            // -r --reverse-complement
                alignReverse = false;
                break;

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
    
    if (useFile == true && fastq_file == "") {
        cerr << "Please specify a fastq input file." << endl;
    }
    
}

