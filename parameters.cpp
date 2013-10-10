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
	<< "Please see README at http://www.github.com/denizkural/clia" << endl
        << "Use --help to read detailed options." << endl
	<< "authors:   Deniz Kural <denizkural@gmail.com>" << endl
	<< "           Erik Garrison <erik.garrison@gmail.com>" << endl;
}

void Parameters::usage(char** argv) {
    cout
        << "usage: " << argv[0] << " -s [SEQUENCE] -f [REFERENCE] -t [TARGET] -v [VCF-FILE] > [OUTPUT]" << endl
        << "Please see README at http://www.github.com/denizkural/clia" << endl
        << "authors:   Deniz Kural <denizkural@gmail.com>" << endl
        << "           Erik Garrison <erik.garrison@gmail.com>" << endl
        << endl
        << "general parameters:" << endl
        << endl
        << "    -h --help                  This dialog." << endl
        //<< "    -q --fastq-file FILE       The fastq file from which to draw reads" << endl
        << "    -f --fasta-reference FILE  The reference sequence for alignment." << endl
        << "    -v --vcf-file FILE         The genome DAG, BGZIP'ed (.vcf.gz) and tabix-indexed (.tbi)" << endl
        //<< "    -o --output-file FILE      Write alignments in BAM format to FILE." << endl
        << "    -m --match N               The alignment match score (integer, default 10)." << endl
        << "    -M --mismatch N            The alignment mismatch score (integer, default -10)." << endl
        << "    -g --gap N                 The alignment gap score (integer, default -10)." << endl
        << endl
        << "single-read alignment:" << endl
        << endl
        << "    -s --sequence SEQ          The sequence to align." << endl
        << "    -t --target TARGET         Target genomic region for alignment, e.g. chr2:1-20" << endl
        << "    -D --display-dag           Write DAG generated from variants to stdout." << endl
        << "    -B --display-backtrace     Write alignment matrix results to stdout." << endl
        << "    -N --display-all-nodes     Same as -B but also for nodes which are not traced." << endl
        //<< "    -P --display-alignment     Print sequence from DAG and read sequence." << endl
        << "    -d --debug                 Enable debugging output." << endl
        << "    -r --reverse-complement    Report the reverse complement if it provides a better alignment." << endl
        << endl
        << "local realignment:" << endl
        << endl
        << "    -R --realign-bam           Realign the BAM stream on stdin to the VCF file, adjusting" << endl
        << "                               position and flattening alignments back into the reference space" << endl
        << "                               where realignment to the graph provides better quality." << endl
/*
          << "    -F --flatten-alignments    Flatten output into reference-relative alignments where the DAG" << endl
          << "                               alignment provides a better match than the reference.  Use this" << endl
          << "                               if you want to use alignment-based variant detectors on the output." << endl
        */
        << "    -w --realignment-window    Number of bp of window to assemble from VCF for realignment." << endl
        << "    -S --soft-clip-qsum-threshold N"  << endl
        << "                               If sum of qualities of soft clipped bases is > N, realign." << endl
        << "    -Q --mismatch-qsum-threshold N" << endl
        << "                               If sum of qualities of mismatched bases is > N, realign." << endl
        << "    -Z --soft-clip-qsum-max N  Accept realignment if qsum of softclips is < N." << endl
        << "    -E --mismatch-qsum-max N   Accept realignment if qsom of mismatches is < N." << endl
        << "    -G --gap-count-max N       Accept realignment if number of gaps is < N." << endl
        << "    -X --dry-run               If realigning, don't output BAM (helps for debugging)." << endl
        << "    -O --only-realigned        Emit only realigned records (debugging)." << endl
        << "    -n --flatten-flank N       Use this many bp of dummy sequence when flattening large" << endl
        << "                               insertions into reference coordinates (default: 2)" << endl;
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
    alignReverse = false;        // -r --reverse-complement

    dag_window_size = 1500;

    match = 10;
    mism = -10;
    gap = -10;

    display_dag = false;
    display_backtrace = false;
    display_all_nodes = false;
    display_alignment = false;

    debug = false;
    dry_run = false;
    only_realigned = false;

    realign_bam = false;
    softclip_qsum_threshold = 10000000;
    mismatch_qsum_threshold = 10000000;
    softclip_qsum_max = 10000000;
    mismatch_qsum_max = 10000000;
    gap_count_max = 10000000;

    flatten_alignments = true;
    flatten_flank = 2;
    
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
        {"gap", required_argument, 0, 'g'},
        {"match", required_argument, 0, 'm'},
        {"mismatch", required_argument, 0, 'M'},
        {"display-backtrace", no_argument, 0, 'B'},
        {"display-all-nodes", no_argument, 0, 'N'},
        {"display-dag", no_argument, 0, 'D'},
        {"debug", no_argument, 0, 'd'},
        {"dry-run", no_argument, 0, 'X'},
        {"realign-bam", no_argument, 0, 'R'},
        {"soft-clip-qsum-threshold", required_argument, 0, 'S'},
        {"mismatch-qsum-threshold", required_argument, 0, 'Q'},
        {"soft-clip-qsum-max", required_argument, 0, 'Z'},
        {"mismatch-qsum-max", required_argument, 0, 'E'},
        {"gap-count-max", required_argument, 0, 'G'},
        {"realignment-window", required_argument, 0, 'w'},
        {"flatten-alignments", no_argument, 0, 'F'},
        {"only-realigned", no_argument, 0, 'O'},
        {"flatten-flank", required_argument, 0, 'n'},
        {0, 0, 0, 0}
        
    };
    
    while (true) {
        int option_index = 0;
        c = getopt_long(argc, argv,
                        "hxrXONBDFRds:q:f:t:v:o:g:m:M:w:Q:S:Z:E:G:n:",
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
            vcf_file = optarg;
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
	    case 'r':
            alignReverse = true;
            break;

            // -N --display-all-nodes
	    case 'N':
            display_all_nodes = true;
            break;

            // -B --display-backtrace
	    case 'B':
            display_backtrace = true;
            break;

            // -X --dry-run
	    case 'X':
            dry_run = true;
            break;

        case 'O':
            only_realigned = true;
            break;

            // -D --display-dag
	    case 'D':
            display_dag = true;
            break;

            // -d --debug
	    case 'd':
            debug = true;
            break;

            // -m --match
        case 'm':
            match = atoi(optarg);
            break;
                
            // -M --mismatch
        case 'M':
	        mism = atoi(optarg);
            break;
                
            // -g --gap
        case 'g':
            gap = atoi(optarg); 
            break;

            // -R --realign-bam
        case 'R':
            realign_bam = true;
            break;

        case 'F':
            flatten_alignments = true;
            break;

        case 'n':
            flatten_flank = atoi(optarg);
            break;

        case 'Q':
            mismatch_qsum_threshold = atoi(optarg);
            break;

        case 'S':
            softclip_qsum_threshold = atoi(optarg);
            break;

        case 'Z':
            softclip_qsum_max = atoi(optarg);
            break;

        case 'E':
            mismatch_qsum_max = atoi(optarg);
            break;

        case 'G':
            gap_count_max = atoi(optarg);
            break;

            // -w --realignment-window
        case 'w':
            dag_window_size = atoi(optarg);
            break;
                
            // -h --help
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

