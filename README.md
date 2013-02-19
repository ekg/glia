clia-cmd
====

A command line gsw aligner with two modes of use:

- interactive mode: can be piped
- batch mode: accepts a fastq file, a reference, and a vcf

Note that there is a bunch of code for future hashing purposes (currently unused)

Installation Instructions
====

Sparsehash:

* Download and install http://code.google.com/p/sparsehash/  using:
    ./configure
    make install

* Make sure that your include paths correspond/include the installation path for sparsehash (usually /usr/local/include)
  This is not the case by default if you are compiling in XCode! You need to go to project build settings and add it to:
    - library search paths
    - header search paths

* Supports latest version. Latest tested version is Sparsehash 2.0.2

* In Xcode 4.5 and above, much of tr1 has become standard; thus in your install path (/usr/local/include); 
  find sparseconfig.h & edit:

    #define HASH_FUN_H <tr1/functional> ->  #define HASH_FUN_H <functional>
    #define HASH_NAMESPACE std::tr1     ->  #define HASH_NAMESPACE std




JSON:

* For json support currently json-cpp 0.6.0-rc2 is used, available at http://jsoncpp.sourceforge.net/
  however due to a number of bugs in the original source; a modified version is included here - no install necessary.
  Considering switch to a better maintained package (like rapidjson) in the future.



TODO:

* enable command-line options

DONE:

+ have a local .gitignore
+ xcode originating push
+ dependancy docs

