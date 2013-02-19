gsw-cmd
========

A command line gsw aligner with two modes of use:
- interactive mode: can be piped
- batch mode: accepts a fastq file, a reference, and a vcf
Note that there is a bunch of code for future hashing purposes (currently unused)

Installation Instructions
-------------------------

### Sparsehash:

* Download and install http://code.google.com/p/sparsehash/  using:
    ./configure
    make install

* Make sure that your include paths correspond/include the installation path for sparsehash (usually /usr/local/include)
    This is not the case by default if you are compiling in XCode! You need to go to project build settings and add it to:
    - library search paths
    - header search paths

* Supports latest version. Latest tested version is Sparsehash 2.0.2

* Make sure to use libstdc++ (GNU C++ standard library). In Xcode 4.5 and above, much of tr1 has become standard; thus make sure to edit Build Settings > Apple LLVM compiler 4.2 - language and change the C++ Standard Library to above


### JSON:

* For json support currently json-cpp 0.6.0-rc2 is used, available at http://jsoncpp.sourceforge.net/
    however due to a number of bugs in the original source; a modified version is included here.
    - Copy it to /usr/local/include/json/json.h  or another path you recognize.
    - TODO: Considering switch to a better maintained package (like rapidjson) in the future.

TODO:
-----

* enable command-line options

DONE:
-----

+ have a local .gitignore
+ xcode originating push
+ dependancy docs