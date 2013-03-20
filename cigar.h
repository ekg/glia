#ifndef CIGAR_H
#define CIGAR_H

#include <iostream>
#include <utility>
#include <sstream>
#include <vector>
#include "vcflib/Variant.h"
#include "api/BamAux.h"  // CigarOp for conversion

using namespace std;
using namespace vcf;

struct CigarElement {
    int length;
    char type;
    void clear(void);
CigarElement() : length(0), type('M') { }
CigarElement(int l, char t) : length(l), type(t) { }
};

struct Cigar : vector<CigarElement> {
    void append(const Cigar& c);
    int refLen(void);
    bool isReference(void);
    string str(void);
    Cigar(void) { }
    Cigar(int, char);
    Cigar(const string& cigarStr);
    Cigar(vector<VariantAllele>& vav);
    Cigar(vector<BamTools::CigarOp>& cigarData);
    void toCigarData(vector<BamTools::CigarOp>& cigarData);
};

std::ostream& operator<<(std::ostream& o, const CigarElement& e);
std::ostream& operator<<(std::ostream& o, const Cigar& c);
Cigar join(std::vector<Cigar>& cigars);

#endif
