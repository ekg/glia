#include "utility.h"

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

char shortInt2QualityChar(short i) {
    return static_cast<char>(i + 33);
}

bool allATGC(string& s) {
    for (string::iterator c = s.begin(); c != s.end(); ++c) {
        char b = *c;
        if (b != 'A' && b != 'T' && b != 'G' && b != 'C') {
            return false;
        }
    }
    return true;
}

bool allN(string& s) {
    for (string::iterator c = s.begin(); c != s.end(); ++c) {
        if (*c != 'N') {
            return false;
        }
    }
    return true;
}
