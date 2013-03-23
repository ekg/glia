#include "utility.h"

short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

char shortInt2QualityChar(short i) {
    return static_cast<char>(i + 33);
}
