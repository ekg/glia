#include <string>
#include <sstream>
#include "gssw.h"

using namespace std;

short qualityChar2ShortInt(char c);
char shortInt2QualityChar(short i);
bool allATGC(string& s);
bool allN(string& s);
string graph_mapping_to_string(gssw_graph_mapping* gm);
