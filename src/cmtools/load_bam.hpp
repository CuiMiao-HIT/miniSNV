#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <vector>

extern "C"
{
    #include "../htslib/htslib/sam.h"
    #include "../htslib/htslib/bgzf.h"
    #include "../htslib/htslib/hts.h"
}

using namespace std;

void readFasta_wholechr(char* ref_path, const std::string & chrName, std::string & sequence);
std::string splitFasta(std::string& sequence, uint32_t chrLength, const std::string& chrName, uint32_t startPos, int length);
string getCigar(const bam1_t *b);
char fourbits2base(uint8_t val);
string getSeq(const bam1_t *b);
string getQual(const bam1_t *b);