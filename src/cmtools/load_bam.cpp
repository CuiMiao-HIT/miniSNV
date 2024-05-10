
#include "load_bam.hpp"

#include <fstream>

using namespace std;

void readFasta_wholechr(char* ref_path, const std::string & chrName, std::string & sequence) {
    ifstream ifref;
    ifref.open(ref_path, ios::in);
    std::string refline;
    bool foundChr = false;
    // int chrLength = 0;

    while (getline(ifref, refline)) {
        if (refline.empty()) continue;
        if (refline[0] == '>') {
            if (foundChr) break;
            stringstream bed_ss(refline);
            string sp1;
            bed_ss >> sp1;
            if (sp1.substr(1) == chrName) {
                foundChr = true;
            }
        } else if (foundChr) {
            sequence += refline;
            // chrLength += refline.length();
        }
    }
    ifref.close();
    if (!foundChr) {
        std::cerr << "Chromosome [" << chrName << "] not found in the file." << std::endl;
    }
}

std::string splitFasta(std::string& sequence, uint32_t chrLength, const std::string& chrName, uint32_t startPos, int length) {
    uint32_t startPos_next = startPos - 1;
    if (startPos_next >= chrLength) {
        std::cerr << "Invalid start position or length." << std::endl;
        return "";
    }
    else if(startPos_next + length >= chrLength){
        int length_end = chrLength-startPos_next;
        return sequence.substr(startPos_next, length_end);
    }
    return sequence.substr(startPos_next, length);
}

string getCigar(const bam1_t *b) {
    uint32_t *data = (uint32_t *)bam_get_cigar(b);
    int cigarNum = b->core.n_cigar;
    stringstream ss;
    for(int i=0; i<cigarNum; i++) {
        uint32_t val = data[i];
        char op = bam_cigar_opchr(val);
        uint32_t len = bam_cigar_oplen(val);
        ss << len << op;
    }
    return ss.str();
}

//Each base is encoded in 4 bits: 1 for A, 2 for C, 4 for G, 8 for T and 15 for N.
char fourbits2base(uint8_t val) {
    switch(val) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            // cerr << "ERROR: Wrong base with value "<< (int)val << endl ;
            return 'N';
    }
}

string getSeq(const bam1_t *b) {
    uint8_t *data = bam_get_seq(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        char base;
        if(i%2 == 1)
            base = fourbits2base(data[i/2] & 0xF); 
        else
            base = fourbits2base((data[i/2]>>4) & 0xF);
        s[i] = base;
    }
    return s;
}

// get seq quality
string getQual(const bam1_t *b) {
    uint8_t *data = bam_get_qual(b);
    int len = b->core.l_qseq;
    string s(len, '\0');
    for(int i=0; i<len; i++) {
        s[i] = (char)(data[i] + 33); // 转换成打印的ascci
    }
    return s;
}