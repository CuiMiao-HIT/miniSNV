/*
 * SNV_code.hpp
 *
 *  Created on: 2023年02月15日
 *      Author: miaocui
 */

#ifndef CALLERSTEP1A3_HPP_
#define CALLERSTEP1A3_HPP_
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <map>
#include <bitset>
#include <cmath>
#include <regex>
#include <cctype>
#include <numeric>
#include <ctime>
#include <cstdio>

#include "../../src/abpoa/abpoa.h"
#include "../../src/ksw2/ksw2.h"
#include "../cmtools/edlib.h"
#include "../cmtools/load_bam.hpp"
#include "../cmtools/math_func.hpp"
#include "../cmtools/mantaAssembler.hpp"
#include "../cmtools/haplotype.hpp"
#include "../cmtools/form.hpp"

extern "C"
{
    #include "../clib/utils.h"
}

using namespace std;

#define ERR_snp 0.04;
// #define ERR_indel 0.1;

typedef struct refid {
    std :: string chrid;
    uint32_t length;
};

typedef struct range_homo {
    uint32_t start;
    uint32_t end;
};

typedef struct READID {
    std :: string readid;
    int stand;
    // int qual;
    bool operator== (const READID &a) const{
        return (readid==a.readid) && (stand == a.stand);
    }
};

typedef struct VRT_info {
    int isknown = 1; //0:konwn 1:unknown and QUA^ 2:unknown and consensus
    // int vrttype;//0:X 1:I|D
    std :: string REF = "";
    std :: string ALT = "";
    int length_in_ref;
    int length_in_read;
    int suppvrt_read = 0;
    int suppstand0 = 0;
    int suppstand16 = 0;
    std :: vector<READID> supp_vrt_reads;//low snv used
    bool operator== (const VRT_info &a) const{
        return (length_in_ref == a.length_in_ref) && (length_in_read == a.length_in_read) && (ALT == a.ALT);
    }
    bool operator< (const VRT_info &a) const{//big->small
        if(suppvrt_read != a.suppvrt_read)
            return suppvrt_read > a.suppvrt_read;
    }
};

typedef struct PosInfo {//skeleton block read info
    // int chrID;
    int ref_num = 0;//real snv pos
    int vrt_num = 0; //how many varients in the pos <candidate used>
    // int max_vrtlengthinref = 1;
    std :: vector<READID> ALLreads;
    std :: vector <VRT_info> Vrts;
};


typedef struct snvVRT {
    uint32_t ref_pos;
    int vrt_num = 0;//reads support the vrt
    int vrtAref_num = 0;//reads support the vrt and reference 
    int all_num = 0;//all reads support this pos
    
    FormatInfo outinfo; 
    std :: string REF;
    std :: string ALT;
    bool operator< (const snvVRT &a) const{
        if(ref_pos != a.ref_pos)
            return ref_pos < a.ref_pos;
    }
    bool operator== (const snvVRT &a) const{
        return (ref_pos==a.ref_pos) && (REF==a.REF) && (ALT==a.ALT);
    }
};

typedef struct lowsnvVRT {
    uint32_t ref_pos;
    std :: string REF;
    std :: string ALT;
    FormatInfo outinfo;
    int all_num;
    int vrt_num;
    int supportnum = 0;
    int HP_num = 0;
    int inhomo = 0;//1:in 0:out
    int vrt_num_slop10 = 0;
    bool operator< (const lowsnvVRT &a) const{
        if(ref_pos != a.ref_pos)
            return ref_pos < a.ref_pos;
    }
    bool operator== (const lowsnvVRT &a) const{
        return (ref_pos==a.ref_pos) && (REF==a.REF) && (ALT==a.ALT);
    }
};

//----------consensus info BEGIN----------

typedef struct P2vrt {
    uint32_t refpos;//real snv pos
    std :: string REF;
    std :: string ALT;
    FormatInfo outinfo;
    int vrt_num = 0;//reads support the vrt
    // int vrtAref_num = 0;//reads support the vrt and reference 
    int all_num = 0;//all reads support this pos 
    int inhomo = 0;
    int vrt_num_slop10 = 0;//0:only the vrt 1:10bp range has more vrts
    bool operator== (const P2vrt &a) const{
        return (refpos==a.refpos) && (REF==a.REF) && (ALT==a.ALT);
    }
};

typedef struct vrtPOS_bam {
    int stand;
    std :: vector<uint32_t> ref_pos_idx;
    bool operator== (const vrtPOS_bam &a) const{
        return (stand == a.stand);
    }
};

typedef struct readstr_name {
    // int chrID;
    string name;
    string str;
    bool operator== (const readstr_name &a) const{
        return (name == a.name);
    }
};

typedef struct CON_suppR{
    string consensus;
    int suppsnv_readnum;
};

typedef struct ConsensusINFO {
    int all_cov_readnum;
    uint32_t ref_pos;//idx pos
    std :: vector <CON_suppR> con_suppr;
    bool operator== (const ConsensusINFO &a) const{
        return ref_pos == a.ref_pos;
    }
};

//----------consensus info END----------

int run_callerStep1(int argc, char *argv[]);
int run_callerStep3(int argc, char *argv[]);

void addreadposin_code(std::string & chr_reference, std :: vector<string> &cigar_sp,std :: string & readname, int stand, uint32_t start_readalnref_pos, std :: string & readseq, std::string & basequal, std :: map<uint32_t , PosInfo> &pos_infos);
void create_candidate_code(uint64_t readbam_start, uint64_t readbam_end, std :: map<uint32_t , PosInfo> & candidata_p1);
void load_bed(char *BED_PATH, char *chrID, std :: vector <range_homo> & homo_ran);
int find_homo(std :: vector <range_homo> homo_ran, uint32_t refpos);
void create_highvrt_p2vrt(std :: map<uint32_t , PosInfo> candidata_p1, std :: vector <snvVRT> &ResultSnv, std :: map <string, vector<vrtPOS_bam>> &p2vrtpos, std :: map<uint32_t, std :: vector <P2vrt>> &p2vrt);
void load_read_part2_bam(std :: map <string, std :: vector<vrtPOS_bam>>  p2vrtpos, std :: map <int,std :: vector <readstr_name>> &clupos_str, uint32_t readbam_start, uint32_t readbam_end);
void mantaAss_m(std :: map <int,std :: vector <readstr_name>>  clupos_str, std :: vector <ConsensusINFO> &consensus_str);
void mantaAss_mSTEP3(std :: map<uint32_t, std :: vector <P2vrt>> p2vrt, std :: map <int,std :: vector <readstr_name>>  clupos_str, std :: map<string,int> read_phased, std :: vector <ConsensusINFO> &consensus_str, std :: vector <lowsnvVRT> &ResultSnv);
void abpoa_consensus(int chrID, std :: map <int,std :: vector <readstr_name>> clupos_str, std :: vector <ConsensusINFO> & consensus_str);
void ksw_score(std :: vector <readstr_name>  reads_str, std :: string  refstr, std :: vector <CON_suppR> &consensus_one);
int ksw_align(std :: string  refstr, std :: string  readstr, int len_tl, int len_ql, int sc_mch, int sc_mis, int gapo, int gape, std :: string &cigar);
int edlib_used(std :: string  qstr, int len_ql, std :: string  tstr, int len_tl);
int edlib_distance(std :: vector <readstr_name>  reads_str, std :: string  refstr, std :: vector <CON_suppR> &consensus_one);
void analyze_K2snv(std :: string refstr, std :: string consensus, std :: string cigar, uint32_t startrefpos, int suppallnum, int suppsnvnum, std :: vector <P2vrt> posp2vrt , std :: vector <lowsnvVRT> &ReaultV_add);
void Re_align(std :: map<uint32_t, std :: vector <P2vrt>> p2vrt, std :: vector <ConsensusINFO> consensus_str, std :: map <int,std :: vector <readstr_name>> clupos_str, std :: vector <lowsnvVRT> &ResultSnv);
int printS1vcf(std :: vector <snvVRT> ResultSnv);
int printS1Resultvcf(char *SNV_PATH, std :: vector <snvVRT> ResultS1snv);
int printALLsnvvcf(std :: vector <lowsnvVRT> &ResultLowsnv);
#endif /* CALLERSTEP1A3_HPP_ */