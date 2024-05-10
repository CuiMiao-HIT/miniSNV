#ifndef FORM_HPP_
#define FORM_HPP_

#include <fstream>
#include <getopt.h>
#include <string.h>

#include "desc.hpp"

#define PATH_LEN 1024

typedef struct{
    int argc;
    char **argv;
    
// Parameters of input
    char ref_path[PATH_LEN];
    char bam_path[PATH_LEN];
    char tmp_outdir[PATH_LEN];
    char snvvcf_path[PATH_LEN];
    
// Intermediate file
    char p2vrt_infopath[PATH_LEN];
    char p2vrtpos_infopath[PATH_LEN];
    char STEP1vcfInfo_path[PATH_LEN];
    char phasedtable_path[PATH_LEN];
    char highvcf_path[PATH_LEN];

// Other Parameters
    // int ifend = 0;//0:give the whole chr length 1:endpos is given
    char index_path[PATH_LEN];
    char bed_path[PATH_LEN];
    char dup_bed_path[PATH_LEN];
    char chrID[PATH_LEN];
    uint32_t length;
    uint32_t readbam_start;
    uint32_t readbam_end;
    float reads_R;
    int baseQ;
    
}opt_t;

void opt_init(opt_t* P);
int opt_parse(int argc, char *argv[], opt_t* p);
void opt_log(opt_t* p);


#endif /* FORM_HPP_ */