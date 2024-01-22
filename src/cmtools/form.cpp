/*
 * form.cpp
 *
 *  Created on: 2023年08月01日
 *      Author: miaocui
 */
#include "form.hpp"

void opt_init(opt_t *opt) {
    strcpy(opt->index_path, "NULL");
    strcpy(opt->bam_path, "NULL");
    strcpy(opt->snvvcf_path, "./merge.vcf");
    opt->readbam_start = 0;
    opt->reads_R = 0.98;
}

int help_usage() {
    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   %s, variant calling\n", PACKAGE_NAME); 
    fprintf(stderr, "Version:   %s\n", PACKAGE_VERSION); 
    fprintf(stderr, "Contact:   %s\n", CONTACT);
    fprintf(stderr, "Usage:     %s <command> [Options]\n\n", PACKAGE_NAME); 
    fprintf(stderr, "Command: \n");
	fprintf(stderr, "    index      index reference sequence\n");
	fprintf(stderr, "    caller	    call snp and indel\n");
	fprintf(stderr, "\n");

    fprintf(stderr, "Usage:	%s index <ref.fa> <index_route>\n", PACKAGE_NAME);
	fprintf(stderr, "\n\n");
	fprintf(stderr, "Usage:	%s caller [options] -i <index_route> -b <bam> -d <bed> -o <output.vcf>\n\n", PACKAGE_NAME);

    fprintf(stderr, "    <index_route>                 The path of RdBG index.\n");
    fprintf(stderr, "    <bam>                         The input bam format.\n");
    fprintf(stderr, "    <bed>                         The input bed format.\n\n");
    
    fprintf(stderr, "\nOutput options:\n\n");
    fprintf(stderr, "    <output.vcf>                  The output vcf file.\n");
    fprintf(stderr, "    -h --help                      Show detailed usage.\n");
    fprintf(stderr, "\n");
    
    return 1; 
}

int caller_usage() {
    fprintf(stderr, "\n"); 
    fprintf(stderr, "Program:   %s, snv calling\n", PACKAGE_NAME); 
	fprintf(stderr, "Usage:	%s caller [options] -i <index_route> -b <bam> -d <bed> -o <output.vcf>\n\n", PACKAGE_NAME);

    fprintf(stderr, "    <index_route>                 The path of RdBG index.\n");
    fprintf(stderr, "    <bam>                         The input bam format.\n");
    fprintf(stderr, "    <bed>                         The input bed format.\n\n");
    
    fprintf(stderr, "Algorithm options:\n\n");
    fprintf(stderr, "    -c --chr                    Chromosome name to be processed.\n");
    fprintf(stderr, "    -s --start           [UINT]   The starting position of the chromosome to be processed.\n");
    fprintf(stderr, "    -e --end             [UINT]   The ending position of the chromosome to be processed.\n");
    fprintf(stderr, "    -R --read_Ratio      [FLOAT]  The ratio of reads.\n");


    fprintf(stderr, "\nOutput options:\n\n");
    fprintf(stderr, "    -o --fout_vcf           [VCF]    Output file for vatient.\n");
    fprintf(stderr, "    -h --help                      Show detailed usage.\n");
    fprintf(stderr, "\n");
    
    return 1; 
}

char *const short_options = "i:b:d:u:c:s:e:R:t:o:h:p:P:v:g:a";
struct option long_options[] = {
    { "fin_index", 1, NULL, 'i'},
    { "fin_bam", 1, NULL, 'b'},
    { "fin_bed", 1, NULL, 'd'},
    { "dup_bed", 1, NULL, 'u'},
    { "chr", 1, NULL, 'c'},
    { "start", 1, NULL, 's'},
    { "end", 1, NULL, 'e'},
    { "read_Ratio", 1, NULL, 'R'},
    { "tmp_outdir", 1, NULL, 't'},
    { "p2vrt", 1, NULL, 'p'},
    { "p2vrtpos", 1, NULL, 'P'},
    { "s1vcf", 1, NULL, 'v'},
    { "highvcf", 1, NULL, 'g'},
    { "phased_tsv", 1, NULL, 'a'},
    { "fout_vcf", 1, NULL, 'o'},
    { "help", 0, NULL, 'h'},
    { 0, 0, 0, 0}
};

int opt_parse(int argc, char *argv[], opt_t *opt) {
	int optc; 
    char *p;
    int option_index=0;
    if (argc == 2) return caller_usage();
    while((optc = getopt_long(argc, argv, short_options, long_options, &option_index))>=0){ 
        switch(optc){
            case 'i': strcpy(opt->index_path, optarg); break;
            case 'b': strcpy(opt->bam_path, optarg); break;
            case 'd': strcpy(opt->bed_path, optarg); break;
            case 'u': strcpy(opt->dup_bed_path, optarg); break;
            case 't': strcpy(opt->tmp_outdir, optarg); break;
            case 'o': strcpy(opt->snvvcf_path, optarg); break;
            case 'c': strcpy(opt->chrID, optarg);  break;
            case 's': opt->readbam_start = atoi(optarg); break;
            case 'e': opt->readbam_end = atoi(optarg); opt->ifend = 1; break;
            case 'R': opt->reads_R = atof(optarg); break;
            case 'p': strcpy(opt->p2vrt_infopath, optarg); break;
            case 'P': strcpy(opt->p2vrtpos_infopath, optarg); break;
            case 'v': strcpy(opt->STEP1vcfInfo_path, optarg); break;
            case 'g': strcpy(opt->highvcf_path, optarg); break;
            case 'a': strcpy(opt->phasedtable_path, optarg); break;
            case 'h': return caller_usage(); break;
            default:
                fprintf(stderr,"Wrong parameters\n");
                return caller_usage();
        }
    }
    if(optind +2 != argc) {
        return caller_usage();
    }
    opt->argv = argv;
    opt->argc = argc;
    optind++;
    return 0;  
}