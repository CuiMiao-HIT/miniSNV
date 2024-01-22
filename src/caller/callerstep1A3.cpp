/*
 * callerstep1.cpp
 *
 *  Created on: 2023年02月15日
 *      Author: miaocui
 */

#include "callerstep1A3.hpp"

opt_t *opt;
std :: map < std::string , int > chrID2int;
std :: map <int , refid> int2chrID;
int give_chrint = -1;

// AaCcGgTtNn ==> 0,1,2,3,4
unsigned char nt4_table[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4 /*'-'*/, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T, else=>N
const char nt256_table[256] = {
    'A', 'C', 'G', 'T', 'N', '-', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', '-', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'T', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'A', 'N', 'C', 'N', 'N', 'N', 'G', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'T', 'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};

int run_callerStep1(int argc, char *argv[]){
    double cpu_time_all = cputime_hp();
    struct timeval start_all;
    gettimeofday(&start_all, NULL);

    opt_t opt_ent;   
    opt = &opt_ent; //global opt
    opt_init(opt);
    if (opt_parse(argc, argv, opt) != 0)  return 2;
    // if(strcmp(opt->index_path,"NULL")==0){
    //     fprintf(stderr, "PLEASE GIVE THE INDEX PATH!    -i <index_route>\n");
    //     return 2;
    // } if(strcmp(opt->bam_path,"NULL")==0){
    //     fprintf(stderr, "PLEASE GIVE THE BAM PATH!  -b <bam_path>\n");
    //     return 2;
    // } if(strcmp(opt->bed_path,"NULL")==0){
    //     fprintf(stderr, "PLEASE GIVE THE BED PATH!  -d <bed_path>\n");
    //     return 2;
    // }if(strcmp(opt->tmp_outdir,"NULL")==0){
    //     fprintf(stderr, "PLEASE GIVE THE TMP_DIR PATH!  -d <tmp_dir_path>\n");
    //     return 2;
    // }

    /*save all chrid name and chrid number and chr length*/
    Simple_ref_handler refID;
    refID.load_bin_ref(opt->index_path);
    int ref_N = refID.get_chr_N();
    for (int refi = 0; refi < ref_N; refi++) {
        std :: string chr_name = refID.get_chr_name(refi);
        uint32_t chr_length = refID.get_chr_length(refi);
        chrID2int[chr_name] = refi;
        refid refinfo;
        refinfo.chrid = chr_name; refinfo.length = chr_length;
        int2chrID[refi] = refinfo;
    }
    std :: map < std::string , int > :: iterator itfindchrid = chrID2int.find(opt->chrID);
    if(itfindchrid != chrID2int.end()) give_chrint = itfindchrid->second; 
    if (give_chrint == -1) {
        fprintf(stderr, "CANNOT FIND %s IN BAM FILE!\n",opt->chrID);
        return 2;
    }
    if(opt->ifend == 0) opt->readbam_end = int2chrID[give_chrint].length;

    // result
    std :: vector <snvVRT> ResultSnv;

    // std :: clock_t c_start1 = std :: clock();
    std :: map<uint32_t , PosInfo> candidata_p1;//real_refpos , posINFO : num&&vrt
    create_candidate_code(opt->dup_bed_path, opt->bam_path, opt->index_path, give_chrint, opt->readbam_start, opt->readbam_end, candidata_p1);
    // std::clock_t c_end1 = std::clock();
    // double time_elapsed_ms1 = 1000.0 * (c_end1-c_start1) / CLOCKS_PER_SEC;
    // fprintf(stderr, "\ncreat_candidate_code CPU : %.5f sec\n", time_elapsed_ms1 / 1000.0);
    
    if (candidata_p1.empty()){
        fprintf(stderr, "[%s_%d] candidate file is empty\n\n", opt->chrID, opt->readbam_start);
        return 2;
    }

    // std :: clock_t c_start2 = std :: clock();
    std :: map<uint32_t, std :: vector <P2vrt>> p2vrt;//idx_refpos , vrt
    std :: map <string, std :: vector<vrtPOS_bam>> p2vrtpos;//readID , idx_refpos
    create_highvrt_p2vrt(opt->bed_path, candidata_p1, ResultSnv, p2vrtpos, p2vrt);
    candidata_p1.clear();
    // std::clock_t c_end2 = std::clock();
    // double time_elapsed_ms2 = 1000.0 * (c_end2-c_start2) / CLOCKS_PER_SEC;
    // fprintf(stderr, "\ncreat_highvrt_p2vrt CPU : %.5f sec\n", time_elapsed_ms2 / 1000.0);
    
    if (p2vrt.empty()){
        int Pvcf = printS1vcf(opt->STEP1vcfInfo_path, opt->highvcf_path, ResultSnv);
        ResultSnv.clear();
        fprintf(stderr, "[%s_%d] low varient file is empty\n\n", opt->chrID, opt->readbam_start);
        return 2;
    }
    ofstream p2vrt_info;
    p2vrt_info.open(opt->p2vrt_infopath);
    std :: map<uint32_t, std :: vector <P2vrt>> :: iterator p2vrt_it = p2vrt.begin();
    while (p2vrt_it != p2vrt.end()){
        p2vrt_info << p2vrt_it->first << endl;
        for(auto aa : p2vrt_it->second) p2vrt_info << aa.refpos  << "_" << aa.REF << "_" << aa.ALT << "_" << aa.outinfo.GQ<< "_" << aa.outinfo.GT << "_" << aa.outinfo.QUAL << "_" << aa.all_num << "_" << aa.vrt_num << "_" << aa.inhomo << "_" << aa.vrt_num_slop10 << "\t";
        p2vrt_info << "\n";
        p2vrt_it++;
    }
    p2vrt_info.close();

    ofstream p2vrtpos_info;
    p2vrtpos_info.open(opt->p2vrtpos_infopath);
    std :: map <string, std :: vector<vrtPOS_bam>>:: iterator p2vrtpos_it = p2vrtpos.begin();
    while (p2vrtpos_it != p2vrtpos.end()){
        p2vrtpos_info << "#" << p2vrtpos_it->first << endl;
        for(auto aa : p2vrtpos_it->second){
            p2vrtpos_info << aa.stand  << "\t";
            for(auto oneS : aa.ref_pos_idx) {
                p2vrtpos_info << oneS << "\t";
            }
            p2vrtpos_info << "\n";
        }
        p2vrtpos_it++;
    }
    p2vrtpos_info.close();

    int Pvcf = printS1vcf(opt->STEP1vcfInfo_path, opt->highvcf_path, ResultSnv);
    ResultSnv.clear();

    fprintf(stderr, "Classify CPU STEP1 of [%s_%d] : %.3f sec\n", opt->chrID, opt->readbam_start, cputime_hp() - cpu_time_all);
    return 0;
}


int run_callerStep3(int argc, char *argv[]){

    double cpu_time_all = cputime_hp();
    struct timeval start_all;
    gettimeofday(&start_all, NULL);

    opt_t opt_ent;   
    opt = &opt_ent; //global opt
    opt_init(opt);
    if (opt_parse(argc, argv, opt) != 0)  return 2;
    // if(strcmp(opt->index_path,"NULL")==0){
    //     fprintf(stderr, "PLEASE GIVE THE INDEX PATH!    -i <index_route>\n");
    //     return 2;
    // } if(strcmp(opt->bam_path,"NULL")==0){
    //     fprintf(stderr, "PLEASE GIVE THE BAM PATH!  -b <bam_path>\n");
    //     return 2;
    // } if(strcmp(opt->bed_path,"NULL")==0){
    //     fprintf(stderr, "PLEASE GIVE THE BED PATH!  -d <bed_path>\n");
    //     return 2;
    // } if(strcmp(opt->tmp_outdir,"NULL")==0){
    //     fprintf(stderr, "PLEASE GIVE THE TMP_DIR PATH!  -d <tmp_dir_path>\n");
    //     return 2;
    // }

    /*save all chrid name and chrid number and chr length*/
    Simple_ref_handler refID;
    refID.load_bin_ref(opt->index_path);
    int ref_N = refID.get_chr_N();
    for (int refi = 0; refi < ref_N; refi++) {
        std :: string chr_name = refID.get_chr_name(refi);
        uint32_t chr_length = refID.get_chr_length(refi);
        chrID2int[chr_name] = refi;
        refid refinfo;
        refinfo.chrid = chr_name; refinfo.length = chr_length;
        int2chrID[refi] = refinfo;
    }
    std :: map < std::string , int > :: iterator itfindchrid = chrID2int.find(opt->chrID);
    if(itfindchrid != chrID2int.end()) give_chrint = itfindchrid->second; 
    if (give_chrint == -1) {
        fprintf(stderr, "CANNOT FIND %s IN BAM FILE!\n",opt->chrID);
        return 2;
    }
    if(opt->ifend == 0) opt->readbam_end = int2chrID[give_chrint].length;

    std :: map<string,int> read_phased;//0:none 1/2:hp1/hp2
    std :: map<uint32_t, std :: vector <P2vrt>> p2vrt;//idx_refpos , vrt
    std :: map <string, std :: vector<vrtPOS_bam>> p2vrtpos;//readID , idx_refpos

    std :: vector <lowsnvVRT> ResultLowsnv;

    /*
    load phased read info
    */
    ifstream ifphased;
    ifphased.open(opt->phasedtable_path, ios::in);
    std :: string phasedline;
    while (getline(ifphased , phasedline)) {
        if(phasedline[0] == '#') continue;
        string readname;
        int phase_flag = 0;
        stringstream bed_ss(phasedline);
        string sp1;
        bed_ss >> readname;
        bed_ss >> sp1;
        if(sp1 == "H1") phase_flag = 1;
        else if(sp1 == "H2") phase_flag = 2;
        read_phased[readname] = phase_flag;
    }
    ifphased.close();
    /*
    load p2vrt info from step1 
    std :: map<uint32_t, std :: vector <P2vrt>> p2vrt     idx_refpos , vrt 
    */
    ifstream ifp2vrtinfo;
    ifp2vrtinfo.open(opt->p2vrt_infopath, ios::in);
    std :: string p2vrtline;
    while (getline(ifp2vrtinfo , p2vrtline)) {
        stringstream bed_ss1(p2vrtline);
        string sp1;
        bed_ss1 >> sp1;
        uint32_t idx_refpos = stoi(sp1);

        P2vrt add_one;
        getline(ifp2vrtinfo , p2vrtline);
        stringstream bed_ss2(p2vrtline);
        string sp2;
        while(bed_ss2 >> sp2){
            std :: vector<string> SP = split2(sp2,"_");
            add_one.refpos = stoi(SP[0]);
            add_one.REF = SP[1];
            add_one.ALT = SP[2];
            add_one.outinfo.GQ = stoi(SP[3]);
            add_one.outinfo.GT = SP[4];
            add_one.outinfo.QUAL = stoi(SP[5]);
            add_one.all_num = stoi(SP[6]);
            add_one.vrt_num = stoi(SP[7]);
            add_one.inhomo = stoi(SP[8]);
            add_one.vrt_num_slop10 = stoi(SP[9]);
            p2vrt[idx_refpos].push_back(add_one);
        }
    }
    ifp2vrtinfo.close();
    /*
    load p2vrtpos info from step1 
    std :: map <string, std :: vector<vrtPOS_bam>> p2vrtpos     readID , idx_refpos
    */
    string readname;
    ifstream ifp2vrtposinfo;
    ifp2vrtposinfo.open(opt->p2vrtpos_infopath, ios::in);
    std :: string p2vrtposline;
    while (getline(ifp2vrtposinfo , p2vrtposline)) {
        if(p2vrtposline[0] == '#'){
            readname = p2vrtposline.substr(1);
            continue;
        } else{
            vrtPOS_bam add_one;
            stringstream bed_ss1(p2vrtposline);
            string sp1;
            bed_ss1 >> sp1;
            add_one.stand = stoi(sp1);
            while(bed_ss1 >> sp1){
                add_one.ref_pos_idx.push_back(stoi(sp1));
            }
            p2vrtpos[readname].push_back(add_one);
        }
    }
    ifp2vrtposinfo.close();

    /*     
    consensus code     
    */
    std :: map <int,std :: vector <readstr_name>> clupos_str; // <start_idx,strings>
    load_read_part2_bam(opt->bam_path, p2vrtpos, clupos_str, give_chrint, opt->readbam_start, opt->readbam_end);
    p2vrtpos.clear();

    std :: vector <ConsensusINFO> consensus_str;
    mantaAss_mSTEP3(p2vrt, clupos_str, read_phased, consensus_str, ResultLowsnv);

    Re_align(p2vrt, consensus_str, clupos_str, ResultLowsnv);
    clupos_str.clear(); consensus_str.clear();
    /*     
    write all vcf file code     
    */
    int Pvcf = printALLsnvvcf(opt->snvvcf_path, opt->STEP1vcfInfo_path, ResultLowsnv);
    fprintf(stderr, "\nClassify CPU STEP3 of [%s] : %.3f sec\n\n", opt->chrID, cputime_hp() - cpu_time_all);
    return 0;
}

/*
--------------------------------------Dividing line----------------------------------------
*/

void addreadposin_code(Simple_ref_handler & ref, std :: vector<string> &cigar_sp, int chrID, std :: string & readname, int stand, uint32_t start_readalnref_pos, std :: string & readseq, std::string & basequal, std :: map<uint32_t , PosInfo> &pos_infos){
    READID addread;
    addread.readid = readname;
    addread.stand = stand;
    int basequal_int;
    uint32_t current_refpos = start_readalnref_pos;
    uint32_t current_readpos = 0;
    
    int temp_size , change_size;
    int if_preM = 0;
    for(auto temp : cigar_sp){
        temp_size = temp.size();
        change_size = stoi(temp.substr(0,temp_size-1));
        if (temp[temp_size-1] == 'S'){
            current_readpos = current_readpos + change_size;
        } else if (temp[temp_size - 1] == 'M'){
            std :: string refpart;
            ref.load_ref_from_buff(chrID, current_refpos, change_size, refpart);
            std :: string readpart = readseq.substr(current_readpos, change_size);

            char Aref,Aread;
            for (int i = 0; i < change_size; i++) {
                char base1[1];
                basequal.copy(base1,1,current_readpos);
                basequal_int = (int)base1[0] - 33 ;
                if(basequal_int >= 13){
                    Aref = refpart[i];
                    Aread = std::toupper(readpart[i]);
                    // uint32_t xpos = current_readpos+i;
                    // Aread = std::toupper(readseq[xpos]);
                    if(Aref == Aread){
                        //=
                        std :: map<uint32_t , PosInfo> :: iterator itArp = pos_infos.find(current_refpos);
                        if(itArp != pos_infos.end()){
                            itArp->second.ref_num = itArp->second.ref_num + 1;
                            itArp->second.ALLreads.push_back(addread);
                        } else{
                            PosInfo addpos1;
                            addpos1.ref_num = 1;
                            addpos1.ALLreads.push_back(addread);
                            pos_infos[current_refpos] = addpos1;
                        }
                    } else{
                        //X
                        VRT_info addsnp;
                        addsnp.suppvrt_read = 1;
                        addsnp.ALT = Aread;
                        addsnp.length_in_ref = 1;
                        addsnp.length_in_read = 1;
                        addsnp.supp_vrt_reads.push_back(addread);
                        if(stand == 0) addsnp.suppstand0 = addsnp.suppstand0 + 1;
                        else addsnp.suppstand16 = addsnp.suppstand16 + 1;
                        std :: map<uint32_t , PosInfo> :: iterator itArp = pos_infos.find(current_refpos);
                        if(itArp != pos_infos.end()){
                            itArp->second.ALLreads.push_back(addread);
                            std :: vector <VRT_info> :: iterator itsnpfind = find(itArp->second.Vrts.begin() , itArp->second.Vrts.end() , addsnp);
                            if(itsnpfind != itArp->second.Vrts.end()){
                                itsnpfind->supp_vrt_reads.push_back(addread);
                                itsnpfind->suppvrt_read = itsnpfind->suppvrt_read + 1;
                                if(stand == 0) itsnpfind->suppstand0 = itsnpfind->suppstand0 + 1;
                                else itsnpfind->suppstand16 = itsnpfind->suppstand16 + 1;
                            } else{
                                itArp->second.Vrts.push_back(addsnp);
                            }
                        } else{
                            PosInfo addpos1;
                            addpos1.Vrts.push_back(addsnp);
                            addpos1.ALLreads.push_back(addread);
                            pos_infos[current_refpos] = addpos1;
                        }
                    }
                }
                current_refpos = current_refpos + 1;
                current_readpos = current_readpos + 1;
            }
            if(Aref == Aread) if_preM = 1;
        } else if (temp[temp_size - 1] == '='){
            int basequal_int;
            for(int i = 0 ; i < change_size ; i++){
                char base1[1];
                basequal.copy(base1,1,current_readpos);
                basequal_int = (int)base1[0] - 33 ;
                if(basequal_int >= 13){
                    std :: map<uint32_t , PosInfo> :: iterator itArp = pos_infos.find(current_refpos);
                    if(itArp != pos_infos.end()){
                        itArp->second.ref_num = itArp->second.ref_num + 1;
                        itArp->second.ALLreads.push_back(addread);
                    } else{
                        PosInfo addpos1;
                        addpos1.ref_num = 1;
                        addpos1.ALLreads.push_back(addread);
                        pos_infos[current_refpos] = addpos1;
                    }
                }
                current_refpos = current_refpos + 1;
                current_readpos = current_readpos + 1;
            }
            if_preM = 1;
        } else if (temp[temp_size - 1] == 'X'){
            for (int mis_i = 0; mis_i < change_size; mis_i++){
                int basequal_int;
                char base1[1];
                basequal.copy(base1,1,current_readpos);
                basequal_int = (int)base1[0] - 33 ;
                if(basequal_int >= 13){
                    // std :: string alt = readseq.substr(current_readpos, 1);
                    // std :: transform(alt.begin(), alt.end(), alt.begin(), ::toupper);
                    char alt = std:: toupper(readseq[current_readpos]);
                    VRT_info addsnp;
                    addsnp.suppvrt_read = 1;
                    addsnp.ALT = alt;
                    addsnp.length_in_ref = 1;
                    addsnp.length_in_read = 1;
                    addsnp.supp_vrt_reads.push_back(addread);
                    if(stand == 0) addsnp.suppstand0 = addsnp.suppstand0 + 1;
                    else addsnp.suppstand16 = addsnp.suppstand16 + 1;
                    std :: map<uint32_t , PosInfo> :: iterator itArp = pos_infos.find(current_refpos);
                    if(itArp != pos_infos.end()){
                        itArp->second.ALLreads.push_back(addread);
                        std :: vector <VRT_info> :: iterator itsnpfind = find(itArp->second.Vrts.begin() , itArp->second.Vrts.end() , addsnp);
                        if(itsnpfind != itArp->second.Vrts.end()){
                            itsnpfind->supp_vrt_reads.push_back(addread);
                            itsnpfind->suppvrt_read = itsnpfind->suppvrt_read + 1;
                            if(stand == 0) itsnpfind->suppstand0 = itsnpfind->suppstand0 + 1;
                            else itsnpfind->suppstand16 = itsnpfind->suppstand16 + 1;
                        } else{
                            itArp->second.Vrts.push_back(addsnp);
                        }
                    } else{
                        PosInfo addpos1;
                        addpos1.Vrts.push_back(addsnp);
                        addpos1.ALLreads.push_back(addread);
                        pos_infos[current_refpos] = addpos1;
                    }
                }
                current_refpos = current_refpos + 1;
                current_readpos = current_readpos + 1;
            }
            if_preM = 0;
        } else if (temp[temp_size - 1] == 'I'){
            uint32_t vrtrefpos_in = current_refpos - 1;
            std :: map<uint32_t , PosInfo> :: iterator itArp = pos_infos.find(vrtrefpos_in);
            if(itArp != pos_infos.end()){
                std :: vector<READID> :: iterator itALLread = find(itArp->second.ALLreads.begin() , itArp->second.ALLreads.end() , addread);
                if(itALLread != itArp->second.ALLreads.end() && if_preM == 1)  itArp->second.ref_num = itArp->second.ref_num - 1;
            }
            if_preM = 0;
            current_readpos = current_readpos + change_size;
        } else if (temp[temp_size - 1] == 'D'){
            uint32_t vrtrefpos_del = current_refpos - 1;
            std :: map<uint32_t , PosInfo> :: iterator itArp = pos_infos.find(vrtrefpos_del);
            if(itArp != pos_infos.end()){
                std :: vector<READID> :: iterator itALLread = find(itArp->second.ALLreads.begin() , itArp->second.ALLreads.end() , addread);
                if(itALLread != itArp->second.ALLreads.end() && if_preM == 1) itArp->second.ref_num = itArp->second.ref_num - 1;
            }
            for(int deladdread = 0; deladdread < change_size; deladdread++){
                uint32_t delotherpos = current_refpos + deladdread;
                std :: map<uint32_t , PosInfo> :: iterator itArp = pos_infos.find(delotherpos);
                if(itArp != pos_infos.end()) itArp->second.ALLreads.push_back(addread);
                else {
                    PosInfo addpos1;
                    addpos1.ALLreads.push_back(addread);
                    pos_infos[delotherpos] = addpos1;
                }
            }
            if_preM = 0;
            current_refpos = current_refpos + change_size;
        }
    }
}

void load_bed(char *BED_PATH, char *chrID, std :: vector <range_homo> & homo_ran){
    ifstream ifbed;
    ifbed.open(BED_PATH, ios::in);
    std :: string bedline;
    while (getline(ifbed , bedline)) {
        stringstream bed_ss(bedline);
        string sp1;
        uint32_t spos , epos;
        bed_ss >> sp1;
        if(sp1 == chrID){
            range_homo push1;
            bed_ss >> sp1;
            push1.start = stoi(sp1)+5;
            bed_ss >> sp1;
            push1.end = stoi(sp1)-4;
            homo_ran.push_back(push1);
            continue;
        }
        else if(homo_ran.size()!= 0) break;
    }
    ifbed.close();
}

int find_homo(std :: vector <range_homo> homo_ran, uint32_t refpos){
    int mid, low = 0, high = homo_ran.size()-1;
    if(refpos<homo_ran[low].start || refpos>homo_ran[high].end) return 0;
    while (low<=high) {
        if((refpos>=homo_ran[low].start && refpos<homo_ran[low].end) || (refpos>=homo_ran[high].start && refpos<homo_ran[high].end)) return 1;

        mid = (low+high)/2;
        if(refpos>=homo_ran[mid].start && refpos<homo_ran[mid].end) return 1;
        if (refpos < homo_ran[mid].start) high = mid-1;
        else if(refpos >= homo_ran[mid].end) low = mid+1;
    }
    return 0;
}

void create_candidate_code(char *DUPBED_PATH, char *BAM_PATH, char *INDEX_PATH, int give_chrid, uint64_t readbam_start, uint64_t readbam_end, std :: map<uint32_t , PosInfo> & candidata_p1) {// analize bam file

    std :: vector <range_homo> dup_ran;
    load_bed(DUPBED_PATH, opt->chrID, dup_ran);
    std :: map <uint32_t , int > secondaryReads;
    std ::map<uint32_t, PosInfo> pos_infos; // real_refpos , posinfo
    Simple_ref_handler ref;
    ref.load_bin_ref(INDEX_PATH);

    // std :: clock_t c_start1 = std :: clock();
    char * region = (char *)malloc(1024);
    snprintf(region, 1024, "%s:%d-%d", opt->chrID, opt->readbam_start, opt->readbam_end);

    samFile *bam_file = sam_open(BAM_PATH, "r");
    samFile *bam_in = sam_open(BAM_PATH, "r"); // open bam file
    bam_hdr_t *bam_header = sam_hdr_read(bam_file);
    hts_idx_t *idxt = hts_idx_load(BAM_PATH, 1);
    hts_itr_t *itrt = bam_itr_querys(idxt, bam_header, region);
    bam1_t *readinbam = bam_init1();
    while (sam_itr_next(bam_in, itrt, readinbam) >= 0) {
        int chrid = readinbam->core.tid;
        int stand = (int)readinbam->core.flag;
        // int qual = (int)readinbam->core.qual;//read MQ
        string cigar = getCigar(readinbam);
        uint32_t start_readalnref_pos = readinbam->core.pos + 1;
        if (chrid != -1 && (stand == 256 || stand == 272)) {
            uint32_t current_refpos = start_readalnref_pos;
            int temp_size, change_size;
            std::string currentToken;
            std::string::iterator iter = cigar.begin();
            while (iter != cigar.end()){
                char c = *iter;
                if (isdigit(c)) {
                    currentToken += c;
                } else if (!currentToken.empty()) {
                    currentToken += c;
                    string temp = currentToken;
                    temp_size = temp.size();
                    change_size = stoi(temp.substr(0,temp_size-1));
                    if (temp[temp_size - 1] == 'M' || temp[temp_size - 1] == '=' || temp[temp_size - 1] == 'X' || temp[temp_size - 1] == 'D') {
                        for (int i = 0; i < change_size; i++) {
                            secondaryReads[current_refpos] = secondaryReads[current_refpos] + 1;
                            current_refpos = current_refpos + 1;
                        }
                    }
                    currentToken.clear();
                }
                ++iter;
            }
        }
        else if (chrid != -1 && (stand == 0 || stand == 16)) {
            std :: string readname = bam_get_qname(readinbam);
            std :: string readseq = getSeq(readinbam);
            std :: string basequal = getQual(readinbam);
            uint32_t current_refpos = start_readalnref_pos;
            uint32_t current_readpos = 0;
            //calculate m/i/d
            int match = 0;
            int mismatch = 0;
            int indel = 0;// length
            int indel_num = 0;// number
            int temp_size , change_size;
            std :: vector<string> cigar_sp;
            std::string currentToken;
            std::string::iterator iter = cigar.begin();
            while (iter != cigar.end()){
                char c = *iter;
                if (isdigit(c)) {
                    currentToken += c;
                } else if (!currentToken.empty()) {
                    currentToken += c;
                    string temp = currentToken;
                    temp_size = temp.size();
                    change_size = stoi(temp.substr(0,temp_size-1));
                    if (temp[temp_size-1] == 'S'){
                        current_readpos = current_readpos + change_size;
                    }else if (temp[temp_size - 1] == 'M'){
                        std :: string refpart;
                        ref.load_ref_from_buff(give_chrid, current_refpos, change_size, refpart);
                        std :: string readpart = readseq.substr(current_readpos, change_size);
                        char Aref,Aread;
                        for (int i = 0; i < change_size; i++) {
                                Aref = refpart[i];
                                Aread = std::toupper(readpart[i]);
                                // uint32_t xpos = current_readpos+i;
                                // Aread = std::toupper(readseq[xpos]);
                                if(Aref == Aread){ // match
                                    match = match + 1;
                                } else{ // mismatch
                                    mismatch = mismatch + 1;
                                }
                            current_refpos = current_refpos + 1;
                            current_readpos = current_readpos + 1;
                        }
                    } else if (temp[temp_size - 1] == '='){
                        match = match + change_size;
                        current_refpos = current_refpos + change_size;
                        current_readpos = current_readpos + change_size;
                    } else if (temp[temp_size - 1] == 'X'){
                        mismatch = mismatch + change_size;
                        current_refpos = current_refpos + change_size;
                        current_readpos = current_readpos + change_size;
                    } else if (temp[temp_size - 1] == 'I'){
                        indel = indel + change_size;
                        indel_num = indel_num + 1;
                        current_readpos = current_readpos + change_size;
                    } else if (temp[temp_size - 1] == 'D'){
                        indel = indel + change_size;
                        indel_num = indel_num + 1;
                        current_refpos = current_refpos + change_size;
                    }
                    cigar_sp.push_back(currentToken);
                    currentToken.clear();
                }
                ++iter;
            }
            if((float)match/(match+mismatch) < opt->reads_R) {
                continue;
            }

            // process poses before the start_readalnref_pos
            std :: map<uint32_t, PosInfo>::iterator itpos = pos_infos.begin();
            while (itpos != pos_infos.end()) { // the read's chrID can find in map
                if (itpos->first > start_readalnref_pos) break;
                if(itpos->first<readbam_start) {
                    itpos = pos_infos.erase(itpos);
                    continue;
                }if(itpos->first>readbam_end) break;
                if (itpos->second.Vrts.size() <= 0) {
                    itpos = pos_infos.erase(itpos);
                    continue;
                }
                int ALLreads_num = itpos->second.ALLreads.size();
                if (ALLreads_num < 2) {
                    itpos = pos_infos.erase(itpos);
                    continue;
                }
                // candidate vrt add
                int inhomo = find_homo(dup_ran, itpos->first);//1:in 0:out
                if(inhomo == 1){
                    int secondnum = 0;
                    std ::map<uint32_t, int>::iterator itsec = secondaryReads.find(itpos->first);
                    if (itsec != secondaryReads.end()) secondnum = itsec->second;
                    int allreads_num = itpos->second.ALLreads.size();
                    if ((float)secondnum/(allreads_num+secondnum) > 0.4){
                        PosInfo addcandidatepos;
                        for (auto onevrt : itpos->second.Vrts) {
                            if ((float)onevrt.suppvrt_read / ALLreads_num < 0.2) continue;
                            int n_alts = onevrt.suppvrt_read, n_total = onevrt.suppvrt_read + itpos->second.ref_num;
                            float err = ERR_snp;
                            std ::vector<float> GL_P;
                            rescale_read_counts(n_alts, n_total);
                            GL_P = cal_GL(n_alts, n_total, err);
                            FormatInfo outinfo = computer_quals(GL_P);
                            if(outinfo.GT == "0/0" || (outinfo.GQ < 80 && outinfo.GT == "0/1")) continue;
                            addcandidatepos.Vrts.push_back(onevrt);
                            addcandidatepos.vrt_num = addcandidatepos.vrt_num + 1;
                        }
                        if (addcandidatepos.vrt_num > 0) {
                            addcandidatepos.ref_num = itpos->second.ref_num;
                            addcandidatepos.ALLreads.assign(itpos->second.ALLreads.begin(), itpos->second.ALLreads.end());
                            candidata_p1[itpos->first] = addcandidatepos;
                        }
                        itpos = pos_infos.erase(itpos);
                        continue;
                    }
                }
                PosInfo addcandidatepos;
                for (auto onevrt : itpos->second.Vrts) {
                    if ((float)onevrt.suppvrt_read / ALLreads_num < 0.2) continue;
                    addcandidatepos.Vrts.push_back(onevrt);
                    addcandidatepos.vrt_num = addcandidatepos.vrt_num + 1;
                } 
                if (addcandidatepos.vrt_num > 0) {
                    addcandidatepos.ref_num = itpos->second.ref_num;
                    addcandidatepos.ALLreads.assign(itpos->second.ALLreads.begin(), itpos->second.ALLreads.end());
                    candidata_p1[itpos->first] = addcandidatepos;
                }
                itpos = pos_infos.erase(itpos);
            }
            // add newread_pos code
            addreadposin_code(ref, cigar_sp, chrid, readname, stand, start_readalnref_pos, readseq, basequal, pos_infos);
        }
    }

    // std::clock_t c_end1 = std::clock();
    // double time_elapsed_ms1 = 1000.0 * (c_end1-c_start1) / CLOCKS_PER_SEC;
    // fprintf(stderr, "\ncreat_candidate_code P1 CPU : %.5f sec\n", time_elapsed_ms1 / 1000.0);

    // process the last part
    std ::map<uint32_t, PosInfo>::iterator itpos = pos_infos.begin();
    while (itpos != pos_infos.end()) { // the read's chrID can find in map
        if(itpos->first<readbam_start){
            itpos = pos_infos.erase(itpos);
            continue;
        }if(itpos->first>readbam_end) break;
        if (itpos->second.Vrts.size() <= 0) {
            itpos++;
            continue;
        }
        int ALLreads_num = itpos->second.ALLreads.size();
        if (ALLreads_num < 2) {
            itpos++;
            continue;
        }
        int inhomo = find_homo(dup_ran, itpos->first);//1:in 0:out
        if(inhomo == 1){
            int secondnum = 0;
            std ::map<uint32_t, int>::iterator itsec = secondaryReads.find(itpos->first);
            if (itsec != secondaryReads.end()) secondnum = itsec->second;
            int allreads_num = itpos->second.ALLreads.size();
            if ((float)secondnum/(allreads_num+secondnum) > 0.4){
                PosInfo addcandidatepos;
                for (auto onevrt : itpos->second.Vrts) {
                    if ((float)onevrt.suppvrt_read / ALLreads_num < 0.2) continue;
                    int n_alts = onevrt.suppvrt_read, n_total = onevrt.suppvrt_read + itpos->second.ref_num;
                    float err = ERR_snp;
                    std ::vector<float> GL_P;
                    rescale_read_counts(n_alts, n_total);
                    GL_P = cal_GL(n_alts, n_total, err);
                    FormatInfo outinfo = computer_quals(GL_P);
                    if(outinfo.GT == "0/0" || (outinfo.GQ < 80 && outinfo.GT == "0/1")) continue;
                    addcandidatepos.Vrts.push_back(onevrt);
                    addcandidatepos.vrt_num = addcandidatepos.vrt_num + 1;
                }
                if (addcandidatepos.vrt_num > 0) {
                    addcandidatepos.ref_num = itpos->second.ref_num;
                    addcandidatepos.ALLreads.assign(itpos->second.ALLreads.begin(), itpos->second.ALLreads.end());
                    candidata_p1[itpos->first] = addcandidatepos;
                }
                itpos++;
                continue;
            }
        }
        // candidate vrt add
        PosInfo addcandidatepos;
        for (auto onevrt : itpos->second.Vrts) {
            if ((float)onevrt.suppvrt_read / ALLreads_num < 0.2) continue;
            addcandidatepos.Vrts.push_back(onevrt);
            addcandidatepos.vrt_num = addcandidatepos.vrt_num + 1;
        }
        if (addcandidatepos.vrt_num > 0) {
            addcandidatepos.ref_num = itpos->second.ref_num;
            addcandidatepos.ALLreads.assign(itpos->second.ALLreads.begin(), itpos->second.ALLreads.end());
            candidata_p1[itpos->first] = addcandidatepos;
        }
        itpos++;
    }
    bam_destroy1(readinbam);
    hts_idx_destroy(idxt);
    bam_itr_destroy(itrt);
    bam_hdr_destroy(bam_header);
    sam_close(bam_in);
    free(region);
    pos_infos.clear();

    // aln pos with haplotype    
    // std :: clock_t c_start2 = std :: clock();

    MM_idx_loader *idx = (MM_idx_loader *)new (MM_idx_loader);
    idx->load_window_ID_idx(INDEX_PATH);
    idx->load_all_index(INDEX_PATH, NULL, true);
    hap_string_loader_single_thread hl_r1;
    int window_ID;
    std::string ref_str;
    String_list_and_var_list hl_str_v;
    uint32_t aln_start_refpos = 0;
    int start_back = 10;
    std ::map<uint32_t, PosInfo>::iterator itcandidateP1 = candidata_p1.begin(); // real_refpos , vrtinfos
    while (itcandidateP1 != candidata_p1.end()) { // which pos
        int maxvrtlength = 1;
        uint32_t realrefpos = itcandidateP1->first;
        // change variable
        if (realrefpos + maxvrtlength >= aln_start_refpos + 149) {
            aln_start_refpos = realrefpos - start_back;
            // start(the first block)
            if (realrefpos + maxvrtlength < start_back)
                aln_start_refpos = 0;
            // load from buff
            window_ID = hl_r1.get_windows_ID(give_chrint, aln_start_refpos, idx);
            // reference 300bp part
            ref.load_ref_from_buff(idx->wb_info[window_ID].chrID, idx->wb_info[window_ID].region_st, idx->wb_info[window_ID].region_length, ref_str);
            // heplotype sv
            hl_str_v = hl_r1.get_string_list_and_var_list(window_ID, ref_str, idx);
        }
        // aln pos with haplotype
        if (realrefpos + maxvrtlength < aln_start_refpos + 149) {
            uint32_t gap_s_refpos = realrefpos - idx->wb_info[window_ID].region_st;
            for (int onevrtF = 0; onevrtF < itcandidateP1->second.Vrts.size(); onevrtF++) { // one vrt
                int findinhpflag = 0;
                std ::string readsv = itcandidateP1->second.Vrts[onevrtF].ALT;
                for (auto hp_one : hl_str_v.var_l) { // one hp
                    for (int i = 0; i < hp_one.size(); i++) { // one hp's one sv
                        if (hp_one[i].ref_pos < gap_s_refpos)
                            continue;
                        if (gap_s_refpos < hp_one[i].ref_pos)
                            break;

                        if (gap_s_refpos == hp_one[i].ref_pos && itcandidateP1->second.Vrts[onevrtF].length_in_ref == hp_one[i].REF.size()) {
                            std ::string hp_vrtstr = hp_one[i].ALT;
                            if (hp_vrtstr == readsv) {
                                itcandidateP1->second.Vrts[onevrtF].isknown = 0;
                                itcandidateP1->second.Vrts[onevrtF].REF = hp_one[i].REF;
                                findinhpflag = 1;
                                break;
                            }
                        }
                    }
                    if (findinhpflag == 1) break;
                }
                if (findinhpflag == 0) {
                    itcandidateP1->second.Vrts[onevrtF].isknown = 1;
                    itcandidateP1->second.Vrts[onevrtF].REF = ref_str.substr(gap_s_refpos, itcandidateP1->second.Vrts[onevrtF].length_in_ref);
                }
            }
        }
        itcandidateP1++;
    }
    delete (idx);
    idx = nullptr;
    // std::clock_t c_end2 = std::clock();
    // double time_elapsed_ms2 = 1000.0 * (c_end2-c_start2) / CLOCKS_PER_SEC;
    // fprintf(stderr, "\ncreat_candidate_code P2 CPU : %.5f sec\n", time_elapsed_ms2 / 1000.0);
}

void create_highvrt_p2vrt(char *BED_PATH, std :: map<uint32_t , PosInfo> candidata_p1, std :: vector <snvVRT> &ResultSnv, std :: map <string, vector<vrtPOS_bam>> &p2vrtpos, std :: map<uint32_t, std :: vector <P2vrt>> &p2vrt){
    std ::vector<uint32_t> snv_posline;
    snv_posline.push_back(0);
    std ::map<uint32_t, PosInfo>::iterator itcan1 = candidata_p1.begin();
    while (itcan1 != candidata_p1.end()) {
        snv_posline.push_back(itcan1->first);
        itcan1++;
    }
    snv_posline.push_back(4294967294);
    sort(snv_posline.begin(), snv_posline.end());

    std :: vector <range_homo> homo_ran;
    load_bed(BED_PATH, opt->chrID, homo_ran);

    uint32_t per_pos, next_pos;
    for (int snvpos_inx = 1; snvpos_inx < snv_posline.size() - 1; snvpos_inx++) { // every pos process
        uint32_t realrefpos_now = snv_posline[snvpos_inx];
        std ::map<uint32_t, PosInfo>::iterator itcan2 = candidata_p1.find(realrefpos_now);
        if (itcan2 != candidata_p1.end()) { // find all vrt in the pos
            per_pos = snv_posline[snvpos_inx - 1];
            next_pos = snv_posline[snvpos_inx + 1];
            int idx = (realrefpos_now - 25) / 50;
            uint32_t idx_refpos = idx * 50;
            int vrtNUM = itcan2->second.vrt_num;
            int supprefnum = itcan2->second.ref_num;
            int ALLreads_num = itcan2->second.ALLreads.size();

            int if2 = 0;
            int inhomo = find_homo(homo_ran, realrefpos_now);//1:in 0:out
            if (ALLreads_num >= 3) {
                for (auto onevrt : itcan2->second.Vrts) {
                    int line0 = onevrt.suppstand0, line16 = onevrt.suppstand16;
                    int del = string_bias(line0, line16);
                    if(del == 1) continue;

                    int ifResult = 0;
                    int VaRnum = onevrt.suppvrt_read + supprefnum;
                    float rate_ = (float)onevrt.suppvrt_read/ALLreads_num;

                    int n_alts = onevrt.suppvrt_read, n_total = VaRnum;
                    float err = ERR_snp;
                    std ::vector<float> GL_P;
                    rescale_read_counts(n_alts, n_total);
                    GL_P = cal_GL(n_alts, n_total, err);
                    FormatInfo outinfo = computer_quals(GL_P);
                    if(outinfo.GT == "0/0") continue;

                    if(inhomo == 1){//in homo
                    }else{
                        if((onevrt.isknown == 0 && rate_ >= 0.3) || (rate_ >= 0.7 && (per_pos <= itcan2->first-10 && next_pos>=itcan2->first+10 && vrtNUM == 1))){
                            ifResult = 1;
                        }
                    }
                     
                    if (ifResult == 1) { // part1 result
                        snvVRT c1;
                        c1.ref_pos = itcan2->first;
                        c1.vrtAref_num = VaRnum;
                        c1.vrt_num = onevrt.suppvrt_read;
                        c1.all_num = ALLreads_num;
                        c1.outinfo = outinfo;
                        c1.REF = onevrt.REF;
                        c1.ALT = onevrt.ALT;
                        ResultSnv.push_back(c1);
                    } else{ // part2 vrt
                        P2vrt snv1;
                        snv1.outinfo = outinfo;
                        snv1.refpos = itcan2->first;
                        snv1.REF = onevrt.REF;
                        snv1.ALT = onevrt.ALT;
                        snv1.vrt_num = onevrt.suppvrt_read;
                        snv1.all_num = ALLreads_num;
                        snv1.inhomo = inhomo;
                        if((per_pos > itcan2->first-10 && next_pos < itcan2->first+10) || vrtNUM > 1) snv1.vrt_num_slop10 = 1;
                        p2vrt[idx_refpos].push_back(snv1);
                        if2 = 1;
                    }
                }
            } else {
                for (auto onevrt : itcan2->second.Vrts) {
                    if(onevrt.suppvrt_read == 1) continue;
                    int ifResult = 0;
                    int VaRnum = onevrt.suppvrt_read + supprefnum;
                    float rate_ = (float)onevrt.suppvrt_read/ALLreads_num;
                    
                    int n_alts = onevrt.suppvrt_read, n_total = VaRnum;
                    float err = ERR_snp;
                    std ::vector<float> GL_P;
                    rescale_read_counts(n_alts, n_total);
                    GL_P = cal_GL(n_alts, n_total, err);
                    FormatInfo outinfo = computer_quals(GL_P);
                    if(outinfo.GT == "0/0") continue;

                    if(inhomo == 1){//in homo
                    }else{
                        if((onevrt.isknown == 0 && rate_ >= 0.3) || (rate_ >= 0.7 && (per_pos <= itcan2->first-10 && next_pos>=itcan2->first+10 && vrtNUM == 1))){
                            ifResult = 1;
                        }
                    }
                    if (ifResult == 1) {
                        snvVRT c1;
                        c1.ref_pos = itcan2->first;
                        c1.vrtAref_num = VaRnum;
                        c1.vrt_num = onevrt.suppvrt_read;
                        c1.all_num = ALLreads_num;
                        c1.outinfo = outinfo;
                        c1.REF = onevrt.REF;
                        c1.ALT = onevrt.ALT;
                        ResultSnv.push_back(c1);
                    } else{ // part2 vrt
                        P2vrt snv1;
                        snv1.outinfo = outinfo;
                        snv1.refpos = itcan2->first;
                        snv1.REF = onevrt.REF;
                        snv1.ALT = onevrt.ALT;
                        snv1.vrt_num = onevrt.suppvrt_read;
                        snv1.all_num = ALLreads_num;
                        snv1.inhomo = inhomo;
                        if((per_pos > itcan2->first-10 && next_pos < itcan2->first+10) || vrtNUM > 1) snv1.vrt_num_slop10 = 1;
                        p2vrt[idx_refpos].push_back(snv1);
                        if2 = 1;
                    }
                }
            }
            // part2 reads pos creat
            if (if2 == 1) {
                for (auto read : itcan2->second.ALLreads) {
                    std ::map<string, vector<vrtPOS_bam>>::iterator it = p2vrtpos.find(read.readid);
                    if (it != p2vrtpos.end()) { // can find readid
                        vrtPOS_bam tmpsnvpos;
                        tmpsnvpos.stand = read.stand;
                        std ::vector<vrtPOS_bam>::iterator itstand = find(it->second.begin(), it->second.end(), tmpsnvpos);
                        if (itstand != it->second.end()) {
                            std ::vector<uint32_t>::iterator itsvpos = find(itstand->ref_pos_idx.begin(), itstand->ref_pos_idx.end(), idx_refpos);
                            if (itsvpos != itstand->ref_pos_idx.end()){}
                            else itstand->ref_pos_idx.push_back(idx_refpos);
                        } else {
                            tmpsnvpos.ref_pos_idx.push_back(idx_refpos);
                            it->second.push_back(tmpsnvpos);
                        }
                    } else { // cannot find readid
                        vrtPOS_bam tmpsnvpos;
                        tmpsnvpos.stand = read.stand;
                        tmpsnvpos.ref_pos_idx.push_back(idx_refpos);
                        p2vrtpos[read.readid].push_back(tmpsnvpos);
                    }
                }
            }
        }
    }
}

//part 2--------------------------------------------------------------------------------------


void load_read_part2_bam(char *BAM_PATH, std :: map <string, std :: vector<vrtPOS_bam>>  p2vrtpos, std :: map <int,std :: vector <readstr_name>> &clupos_str, int give_chrid, uint32_t readbam_start, uint32_t readbam_end){
    char * region = (char *)malloc(1024);
    snprintf(region, 1024, "%s:%d-%d", opt->chrID, opt->readbam_start, opt->readbam_end);

    samFile *bam_file = sam_open(BAM_PATH, "r");
    samFile *bam_in = sam_open(BAM_PATH, "r"); // open bam file
    bam_hdr_t *bam_header = sam_hdr_read(bam_file);
    hts_idx_t *idxt = hts_idx_load(BAM_PATH, 1);
    hts_itr_t *itrt = bam_itr_querys(idxt, bam_header, region);
    bam1_t *readinbam = bam_init1();
    while (sam_itr_next(bam_in, itrt, readinbam) >= 0) {
        string queryname = bam_get_qname(readinbam);
        string seq = getSeq(readinbam);
        uint32_t start_in_ref = readinbam->core.pos + 1;
        std ::string cigar = getCigar(readinbam);
        std ::map<string, vector<vrtPOS_bam>>::iterator itvrtpos;
        itvrtpos = p2vrtpos.find(queryname);
        if (itvrtpos != p2vrtpos.end()) {
            vrtPOS_bam tmpp;
            tmpp.stand = (int)readinbam->core.flag;
            std ::vector<vrtPOS_bam>::iterator itsameP1 = find(itvrtpos->second.begin(), itvrtpos->second.end(), tmpp);
            if (itsameP1 != itvrtpos->second.end()) {
                for (auto svpos : itsameP1->ref_pos_idx) { // one id_pos
                    if(svpos + 100 <= start_in_ref) continue;
                    int ref_diff_start = svpos - start_in_ref;
                    int ref_diff_end = ref_diff_start + 100;
                    int read_diff = 0, ref_tmpdiff = 0;
                    int s_read_pos = -1, e_read_pos = -1;
                    std::string tmp_curr;
                    int temp_size, change_size;

                    std::string currentToken;
                    std::string::iterator iter = cigar.begin();
                    while (iter != cigar.end()){
                        char c = *iter;
                        if (isdigit(c)) {
                            currentToken += c;
                        } else if (!currentToken.empty()) {
                            currentToken += c;
                            string temp = currentToken;
                            if (ref_diff_start > 0) {
                                tmp_curr = temp;
                                temp_size = temp.size();
                                change_size = stoi(temp.substr(0, temp_size - 1));
                                if (temp[temp_size - 1] == 'S') {
                                    read_diff = read_diff + change_size;
                                    if (ref_tmpdiff >= ref_diff_start && s_read_pos == -1) s_read_pos = read_diff;
                                    if (ref_tmpdiff >= ref_diff_end) break;
                                } else if (temp[temp_size - 1] == 'I') {
                                    read_diff = read_diff + change_size;
                                    if (ref_tmpdiff >= ref_diff_start && s_read_pos == -1) s_read_pos = read_diff;
                                    if (ref_tmpdiff >= ref_diff_end) break;
                                } else if (temp[temp_size - 1] == 'D') {
                                    ref_tmpdiff = ref_tmpdiff + change_size;
                                    if (ref_tmpdiff >= ref_diff_start && s_read_pos == -1) s_read_pos = read_diff;
                                    if (ref_tmpdiff >= ref_diff_end) break;
                                } else if (temp[temp_size - 1] == 'M' || temp[temp_size - 1] == '=' || temp[temp_size - 1] == 'X') {
                                    ref_tmpdiff = ref_tmpdiff + change_size;
                                    if (ref_tmpdiff >= ref_diff_start && s_read_pos == -1) {
                                        int change_size_p = change_size - (ref_tmpdiff - ref_diff_start);
                                        s_read_pos = read_diff + change_size_p;
                                    } if (ref_tmpdiff > ref_diff_end) {
                                        change_size = change_size - (ref_tmpdiff - ref_diff_end);
                                        read_diff = read_diff + change_size;
                                        break;
                                    }
                                    read_diff = read_diff + change_size;
                                }
                            } else{
                                tmp_curr = temp;
                                temp_size = temp.size();
                                change_size = stoi(temp.substr(0, temp_size - 1));
                                if (temp[temp_size - 1] == 'S') {
                                    read_diff = read_diff + change_size;
                                    if (s_read_pos == -1) s_read_pos = read_diff;
                                    if (ref_tmpdiff >= ref_diff_end) break;
                                } else if (temp[temp_size - 1] == 'I') {
                                    read_diff = read_diff + change_size;
                                    if (s_read_pos == -1) s_read_pos = read_diff;
                                    if (ref_tmpdiff >= ref_diff_end) break;
                                } else if (temp[temp_size - 1] == 'D') {
                                    ref_tmpdiff = ref_tmpdiff + change_size;
                                    if (s_read_pos == -1) s_read_pos = read_diff;
                                    if (ref_tmpdiff >= ref_diff_end) break;
                                } else if (temp[temp_size - 1] == 'M' || temp[temp_size - 1] == '=' || temp[temp_size - 1] == 'X') {
                                    ref_tmpdiff = ref_tmpdiff + change_size;
                                    if (s_read_pos == -1) s_read_pos = read_diff;
                                    if (ref_tmpdiff > ref_diff_end) {
                                        change_size = change_size - (ref_tmpdiff - ref_diff_end);
                                        read_diff = read_diff + change_size;
                                        break;
                                    }
                                    read_diff = read_diff + change_size;
                                }
                            }
                            currentToken.clear();
                        }
                        ++iter;
                    }
                    if(ref_tmpdiff <= ref_diff_start) continue;

                    e_read_pos = read_diff;
                    int sublen = e_read_pos - s_read_pos;
                    int linelen = seq.size();
                    if (tmp_curr[temp_size - 1] == 'S') linelen = linelen - change_size;
                    if (e_read_pos > linelen) sublen = linelen - s_read_pos;
                    if (sublen < 10) continue;
                    int clu_index = svpos / 50;
                    std ::string str1 = seq.substr(s_read_pos, sublen);
                    readstr_name r1;
                    r1.str = str1;
                    r1.name = queryname;
                    clupos_str[clu_index].push_back(r1);
                }
            }
        }
    }
    bam_destroy1(readinbam);
    hts_idx_destroy(idxt);
    bam_itr_destroy(itrt);
    bam_hdr_destroy(bam_header);
    sam_close(bam_in);
    free(region);
}

void mantaAss_m(std :: map <int,std :: vector <readstr_name>>  clupos_str, std :: vector <ConsensusINFO> &consensus_str) {
	std :: map <int,std :: vector <readstr_name>> :: iterator itpos = clupos_str.begin();
	while(itpos != clupos_str.end()){//which pos
		//length of the list
		int seq_n = itpos->second.size();
		AssemblyManager *am_ = new AssemblyManager[1];
		AssemblyManager &am = *am_;
		am.o.minWordLength = 30;
		am.o.maxWordLength = 120;
		for(int i = 0; i < seq_n; ++i) am.reads.emplace_back(itpos->second[i].str);
		am.assembley();
		auto &contigs = am.getResults();
		int cons_n = contigs.size();
		if(cons_n > 0){
			ConsensusINFO css_;
			// css_.chrID = itpos->second[0].chrID;
			css_.all_cov_readnum = seq_n;
			css_.ref_pos = itpos->first * 50;
			for(auto &contig:contigs) {
				CON_suppR s1;
				s1.consensus = contig.seq;
				s1.suppsnv_readnum = contig.supportReads.size();
				css_.con_suppr.push_back(s1);
			}
			consensus_str.push_back(css_);
		}
		delete []am_;
		itpos++;
	}
}

void mantaAss_mSTEP3(std :: map<uint32_t, std :: vector <P2vrt>> p2vrt, std :: map <int,std :: vector <readstr_name>>  clupos_str, std :: map<string,int> read_phased, std :: vector <ConsensusINFO> &consensus_str, std :: vector <lowsnvVRT> &ResultSnv) {
	std :: map <int,std :: vector <readstr_name>> :: iterator itpos_idx = clupos_str.begin();
	while(itpos_idx != clupos_str.end()){//which pos
        int HP1 = 0, HP2 = 0;
        std :: vector <readstr_name> HP1info,HP2info;
        for(auto oneread : itpos_idx->second){
            std :: map<string,int> :: iterator itreadP = read_phased.find(oneread.name);
            if(itreadP != read_phased.end()){
                if(itreadP->second == 1) {
                    HP1info.push_back(oneread);
                    HP1 = HP1 +1;
                }
                if(itreadP->second == 2) {
                    HP2info.push_back(oneread);
                    HP2 = HP2 +1;
                }
            }
        }
        uint32_t idxpos = itpos_idx->first * 50;
        int allnum = p2vrt[idxpos][0].all_num;
        ConsensusINFO css_;
        css_.all_cov_readnum = itpos_idx->second.size();
        css_.ref_pos = idxpos;
        if(HP1+HP2 < allnum*0.1) {
            for (auto onep2vrt: p2vrt[idxpos]) {
                lowsnvVRT addone;
                addone.ref_pos = onep2vrt.refpos;
                addone.REF = onep2vrt.REF;
                addone.ALT = onep2vrt.ALT;
                addone.outinfo = onep2vrt.outinfo;
                ResultSnv.push_back(addone);
            }
            itpos_idx++;
            continue;
        }

        if(HP1 != 0){
            //length of the list
            int seq_n = HP1info.size();
            AssemblyManager *am_ = new AssemblyManager[1];
            AssemblyManager &am = *am_;
            am.o.minWordLength = 30;
            am.o.maxWordLength = 120;
            for(int i = 0; i < seq_n; ++i) am.reads.emplace_back(HP1info[i].str);
            am.assembley();
            auto &contigs = am.getResults();
            int cons_n = contigs.size();
            if(cons_n > 0) {
                for(auto &contig:contigs) {
                    CON_suppR s1;
                    s1.consensus = contig.seq;
                    s1.suppsnv_readnum = HP1;
                    css_.con_suppr.push_back(s1);
                }
            }
            delete []am_;
        }
        if(HP2 != 0){
            //length of the list
            int seq_n = HP2info.size();
            AssemblyManager *am_ = new AssemblyManager[1];
            AssemblyManager &am = *am_;
            am.o.minWordLength = 30;
            am.o.maxWordLength = 120;
            for(int i = 0; i < seq_n; ++i) am.reads.emplace_back(HP2info[i].str);
            am.assembley();
            auto &contigs = am.getResults();
            int cons_n = contigs.size();
            if(cons_n > 0) {
                for(auto &contig:contigs) {
                    CON_suppR s1;
                    s1.consensus = contig.seq;
                    s1.suppsnv_readnum = HP2;
                    css_.con_suppr.push_back(s1);
                }
            }
            delete []am_;
        }
        consensus_str.push_back(css_);
		itpos_idx++;
	}
}

void abpoa_consensus(int chrID, std :: map <int,std :: vector <readstr_name>> clupos_str, std ::vector<ConsensusINFO> &consensus_str) {
    abpoa_para_t *abpt = abpoa_init_para();
    // output options
    abpt->out_msa = 0;  // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    abpt->w = 6, abpt->k = 9;
    abpt->min_w = 10; // minimizer-based seeding and partition
    // abpt->min_freq
    abpt->progressive_poa = 1;
    abpt->max_n_cons = 2; // to generate 2 consensus sequences
    abpoa_post_set_para(abpt);

    std :: map <int,std :: vector <readstr_name>>::iterator itpos = clupos_str.begin();
    while (itpos != clupos_str.end()) { // which pos
        // initialize variables
        abpoa_t *ab = abpoa_init();
        int n_seqs = itpos->second.size();
        int *seq_lens = (int *)malloc(sizeof(int) * n_seqs);
        uint8_t **bseqs = (uint8_t **)malloc(sizeof(uint8_t *) * n_seqs);
        int **weights = (int **)malloc(sizeof(int *) * n_seqs);
        int i = 0;
        for (auto seq : itpos->second) {
            seq_lens[i] = seq.str.length();
            bseqs[i] = (uint8_t *)malloc(sizeof(uint8_t) * seq_lens[i]);
            weights[i] = (int *)malloc(sizeof(int) * seq_lens[i]);
            for (int j = 0; j < seq_lens[i]; ++j) {
                bseqs[i][j] = nt4_table[(int)seq.str[j]];
                if (j >= 12) weights[i][j] = 2;
                else weights[i][j] = 0;
            }
            i = i + 1;
        }
        abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, weights);
        abpoa_cons_t *abc = ab->abc;

        ConsensusINFO css_;
        // css_.chrID = chrID;
        css_.all_cov_readnum = n_seqs;
        css_.ref_pos = itpos->first * 50;
        for (i = 0; i < abc->n_cons; ++i) {
            std ::string css_str = "";
            for (int j = 0; j < abc->cons_len[i]; ++j) css_str.push_back(nt256_table[abc->cons_base[i][j]]);
            string oneconsensus = css_str.substr(0, css_str.length());
            CON_suppR s1;
            s1.consensus = oneconsensus;
            s1.suppsnv_readnum = abc->clu_n_seq[i];
            css_.con_suppr.push_back(s1);
        }
        consensus_str.push_back(css_);
        // free seq-related variables
        for (i = 0; i < n_seqs; ++i) {
            free(bseqs[i]);
            free(weights[i]);
        }
        free(bseqs);
        free(seq_lens);
        free(weights);
        // free abpoa-related variables
        abpoa_free(ab);
        itpos++;
    }
    abpoa_free_para(abpt);
}

void ksw_score(std :: vector <readstr_name> reads_str, std :: string  refstr, std :: vector <CON_suppR> &consensus_one) {
    int suppnum1 = 0, suppnum2 = 0;
    int score1 = 0, score2 = 0;
    string cigar = "";
    int con_num = consensus_one.size();
    if (con_num == 1) {
        for (auto read_one : reads_str) {
            score1 = ksw_align(consensus_one[0].consensus, read_one.str, consensus_one[0].consensus.length(), read_one.str.length(), 1, 4, 8, 2, cigar);
            score2 = ksw_align(refstr, read_one.str, refstr.length(), read_one.str.length(), 1, 4, 8, 2, cigar);
            if (score1 > 0 && score1 > score2)
                suppnum1 = suppnum1 + 1;
        }
        consensus_one[0].suppsnv_readnum = suppnum1;
    } else if (con_num == 2) {
        for (auto read_one : reads_str) {
            score1 = ksw_align(consensus_one[0].consensus, read_one.str, consensus_one[0].consensus.length(), read_one.str.length(), 1, 4, 8, 2, cigar);
            score2 = ksw_align(consensus_one[1].consensus, read_one.str, consensus_one[1].consensus.length(), read_one.str.length(), 1, 4, 8, 2, cigar);
            if (score1 > 0 && score1 > score2)
                suppnum1 = suppnum1 + 1;
            if (score2 > 0 && score1 < score2)
                suppnum2 = suppnum2 + 1;
        }
        consensus_one[0].suppsnv_readnum = suppnum1;
        consensus_one[1].suppsnv_readnum = suppnum2;
    }
}

int ksw_align(std :: string refstr, std :: string  readstr, int len_tl, int len_ql, int sc_mch, int sc_mis, int gapo, int gape, std :: string &cigar) {
    int i, a = sc_mch, b = sc_mis < 0 ? sc_mis : -sc_mis; // a>0 and b<0
    int8_t mat[25] = {a, b, b, b, 0, b, a, b, b, 0, b, b, a, b, 0, b, b, b, a, 0, 0, 0, 0, 0, 0};
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0;
    c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2;
    c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t *)malloc(len_tl);
    qs = (uint8_t *)malloc(len_ql);
    for (i = 0; i < len_tl; ++i) ts[i] = c[(uint8_t)refstr[i]]; // encode to 0/1/2/3
    for (i = 0; i < len_ql; ++i) qs[i] = c[(uint8_t)readstr[i]];
    ksw_extz2_sse(0, len_ql, qs, len_tl, ts, 5, mat, gapo, gape, -1, -1, -1, 0, &ez);
    for (i = 0; i < ez.n_cigar; ++i) { // creat CIGAR
        int a = ez.cigar[i] >> 4;
        char b = "MID"[ez.cigar[i] & 0xf];
        cigar += to_string(a);
        cigar.push_back(b);
    }
    free(ez.cigar);
    free(ts);
    free(qs);
    return ez.score;
}

int edlib_used(std :: string  qstr, int len_ql, std :: string  tstr, int len_tl) {
    char *qC = (char *)qstr.c_str();
    char *tC = (char *)tstr.c_str();
    int distance = -1;
    EdlibAlignResult result = edlibAlign(qC, len_ql, tC, len_tl, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) distance = result.editDistance;
    edlibFreeAlignResult(result);
    return distance;
}

int edlib_distance(std :: vector <readstr_name> reads_str, std :: string refstr, std :: vector <CON_suppR> &consensus_one)
{
    int suppnum1 = 0, suppnum2 = 0, allnum = 0;
    int dis1 = -1, dis2 = -1;
    int con_num = consensus_one.size();
    if (con_num == 1) {
        for (auto read_one : reads_str) {
            dis1 = edlib_used(read_one.str, read_one.str.length(), consensus_one[0].consensus, consensus_one[0].consensus.length());
            dis2 = edlib_used(read_one.str, read_one.str.length(), refstr, refstr.length());
            if (dis2 > dis1 && dis1 >= 0) suppnum1 = suppnum1 + 1;
            if (dis2 < dis1 && dis2 >= 0) suppnum2 = suppnum2 + 1;
        }
        consensus_one[0].suppsnv_readnum = suppnum1;
    } else if (con_num == 2) {
        for (auto read_one : reads_str) {
            dis1 = edlib_used(read_one.str, read_one.str.length(), consensus_one[0].consensus, consensus_one[0].consensus.length());
            dis2 = edlib_used(read_one.str, read_one.str.length(), consensus_one[1].consensus, consensus_one[1].consensus.length());
            if (dis2 > dis1 && dis1 >= 0) suppnum1 = suppnum1 + 1;
            if (dis2 < dis1 && dis2 >= 0) suppnum2 = suppnum2 + 1;
        }
        consensus_one[0].suppsnv_readnum = suppnum1;
        consensus_one[1].suppsnv_readnum = suppnum2;
    }
    allnum = suppnum1 + suppnum2;
    return allnum;
}

void analyze_K2snv(std :: string refstr, std :: string consensus, std :: string cigar, uint32_t startrefpos, int suppallnum, int suppsnvnum, int chrID, std :: vector <P2vrt> posp2vrt, std :: vector <lowsnvVRT> &ReaultV_add) {

    int pos_start_ref = 0, pos_start_read = 0;
    std ::string cigarX = "";

    std :: vector<string> cigar_sp1;
    std::string currentToken;  // 保存当前正在提取的数字或字符
    std::string::iterator iter1 = cigar.begin();
    while (iter1 != cigar.end()){
        char c = *iter1;
        if (isdigit(c)) {
            currentToken += c;
        } else if (!currentToken.empty()) {
            currentToken += c;
            cigar_sp1.push_back(currentToken);
            currentToken.clear();
        }
        ++iter1;
    }
    
    for(auto temp : cigar_sp1){
        int temp_size1 = temp.size();
        int change_size1 = stoi(temp.substr(0, temp_size1 - 1));
        if (temp[temp_size1 - 1] == 'M') {
            string refs = refstr.substr(pos_start_ref, change_size1);
            string cons = consensus.substr(pos_start_read, change_size1);
            int match = 0, mismatch = 0;
            for (int i = 0; i < change_size1; i++) {
                if (refs[i] == cons[i]) {
                    if (mismatch != 0) {
                        cigarX += to_string(mismatch);
                        cigarX.push_back('X');
                        mismatch = 0;
                    }
                    match = match + 1;
                } else {
                    if (match != 0) {
                        cigarX += to_string(match);
                        cigarX.push_back('=');
                        match = 0;
                    }
                    mismatch = mismatch + 1;
                }
            }
            if (match != 0) {
                cigarX += to_string(match);
                cigarX.push_back('=');
            } else if (mismatch != 0) {
                cigarX += to_string(mismatch);
                cigarX.push_back('X');
            }
            pos_start_ref = pos_start_ref + change_size1;
            pos_start_read = pos_start_read + change_size1;
        } else if (temp[temp_size1 - 1] == 'I') {
            cigarX += temp;
            pos_start_read = pos_start_read + change_size1;
        } else if (temp[temp_size1 - 1] == 'D') {
            cigarX += temp;
            pos_start_ref = pos_start_ref + change_size1;
        }
    }

    int length = refstr.size();
    int delrangeB = length * 0.2, delrangeE = length * 0.8;
    pos_start_ref = 0, pos_start_read = 0;
    std ::string cigarp2 = "";
    int addif = 0, st_ref, st_read;

    std :: vector<string> cigar_sp2;
    std::string::iterator iter2 = cigarX.begin();
    while (iter2 != cigarX.end()){
        char c = *iter2;
        if (isdigit(c)) {
            currentToken += c;
        } else if (!currentToken.empty()) {
            currentToken += c;
            cigar_sp2.push_back(currentToken);
            currentToken.clear();
        }
        ++iter2;
    }
    
    for(auto temp : cigar_sp2){
        int temp_size = temp.size();
        int change_size = stoi(temp.substr(0, temp_size - 1));
        if (temp[temp_size - 1] == '=' || temp[temp_size - 1] == 'X') {
            pos_start_ref = pos_start_ref + change_size;
            pos_start_read = pos_start_read + change_size;
        } else if (temp[temp_size - 1] == 'I') {
            pos_start_read = pos_start_read + change_size;
        } else if (temp[temp_size - 1] == 'D') {
            pos_start_ref = pos_start_ref + change_size;
        }
        if (pos_start_ref < delrangeB) {
            continue;
        } if (pos_start_ref > delrangeE) break;
        if (addif == 1) cigarp2 += temp;
        if (pos_start_ref >= delrangeB) {
            if (addif == 0) {
                st_ref = pos_start_ref;
                st_read = pos_start_read;
            }
            addif = 1;
        }
    }

    lowsnvVRT addReaultV;
    std :: vector<string> cigar_sp3;
    std::string::iterator iter3 = cigarp2.begin();
    while (iter3 != cigarp2.end()){
        char c = *iter3;
        if (isdigit(c)) {
            currentToken += c;
        } else if (!currentToken.empty()) {
            currentToken += c;
            cigar_sp3.push_back(currentToken);
            currentToken.clear();
        }
        ++iter3;
    }
    for(auto temp : cigar_sp3){
        int temp_size = temp.size();
        int change_size = stoi(temp.substr(0, temp_size - 1));
        P2vrt findflag;

        if (temp[temp_size - 1] == '=') {
            st_ref = st_ref + change_size;
            st_read = st_read + change_size;
            continue;
        } else if (temp[temp_size - 1] == 'X') {
            for (int mis_i = 0; mis_i < change_size; mis_i++) {
                // string aref = refstr.substr(st_ref, 1);
                // string aalt = consensus.substr(st_read, 1);
                char aref = refstr[st_ref];
                char aalt = consensus[st_read];

                findflag.refpos = startrefpos + st_ref;
                findflag.REF = aref;
                findflag.ALT = aalt;
                std ::vector<P2vrt>::iterator itk2 = find(posp2vrt.begin(), posp2vrt.end(), findflag);
                if (itk2 != posp2vrt.end()) {
                    addReaultV.REF = aref;
                    addReaultV.ALT = aalt;
                    addReaultV.ref_pos = startrefpos + st_ref;
                    addReaultV.outinfo.GQ = itk2->outinfo.GQ;
                    addReaultV.outinfo.GT = itk2->outinfo.GT;
                    addReaultV.outinfo.QUAL = itk2->outinfo.QUAL;
                    addReaultV.all_num = itk2->all_num;
                    addReaultV.vrt_num = itk2->vrt_num;
                    addReaultV.inhomo = itk2->inhomo;
                    addReaultV.vrt_num_slop10 = itk2->vrt_num_slop10;
                    std :: vector <lowsnvVRT> ::iterator findinReaultV = find(ReaultV_add.begin(), ReaultV_add.end(), addReaultV);
                    if (findinReaultV != ReaultV_add.end()){
                        findinReaultV->supportnum = findinReaultV->supportnum + suppsnvnum;
                        findinReaultV->HP_num = findinReaultV->HP_num + 1;
                    } else {
                        addReaultV.supportnum = suppsnvnum;
                        addReaultV.HP_num = 1;
                        ReaultV_add.push_back(addReaultV);
                    }
                }
                st_ref = st_ref + 1;
                st_read = st_read + 1;
            }
        } else if (temp[temp_size - 1] == 'I') {
            st_read = st_read + change_size;
        } else if (temp[temp_size - 1] == 'D') {
            st_ref = st_ref + change_size;
        }
    }
}

void Re_align(std :: map<uint32_t, std :: vector <P2vrt>> p2vrt, std :: vector <ConsensusINFO> consensus_str, std :: map <int,std :: vector <readstr_name>> clupos_str, std :: vector <lowsnvVRT> &ResultSnv){
    Simple_ref_handler ref;
    ref.load_bin_ref(opt->index_path);
    MM_idx_loader *idx = (MM_idx_loader *)new (MM_idx_loader);
    idx -> load_window_ID_idx(opt->index_path);
    idx -> load_all_index(opt->index_path, NULL, true);
    for(auto constr : consensus_str){
        //load sub reference by consensus length
        hap_string_loader_single_thread hl_r1; 
        int window_ID = hl_r1.get_windows_ID(give_chrint, constr.ref_pos, idx);
        std::string ref_str;
        ref.load_ref_from_buff(idx->wb_info[window_ID].chrID, idx->wb_info[window_ID].region_st, idx->wb_info[window_ID].region_length, ref_str);
        std :: map<uint32_t, std :: vector <P2vrt>> :: iterator itp2vrt = p2vrt.find(constr.ref_pos);
        if (itp2vrt != p2vrt.end()) {
            std :: vector <string> subref_v;
            int start = constr.ref_pos - idx->wb_info[window_ID].region_st;
            for(auto oneconsensus : constr.con_suppr){
                int strlength = 100;
                string subref;
                if (ref_str.length() >= strlength + start){
                    subref = ref_str.substr(start,strlength);
                } else{
                    uint32_t refpostmp = constr.ref_pos + 25;
                    window_ID = hl_r1.get_windows_ID(give_chrint, refpostmp, idx);
                    ref.load_ref_from_buff(idx->wb_info[window_ID].chrID, idx->wb_info[window_ID].region_st, idx->wb_info[window_ID].region_length, ref_str);
                    start = constr.ref_pos + 1 - idx->wb_info[window_ID].region_st;
                    subref = ref_str.substr(start,strlength);
                }
                subref_v.push_back(subref);
            }

            int idxx = constr.ref_pos/50;
            // ksw_score(readstr_clu[give_chrint].clupos_str[idxx], constr.con_suppr, subref_v[0]);
            // int allnum = edlib_distance(clupos_str[idxx], subref_v[0], constr.con_suppr);
            int concluid = 0;
            std :: vector <lowsnvVRT> ResultV_add;
            for(auto oneconsensus : constr.con_suppr){
                int strlength = subref_v[concluid].length();
                int strlengthcon = oneconsensus.consensus.length();
                std :: string cigar = "";
                int score= ksw_align(subref_v[concluid], oneconsensus.consensus, strlength, strlengthcon, 1, 4, 8, 2, cigar);
                if(score > 0){
                    analyze_K2snv(subref_v[concluid], oneconsensus.consensus, cigar, constr.ref_pos, constr.all_cov_readnum, oneconsensus.suppsnv_readnum, give_chrint , itp2vrt->second , ResultV_add);
                }
                concluid = concluid + 1;
            }
            for(auto ResultV_one : ResultV_add) {
                if (ResultV_one.outinfo.GT == "1/1" && (ResultV_one.HP_num == 2 || ResultV_one.HP_num == 1)) {
                    ResultSnv.push_back(ResultV_one);
                    continue;
                }
                if(ResultV_one.outinfo.GT == "0/1" && ResultV_one.HP_num == 1){
                    ResultSnv.push_back(ResultV_one);
                    continue;
                }
                if(ResultV_one.inhomo == 1){
                    if (ResultV_one.outinfo.GT == "0/1" && ResultV_one.vrt_num_slop10==0 && ResultV_one.HP_num == 2) {
                        lowsnvVRT changeone;
                        changeone.ref_pos = ResultV_one.ref_pos;
                        changeone.REF = ResultV_one.REF;
                        changeone.ALT = ResultV_one.ALT;
                        changeone.outinfo.GQ = ResultV_one.outinfo.GQ;
                        changeone.outinfo.GT = "1/1";
                        changeone.outinfo.QUAL = ResultV_one.outinfo.QUAL;
                        ResultSnv.push_back(changeone);
                    }
                }
            }
        }
    }
    delete(idx);
    idx = nullptr;
}

int printS1vcf(char *STEP1vcfInfo_path, char *highvcf_path, std :: vector <snvVRT> ResultSnv) {
    if(ResultSnv.empty()) return 2;
    std :: vector <uint32_t> ResultSnv_refpos;
    sort(ResultSnv.begin() , ResultSnv.end());
    for (auto cc : ResultSnv) ResultSnv_refpos.push_back(cc.ref_pos);
    ResultSnv_refpos.push_back(4294967294);
    sort(ResultSnv_refpos.begin() , ResultSnv_refpos.end());

    time_t data = time(0);
    char tmp[32] = {""};
    strftime(tmp, sizeof(tmp), "%Y%m%d", localtime(&data));

    ofstream STEP1snvinfo_file;
    STEP1snvinfo_file.open(STEP1vcfInfo_path);

    ofstream snvhigh_file;
    snvhigh_file.open(highvcf_path);
    // snvhigh_file << "##fileformat=VCFv4.2" << endl;
    // snvhigh_file << "##fileData=" << tmp << endl;
    // snvhigh_file << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    // snvhigh_file << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
    // snvhigh_file << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    // snvhigh_file << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";

    // std :: map <int , refid> :: iterator itrefinfo = int2chrID.begin();
    // while (itrefinfo != int2chrID.end()) {;
    //     snvhigh_file << "##contig=<ID=" << itrefinfo->second.chrid << ",length=" << itrefinfo->second.length << ">\n";
    //     itrefinfo++;
    // }
    // snvhigh_file << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002" << endl;

    // body
    int snvpos_next = 1,pass = 0;
    for (auto oneSnv : ResultSnv) {
        if (oneSnv.ref_pos < opt->readbam_start || oneSnv.ref_pos > opt->readbam_end) continue;
        if (pass == 1) {
            snvpos_next = snvpos_next + 1;
            pass = 0;
            continue;
        }
        uint32_t nextpos = ResultSnv_refpos[snvpos_next];
        if (oneSnv.ref_pos == nextpos) {//the same pos has more than one vrt
            STEP1snvinfo_file << oneSnv.ref_pos << "\t" << oneSnv.REF << "\t" << oneSnv.ALT << "\t" << oneSnv.outinfo.QUAL << "\t0/1" << "\t" << oneSnv.outinfo.GQ << endl;
            STEP1snvinfo_file << oneSnv.ref_pos << "\t" << ResultSnv[snvpos_next].REF << "\t" << ResultSnv[snvpos_next].ALT << "\t" << ResultSnv[snvpos_next].outinfo.QUAL << "\t1/0" << "\t" << ResultSnv[snvpos_next].outinfo.GQ << endl;
            if(oneSnv.outinfo.GQ > 15){
                snvhigh_file << "chr" << give_chrint + 1 << "\t" << oneSnv.ref_pos << "\t"
                    << ".\t" << oneSnv.REF << "\t" << oneSnv.ALT << "\t" << oneSnv.outinfo.QUAL << "\t"
                    << "PASS"
                    << "\t"
                    << "."
                    << "\t"
                    << "GT:GQ"
                    << "\t"
                    << "0/1"
                    << ":" << oneSnv.outinfo.GQ << endl;
            } if(ResultSnv[snvpos_next].outinfo.GQ > 15){
                snvhigh_file << "chr" << give_chrint + 1 << "\t" << oneSnv.ref_pos << "\t"
                    << ".\t" << ResultSnv[snvpos_next].REF << "\t" << ResultSnv[snvpos_next].ALT << "\t" << ResultSnv[snvpos_next].outinfo.QUAL << "\t"
                    << "PASS"
                    << "\t"
                    << "."
                    << "\t"
                    << "GT:GQ"
                    << "\t"
                    << "0/1"
                    << ":" << ResultSnv[snvpos_next].outinfo.GQ << endl;
            }
            pass = 1;
            snvpos_next = snvpos_next + 1;

        } else {//0/1 or 1/1
            STEP1snvinfo_file << oneSnv.ref_pos << "\t" << oneSnv.REF << "\t" << oneSnv.ALT << "\t" << oneSnv.outinfo.QUAL << "\t" << oneSnv.outinfo.GT << "\t" << oneSnv.outinfo.GQ << endl;
            if((oneSnv.outinfo.GT == "0/1" || oneSnv.outinfo.GT == "1/0") && oneSnv.outinfo.GQ > 15){
                snvhigh_file << "chr" << give_chrint + 1 << "\t" << oneSnv.ref_pos << "\t"
                    << ".\t" << oneSnv.REF << "\t" << oneSnv.ALT << "\t" << oneSnv.outinfo.QUAL << "\t"
                    << "PASS"
                    << "\t"
                    << "."
                    << "\t"
                    << "GT:GQ"
                    << "\t" << oneSnv.outinfo.GT << ":" << oneSnv.outinfo.GQ << endl;
            }
            snvpos_next = snvpos_next + 1;
        }
    }
    ResultSnv_refpos.clear();
    STEP1snvinfo_file.close();
    snvhigh_file.close();
    return 0;
}


int printALLsnvvcf(char *SNV_PATH, char *STEP1vcfInfo_path, std :: vector <lowsnvVRT> &ResultLowsnv){
    /*
    load vcf from step1
    */
    ifstream ifs1vrtinfo;
    ifs1vrtinfo.open(STEP1vcfInfo_path, ios::in);
    std :: string S1vrtline;
    while (getline(ifs1vrtinfo , S1vrtline)) {
        lowsnvVRT add_one;
        stringstream bed_ss1(S1vrtline);
        string sp1;
        bed_ss1 >> sp1;
        add_one.ref_pos = stoi(sp1);
        bed_ss1 >> sp1;
        add_one.REF = sp1;
        bed_ss1 >> sp1;
        add_one.ALT = sp1;
        bed_ss1 >> sp1;
        add_one.outinfo.QUAL = stoi(sp1);
        bed_ss1 >> sp1;
        add_one.outinfo.GT = sp1;
        bed_ss1 >> sp1;
        add_one.outinfo.GQ = stoi(sp1);
        ResultLowsnv.push_back(add_one);
    }
    ifs1vrtinfo.close();

    sort(ResultLowsnv.begin() , ResultLowsnv.end());
    std :: vector <uint32_t> ResultSnv_refpos;
    for (auto cc : ResultLowsnv) ResultSnv_refpos.push_back(cc.ref_pos);
    ResultSnv_refpos.push_back(4294967294);
    sort(ResultSnv_refpos.begin() , ResultSnv_refpos.end());

    ofstream ALLsnvfile;
    ALLsnvfile.open(SNV_PATH);
    
    //header

    // ALLsnvfile << "##fileformat=VCFv4.2" << endl;
    // time_t data = time(0);
    // char tmp[32] = {""};
    // strftime(tmp, sizeof(tmp), "%Y%m%d", localtime(&data));
    // ALLsnvfile << "##fileData=" << tmp << endl;
    // ALLsnvfile << "##FILTER=<ID=PASS,Description=\"All filters passed\">\n";
    // ALLsnvfile << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n";
    // ALLsnvfile << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
    // ALLsnvfile << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
    // ALLsnvfile << "##ALT=<ID=DEL,Description=\"Deletion\">\n";
    // ALLsnvfile << "##ALT=<ID=INS,Description=\"Insertion\">\n";
    // std :: map <int , refid> :: iterator itrefinfo = int2chrID.begin();
    // while (itrefinfo != int2chrID.end()) {
    //     ALLsnvfile << "##contig=<ID=" << itrefinfo->second.chrid << ",length=" << itrefinfo->second.length << ">\n";
    //     itrefinfo++;
    // }
    // ALLsnvfile << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" << endl;

    // body
    int snvpos_next = 1,pass = 0;
    for(auto onelow : ResultLowsnv){
        if (onelow.ref_pos < opt->readbam_start || onelow.ref_pos > opt->readbam_end) continue;
        if (pass == 1) {
            snvpos_next = snvpos_next + 1;
            pass = 0;
            continue;
        }
        uint32_t nextpos = ResultSnv_refpos[snvpos_next];
        if (onelow.outinfo.GT == "0/0") { // the pos GT now
            snvpos_next = snvpos_next + 1;
            continue;
        }
        if (onelow.ref_pos == nextpos) {//the same pos has more than one vrt
            if(ResultLowsnv[snvpos_next].outinfo.GT == "0/0"){//print this snv
                ALLsnvfile << "chr" << give_chrint + 1 << "\t" << onelow.ref_pos << "\t"
                    << ".\t" << onelow.REF << "\t" << onelow.ALT << "\t" << onelow.outinfo.QUAL << "\t"
                    << "PASS"
                    << "\t"
                    << "."
                    << "\t"
                    << "GT:GQ"
                    << "\t" << onelow.outinfo.GT << ":" << onelow.outinfo.GQ << endl;
            }
            else{// 1/2
                ALLsnvfile << "chr" << give_chrint + 1 << "\t" << onelow.ref_pos << "\t"
                    << ".\t" << onelow.REF << "\t" << onelow.ALT << "\t" << onelow.outinfo.QUAL << "\t"
                    << "PASS"
                    << "\t"
                    << "."
                    << "\t"
                    << "GT:GQ"
                    << "\t" << "0/1" << ":" << onelow.outinfo.GQ << endl;
                ALLsnvfile << "chr" << give_chrint + 1 << "\t" << onelow.ref_pos << "\t"
                    << ".\t" << ResultLowsnv[snvpos_next].REF << "\t" << ResultLowsnv[snvpos_next].ALT << "\t" << ResultLowsnv[snvpos_next].outinfo.QUAL << "\t"
                    << "PASS"
                    << "\t"
                    << "."
                    << "\t"
                    << "GT:GQ"
                    << "\t" << "1/0" << ":" << ResultLowsnv[snvpos_next].outinfo.GQ << endl;
            }
            pass = 1;
            snvpos_next = snvpos_next + 1;
        } else{// 0/1 or 1/1
            ALLsnvfile << "chr" << give_chrint + 1 << "\t" << onelow.ref_pos << "\t"
                << ".\t" << onelow.REF << "\t" << onelow.ALT << "\t" << onelow.outinfo.QUAL << "\t"
                << "PASS"
                << "\t"
                << "."
                << "\t"
                << "GT:GQ"
                << "\t" << onelow.outinfo.GT << ":" << onelow.outinfo.GQ << endl;
            snvpos_next = snvpos_next + 1;
        }
    }
    ResultSnv_refpos.clear();
    ALLsnvfile.close();
    return 0;
}