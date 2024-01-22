#!/usr/bin/env python
# coding=utf-8
''' 
 * @Title:  psi-pipe.py
 * @author: Cui Miao
 * @version V1.0   
'''
import os
import sys
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import logging
import argparse
from time import time

# import parameters as param
from utils import setup_logging ,PrintVCFHeader
from command_options import *
from check import *

#majorContigs = {"chr"+str(a) for a in list(range(0,23))}.union({str(a) for a in list(range(0,23))})
majorContigs = {"chr"+str(a) for a in list(range(1,23))}.union({str(a) for a in list(range(1,23))})

VERSION = 'V1.0.0'
PROG = "XXX"

class maindp(object):
	'''
	Detailed descriptions of SNP-caller version and its parameters.

	'''

	USAGE="""\
	Fast and accurate SNP for ONT Reads with %s.
	
	Current version: %s
	Author: Cui Miao
	Contact: XXXXX

	"""%(PROG, VERSION)

def run(args):
    basedir = os.path.dirname(__file__)
    callBin = basedir + '/Release/' + "VRT_ caller"

    fin_bam = CheckFileExist(args.fin_bam)
    fin_ref = CheckFileExist(args.fin_ref)
    fin_ref_fai = CheckFileExist(fin_ref, sfx=".fai")
    fin_bed = CheckFileExist(args.fin_bed)
    dup_bed = CheckFileExist(args.dup_bed)
    fin_index = CheckDirExist(args.fin_index)
    
    region_chunk_size = args.chunkWidth
    read_Ratio_ = args.read_Ratio
    thread_ = args.thread


    if not os.path.exists(args.workDir):
        os.makedirs(args.workDir)

    workDir = os.path.abspath(args.workDir)
    # workDir = basedir + '/' + args.workDir
    output_prefix = workDir + '/' if workDir[-1] != '/' else workDir

    header_fn = output_prefix + "TMPH_header.vcf"
    f_header = open(header_fn, 'w')
    PrintVCFHeader(f_header, fin_ref, args.sample)

    callVar_command_optionsS1 = [
        ExecuteCommand(callBin, "callerStep1"),
        CommandOption('fin_bam', fin_bam),
        CommandOption('fin_bed', fin_bed),
        CommandOption('dup_bed', dup_bed),
        CommandOption('fin_index', fin_index),
    ]
    callVar_command_optionsS3 = [
        ExecuteCommand(callBin, "callerStep3"),
        CommandOption('fin_bam', fin_bam),
        CommandOption('fin_bed', fin_bed),
        CommandOption('dup_bed', dup_bed),
        CommandOption('fin_index', fin_index),
    ]
    
    chrName = []
    if args.chrName is not None:
        chrName = [c.strip() for c in args.chrName.strip().split(',')]
    
    # 定义多个Bash脚本或者命令
    S1scripts = []
    S23scripts = []
    with open(fin_ref_fai, 'r') as fai_fp:
        for row in fai_fp:
            columns = row.strip().split("\t")
            contig_name = columns[0]
            if args.chrName is not None and contig_name not in chrName: 
                continue 

            if str(contig_name) not in majorContigs:
                continue

            tmp_dir = "_tmp_" + contig_name + "/"#_tmp_chr22/
            chr_outdir = output_prefix + tmp_dir

            mkdir_command ="mkdir " + chr_outdir
            result = subprocess.run(mkdir_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            # print(mkdir_command)

#          ---         STEP1         ---

            region_start, region_end = 0, 0
            contig_length = int(columns[1])
            while region_end < contig_length:
                region_start = region_end
                region_end = region_start + region_chunk_size
                if region_end > contig_length:
                    region_end = contig_length

                chunk_name = "_%s_%d_%d" % (contig_name, region_start, region_end)#chr22_0_100
                S1sh = "%s.sh" % (output_prefix+"TMPS1"+chunk_name)#TMP_chr22_0_50818468.sh
                S1sh_command = "bash "+S1sh
                S1scripts.append(S1sh_command)
                S1p2vrt_fn = chr_outdir+"tmpp2vrtinfo"+chunk_name
                S1p2vrtpos_fn = chr_outdir+"tmpp2vrtpos"+chunk_name
                S1s1vcf_fn = chr_outdir+"tmps1vcf"+chunk_name
                S1highvcf_fn = chr_outdir+"tmphighvcf"+chunk_name
                additional_command_optionsS1 = [
                    CommandOption('chr', contig_name),
                    CommandOption('start', region_start),
                    CommandOption('end', region_end),
                    CommandOption('tmp_outdir', chr_outdir),
                    CommandOption('p2vrt', S1p2vrt_fn),
                    CommandOption('p2vrtpos', S1p2vrtpos_fn),
                    CommandOption('s1vcf', S1s1vcf_fn),
                    CommandOption('highvcf', S1highvcf_fn),
                    CommandOption('read_Ratio', read_Ratio_),
                ]
                if contig_name == "chr6" :
                    additional_command_optionsS1 = [
                        CommandOption('chr', contig_name),
                        CommandOption('start', region_start),
                        CommandOption('end', region_end),
                        CommandOption('tmp_outdir', chr_outdir),
                        CommandOption('p2vrt', S1p2vrt_fn),
                        CommandOption('p2vrtpos', S1p2vrtpos_fn),
                        CommandOption('s1vcf', S1s1vcf_fn),
                        CommandOption('highvcf', S1highvcf_fn),
                        CommandOption('read_Ratio', 0),
                    ]
                with open(S1sh, 'w') as S1f:
                    print("#!/bin/bash", file=S1f)
                    print("",file=S1f)
                    print("\n#    ---- STEP1 ----\n" , file=S1f)
                    print(command_string_from(callVar_command_optionsS1) + " " + command_string_from(additional_command_optionsS1), file=S1f)


#          ---         STEP2&&3         ---

            S23sh = "%s.sh" % (output_prefix+"TMPS23_"+contig_name)#TMP_chr22_0_50818468.sh
            S23sh_command = "bash "+S23sh
            S23scripts.append(S23sh_command)
            S2highvcf_fn = chr_outdir+"Highsnv.vcf"
            S3p2vrt_fn = chr_outdir+"P2vrtInfo"
            S3p2vrtpos_fn = chr_outdir+"P2vrtposInfo"
            S3s1vcf_fn = chr_outdir+"S1vcfInfo"
            S3phasedtsv_fn = chr_outdir+"phasedinfo.tsv"
            output_fn = "%s.vcf" % (output_prefix+"TMPV_"+contig_name)
            additional_command_options3 = [
                CommandOption('chr', contig_name),
                CommandOption('tmp_outdir', chr_outdir),
                CommandOption('p2vrt', S3p2vrt_fn),
                CommandOption('p2vrtpos', S3p2vrtpos_fn),
                CommandOption('s1vcf', S3s1vcf_fn),
                CommandOption('phased_tsv', S3phasedtsv_fn),
                CommandOption('fout_vcf', output_fn),
            ]

            with open(S23sh, 'w') as S23f:
                print("#!/bin/bash", file=S23f)
                print("",file=S23f)

                print("cat " + chr_outdir+"tmpp2vrtinfo* > " + S3p2vrt_fn , file=S23f)
                print("cat " + chr_outdir+"tmpp2vrtpos* > " + S3p2vrtpos_fn , file=S23f)
                print("cat " + chr_outdir+"tmps1vcf* > " + S3s1vcf_fn , file=S23f)
                print("cat " + chr_outdir+"tmphighvcf* > " + chr_outdir+"tmpHighbody.vcf" , file=S23f)
                print("sort -k1,1 -k2,2g " + chr_outdir+"tmpHighbody.vcf > " + chr_outdir+"tmpHighbody_sort.vcf",file=S23f)
                print("cat " + header_fn + " " + chr_outdir+"tmpHighbody_sort.vcf > " + S2highvcf_fn,file=S23f)
                print("bgzip " + S2highvcf_fn, file=S23f)
                print("bcftools index " + chr_outdir+"Highsnv.vcf.gz", file=S23f)
                print("rm -r "+chr_outdir+"tmp*", file=S23f)
                print("",file=S23f)
            #     STEP2 whatshap
                print("\n#    ---- STEP2 ----\n" , file=S23f)
                print("whatshap phase " + "--reference=" + fin_ref +" " +chr_outdir+"Highsnv.vcf.gz " + fin_bam + " --chromosome " + contig_name + " --distrust-genotypes --ignore-read-groups -o " + chr_outdir+"phased.vcf.gz ", file=S23f)
                print("bcftools index " + chr_outdir+"phased.vcf.gz", file=S23f)
                print("whatshap haplotag " + "--reference=" + fin_ref +" " +chr_outdir+"phased.vcf.gz " + fin_bam + " --regions " + contig_name + " --ignore-read-groups -o " + chr_outdir+"tag.bam " + "--output-haplotag-list="+S3phasedtsv_fn, file=S23f)
                print("",file=S23f)
            #     STEP3 
                print("\n#    ---- STEP3 ----\n" , file=S23f)
                print(command_string_from(callVar_command_optionsS3) + " " + command_string_from(additional_command_options3), file=S23f)
                # print("rm -r "+chr_outdir, file=S23f)
    

#run_part
    # Step 1 run and wait finished
    t_s1 = time()
    with ThreadPoolExecutor(max_workers=thread_) as executor:
        futures = [executor.submit(subprocess.run, script, shell=True) for script in S1scripts]
        for future in as_completed(futures):
            try:
                result = future.result()
                # print(f"FINISHed : {result.returncode}\n")
            except Exception as e:
                print(f"caller exited with exceptions. Exiting... : {e}")
    rm_TMPS1 ="rm -r " + output_prefix+"TMPS1*"
    result = subprocess.run(rm_TMPS1, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    t_e1 = time()
    logging.info("[chr_chunk.py] Finish SNV calling step1 in %.2f seconds\n" %(t_e1 - t_s1))

# Step 2&&3 run and wait finished
    t_s23 = time()
    with ThreadPoolExecutor(max_workers=thread_) as executor:
        futures = [executor.submit(subprocess.run, script, shell=True) for script in S23scripts]
        for future in as_completed(futures):
            try:
                result = future.result()
                # print(f"FINISHed : {result.returncode}\n")
            except Exception as e:
                print(f"caller exited with exceptions. Exiting... : {e}")
    t_e23 = time()
    logging.info("[chr_chunk.py] Finish SNV calling step2 and step3 in %.2f seconds\n" %(t_e23 - t_s23))

#    merge vcf's body
    catbody_command ="cat " + output_prefix+"TMPV_* > " + output_prefix+"TMP_body.vcf"
    result = subprocess.run(catbody_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # print(catbody_command)
#    sort vcf's body
    sort_command = "sort -k1,1 -k2,2g " + output_prefix+"TMP_body.vcf > " + output_prefix+"TMP_sort.vcf"
    result = subprocess.run(sort_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # print(sort_command)
#    merge vcf's head and sorted body
    cathAb_command ="cat " + header_fn + " " + output_prefix+"TMP_sort.vcf " + "> " + output_prefix+"TMP_all.vcf"
    result = subprocess.run(cathAb_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # print(cathAb_command)

    bgzip_command ="bgzip " + output_prefix+"TMP_all.vcf"
    result = subprocess.run(bgzip_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # print(bgzip_command)

    norm_command ="bcftools norm -m+ " + output_prefix+"TMP_all.vcf.gz -Oz -o "+ output_prefix+"merge.vcf.gz"
    result = subprocess.run(norm_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # print(norm_command)

    tabix_command ="tabix " + output_prefix+"merge.vcf.gz"
    result = subprocess.run(tabix_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # print(tabix_command)

    rm_command ="rm -r " + output_prefix+"TMP*"
    result = subprocess.run(rm_command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    # print(rm_command)




def main():
    parser = argparse.ArgumentParser(prog="python chr_chunk_task.py",
            description=maindp.USAGE,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--version', '-v', 
		action = 'version', 
		version = '%(prog)s {version}'.format(version=VERSION))

    # **************Parameters of input******************

    parser.add_argument("--fin_ref",
            type=str,
            required=True,
            help="The reference genome in fasta format. [FASTA]")
    
    parser.add_argument("--fin_bam", 
            type=str, 
            required=True,
            help="Sorted .bam file fromMinimap2. [BAM]")
    
    parser.add_argument("--fin_index", 
            type=str, 
            required=True,
            help="The path of RdBG index. [DIR_PATH]")
    
    parser.add_argument("--fin_bed", 
            type=str, 
            required=True,
            help="HOMO Bed format input. [BED]")
    
    parser.add_argument("--dup_bed", 
            type=str, 
            required=True,
            help="DUP HOMO Bed format input. [BED]")

    parser.add_argument("--workDir",
            type=str,
            required=True,
            help="Work-directory for distributed job")
    
    # ************** Other Parameters******************
    
    parser.add_argument("--chrName", 
            type=str, 
            default=None,
            help="list of chrName (contig names) to be processed, separeted by comma without any blank space. [default : %(default)s]")

    parser.add_argument('--sample',
            type = str,
            default = "SAMPLE",
            help = "Sample name/id. [default : %(default)s]")
    
    parser.add_argument("--chunkWidth", 
            type=int, 
            default=10000000,
            help="Reference length to detect candidate in one loop. [default : %(default)d]")
    
    parser.add_argument("--read_Ratio", 
            type=float, 
            default=0.98,
            help="read_Ratio. [default : %(default)d]")
    
    parser.add_argument("--thread", 
            type=int, 
            default=16,
            help="Number of threads to use. [default : %(default)d]")

    args = parser.parse_args()
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    setup_logging()
    t_s = time()
    run(args)
    t_e = time()
    logging.info("[calling.py] Finish SNV calling of genomic in %.5f seconds" %(t_e - t_s))

#     run(args)

if __name__ == "__main__":
    main()
