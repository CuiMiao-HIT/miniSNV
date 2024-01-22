#!/usr/bin/env python
# coding=utf-8

import sys
import shlex
import logging
from subprocess import check_output, PIPE, Popen
from check import CheckFileExist
from collections import namedtuple
# from intervaltree import IntervalTree

majorContigs = {"chr"+str(a) for a in list(range(0,23))+["X", "Y"]}.union({str(a) for a in list(range(0,23))+["X", "Y"]})
chrName_to_dict = {"chr1":0, "chr2":1, "chr3":2, "chr4":3, "chr5":4, "chr6":5, "chr7":6, "chr8":7, "chr9":8, "chr10":9, 
                   "chr11":10, "chr12":11, "chr13":12, "chr14":13, "chr15":14, "chr61":15, "chr17":16, "chr18":17, "chr19":18, "chr20":19,
                   "chr21":20, "chr22":21}

base_to_ACGT = dict(zip(
    "ACGTURYSWKMBDHVN",
    ("A", "C", "G", "T", "T", "A", "C", "C", "A", "G", "A", "C", "A", "A", "A", "A")
))

def PrintGVCFHeader(file, fin_ref, sample):
    contigINFO = list() 
    fai_fn = CheckFileExist(fin_ref + ".fai")

    with open(fai_fn, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if line[0] in majorContigs:
                contigINFO.append([line[0], int(line[1])])

    file.write('##fileformat=VCFv4.2\n')

    for i in contigINFO:
        file.write("##contig=<ID=%s,length=%d>\n"%(i[0], i[1]))

    file.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    file.write('##FILTER=<ID=FILTER,Description="Genotyping model thinks this site is reference or lower than a threshold score.">\n')
    #file.write('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">\n')
    file.write('##INFO=<ID=CR,Number=1,Type=String,Description="Variant called from type of candidate, one of HC(high confidence candidate), LC(low confidence candiadte) and TRC(tandem repeat candidate)">\n')
    file.write('##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n')
    file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    file.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
    file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
    file.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read Depth for each allele">\n')
    file.write('##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">\n')
    file.write('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">\n')

    file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % (sample))

def PrintVCFHeader(file, fin_ref, sample):
    contigINFO = list() 
    fai_fn = CheckFileExist(fin_ref + ".fai")

    with open(fai_fn, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if line[0] in majorContigs:
                contigINFO.append([line[0], int(line[1])])

    file.write('##fileformat=VCFv4.2\n')

    for i in contigINFO:
        file.write("##contig=<ID=%s,length=%d>\n"%(i[0], i[1]))

    file.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
    file.write('##FILTER=<ID=FILTER,Description="Genotyping model thinks this site is reference or lower than a threshold score.">\n')
    file.write('##INFO=<ID=CR,Number=1,Type=String,Description="Variant called from type of candidate, one of HC(high confidence candidate), LC(low confidence candiadte) and TRC(tandem repeat candidate)">\n')
    file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    file.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')

    file.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n' % (sample))
    file.close()


# def Generate_VCFOutput(file, Result):
#     idx = 0
#     prePos = 0
#     preTyp = "INS"
#     preAlt = ""
#     for item in Result:
#         if item.ALT == "NON_REF": continue

#         GQ = item.GQ
#         if item.qname == "HC":
#             item = [item.chr, int(item.pos), item.type, item.REF, item.ALT, item.qname, item.MQ, item.GT, item.varlen, item.qual, item.supportA, item.supportT, item.supportT-item.supportA, item.supportA, item.supportT]
#         else:
#             if item.GT == "1/2":
#                 item = [item.chr, int(item.pos), item.type, item.REF, item.ALT, item.qname, item.MQ, item.GT, item.varlen, item.qual, item.supportA, item.supportT, item.con1, item.con2, item.read_T]
#             else:
#                 item = [item.chr, int(item.pos), item.type, item.REF, item.ALT, item.qname, item.MQ, item.GT, item.varlen, item.qual, item.supportA, item.supportT, item.con2, item.con1, item.read_T]
#         if item[1] == prePos and item[2] == preTyp and item[4] == preAlt:
#             prePos = item[1]
#             preTyp = item[2]
#             preAlt = item[4]
#             continue
#         info_list = "CR={CON}".format(
#             CON = item[5])
#         AD = str(item[12])  + "," + str(item[13])
#         file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t{SAMPLE}\n".format(
#             CHR = item[0],
#             POS = item[1],
#             ID = ".",   #ID = "MLSV.%d"%(idx),
#             REF = item[3],
#             ALT = item[4],
#             QUAL = item[9], #item[9],
#             FILTER = "FILTER" if item[7] == "0/0" else "PASS",
#             INFO = info_list,
#             FORMAT = "GT:GQ:DP:AD",
#             SAMPLE = item[7] + ":" + str(GQ) + ":" + str(item[14]) + ":" + AD))
#         prePos = item[1]
#         preTyp = item[2]
#         preAlt = item[4]
#         idx += 1

# def Generate_GVCFOutput(file, Result, gq_binsize):
#     Result = sorted(Result, key = lambda x:int(x.pos))
#     prePos = 0; preTyp = "INS"; preAlt = ""
#     MIN_GQ  = 0; MIN_DP = 0; START = 0;  END = 0; Lref = 0; PL = []
#     for item in Result:
#         #print(item.chr,item.pos, item.REF, item.ALT, item.GQ)
#         if item.ALT == "NON_REF":  # non-reference
#             if preAlt != "NON_REF":
#                 binrange = int(item.GQ / gq_binsize); START = item.pos; END = item.pos; Lref = item.REF
#                 MIN_GQ = item.GQ; MIN_DP = item.supportT; preAlt = item.ALT; PL = item.PL
#                 continue
#             if int(item.GQ / gq_binsize) == binrange:
#                 END = item.pos; preAlt = item.ALT
#                 MIN_GQ = MIN_GQ if MIN_GQ < item.GQ else item.GQ
#                 MIN_DP = MIN_DP if MIN_DP < item.supportT else item.supportT
#             else:
#                 ## output previous block
#                 info_list = "END={ed}".format(ed=END)
#                 file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t{SAMPLE}\n".format(
#                     CHR = item.chr,
#                     POS = START,
#                     ID = ".",   #ID = "MLSV.%d"%(idx),
#                     REF = Lref,
#                     ALT = "<NON_REF>",
#                     QUAL = 0, #item[9],
#                     FILTER = ".",
#                     INFO = info_list,
#                     FORMAT = "GT:GQ:MIN_DP:PL",
#                     SAMPLE = "0/0" + ":" + str(MIN_GQ) + ":" + str(MIN_DP) + ":" + ",".join([str(i) for i in PL])))
#                 ## current value
#                 binrange = int(item.GQ / gq_binsize); START = item.pos; END = item.pos; Lref = item.REF; PL = item.PL
#                 MIN_GQ = item.GQ; MIN_DP = item.supportT; preAlt = item.ALT
#         else: ## have variants
#             if preAlt == "NON_REF":
#                 ## output previous block
#                 info_list = "END={ed}".format(ed=END)
#                 file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t{SAMPLE}\n".format(
#                 CHR = item.chr,
#                 POS = START,
#                 ID = ".",   #ID = "MLSV.%d"%(idx),
#                 REF = Lref,
#                 ALT = "<NON_REF>",
#                 QUAL = 0, #item[9],
#                 FILTER = ".",
#                 INFO = info_list,
#                 FORMAT = "GT:GQ:MIN_DP:PL",
#                 SAMPLE = "0/0" + ":" + str(MIN_GQ) + ":" + str(MIN_DP) + ":" + ",".join([str(i) for i in PL])))

#             GQ = item.GQ; GT = item.GT; PL = item.PL; AD = item.AD + ",0"
#             if item.qname == "HC":
#                 item = [item.chr, int(item.pos), item.type, item.REF, item.ALT, item.qname, item.MQ, item.GT, item.varlen, item.qual, item.supportA, item.supportT, item.supportT-item.supportA, item.supportA, item.supportT]
#             else:
#                 if item.GT == "1/2":
#                     item = [item.chr, int(item.pos), item.type, item.REF, item.ALT, item.qname, item.MQ, item.GT, item.varlen, item.qual, item.supportA, item.supportT, item.con1, item.con2, item.read_T]
#                 else:
#                     item = [item.chr, int(item.pos), item.type, item.REF, item.ALT, item.qname, item.MQ, item.GT, item.varlen, item.qual, item.supportA, item.supportT, item.con2, item.con1, item.read_T]
#             if item[1] == prePos and item[2] == preTyp and item[4] == preAlt:
#                 prePos = item[1]
#                 preTyp = item[2]
#                 preAlt = item[4]
#                 continue
#             info_list = "CR={CON}".format(
#                 CON = item[5])
#             #if GT == "1/2" :
#             #    AD = str(item[14]-item[12]-item[13]) + "," + str(item[12])  + "," + str(item[13]) + ",0"
#             #else:
#             #    AD = str(item[12])  + "," + str(item[13]) + ",0"
#             file.write("{CHR}\t{POS}\t{ID}\t{REF}\t{ALT}\t{QUAL}\t{FILTER}\t{INFO}\t{FORMAT}\t{SAMPLE}\n".format(
#                 CHR = item[0],
#                 POS = item[1],
#                 ID = ".",   #ID = "MLSV.%d"%(idx),
#                 REF = item[3],
#                 ALT = item[4]+",<NON_REF>",
#                 QUAL = item[9], #item[9],
#                 FILTER = "FILTER" if item[7] == "0/0" else "PASS",
#                 INFO = info_list,
#                 FORMAT = "GT:GQ:DP:AD:PL",
#                 SAMPLE = item[7] + ":" + str(GQ) + ":" + str(item[14]) + ":" + AD + ":" + ",".join([str(i) for i in PL])))
#             prePos = item[1]
#             preTyp = item[2]
#             preAlt = item[4]

def setup_logging():
    """
    Default logger
    """
    log_format = "%(asctime)s [%(levelname)s] %(message)s"
    logging.basicConfig(level=logging.INFO, format=log_format)

def get_contig_length(fp_fai):
    Contig_dict = dict() 
    with open(fp_fai, 'r') as fai_fp:
        for row in fai_fp:
            columns = row.strip().split("\t")
            contig_name = columns[0]
            contig_length = int(columns[1])
            Contig_dict[contig_name] = contig_length 
    
    return Contig_dict

def reset_Start_End(chrStart, chrEnd, chrEnd_max, flanking, is_range_given):
    if is_range_given:
        rS = chrStart - flanking; rE = chrEnd + flanking # for debug
        rS = 1 if rS < 1 else rS
        rE = rE if rE < chrEnd_max else chrEnd_max
    else:
        rS= 1; rE = None

    return rS, rE

def evc_base_from(base):
    return base if base == "N" else base_to_ACGT[base]

def subprocess_popen(args, stdin=None, stdout=PIPE, stderr=sys.stderr, bufsize=8388608):
    return Popen(args, stdin=stdin, stdout=stdout, stderr=stderr, bufsize=bufsize, universal_newlines=True)

def get_reference_sequence(ref_path, chrName, chrStart, chrEnd):
    refernce_sequences = []

    if chrStart != None and chrEnd != None:
        region_value = "{}:{}-{}".format(chrName, chrStart, chrEnd)
    else:
        region_value = str(chrName)
    
    samtools_faidx_process = subprocess_popen(
        shlex.split("samtools faidx {} {}".format(ref_path, region_value))
    )

    while True:
        row = samtools_faidx_process.stdout.readline()
        is_finish_reading_output = row == '' and samtools_faidx_process.poll() is not None
        if is_finish_reading_output:
            break
        if row:
            refernce_sequences.append(row.rstrip())

    # first line is reference name ">xxxx", need to be ignored
    reference_sequence = "".join(refernce_sequences[1:])

    # uppercase for masked sequences
    reference_sequence = reference_sequence.upper()

    samtools_faidx_process.stdout.close()
    samtools_faidx_process.wait()
    if samtools_faidx_process.returncode != 0:
        return None

    return reference_sequence

def get_sample_count(bamList):
    bam_file_path = list()
    sample_name_list = list()
    sampleCnt = 0
    with open(bamList, 'r') as f:
        for line in f:
            sample_name, sample_path = line.strip().split('\t')
            bam_file_path.append(sample_path)
            sample_name_list.append(sample_name)
            sampleCnt += 1

    return sampleCnt, bam_file_path, sample_name_list

#def overlap(s1, e1, s2, e2):
#    S = max(s1, s2)
#    E = min(e1, e2)
#
#    if S <= E: return True
#    
#    return False
def is_homopolymer(refSeq, refStart, start, end):
    reference = refSeq[start-refStart:end-refStart]
    homo_cnt = 0
    min_homo_length = 9
    for i in range(1, len(reference)):
        if reference[i] == reference[i-1]:
            homo_cnt += 1
        else:
            if homo_cnt >= min_homo_length:
                return True
            else:
                homo_cnt = 0


    return True if homo_cnt >= min_homo_length else False


# def get_repeat_region(repeat_file_path, ctgName, ctgStart, ctgEnd, refSeq, refStart, parser_homo=False):
#     repeatMask = IntervalTree()
#     with open(repeat_file_path, 'r') as f:
#         for line in f:
#             line = line.strip().split('\t')
#             if line[0] != ctgName: continue
#             s = int(line[1]); e = int(line[2])
#             if e - s > 500: continue
#             if e < ctgStart: continue
#             if s > ctgEnd: break
#             if parser_homo:
#                 is_homo = is_homopolymer(refSeq, refStart, s, e)
#                 repeatMask.addi(s, e, is_homo)
#             else:
#                 repeatMask.addi(s, e)


#     return repeatMask

def chrName_to_chrID(chrName):
    return chrName_to_dict[chrName]

OutputMethods = namedtuple('OutputMethods', [
    'output',
    'output_header',
    'close_opened_files',
    'gen_output_file'
])

def output_methods_from(sample_names,
                        reference_file_path,
                        output_file_path):

    def gen_output_file():
        global output_file
        output_file = open(output_file_path, "w")

    def output(string_value):
        global output_file
        print(string_value, file=output_file)

    def close_opened_files():
        output_file.close()

    def output_header():
        from textwrap import dedent
        output(dedent("""\
            ##fileformat=VCFv4.2
            ##FILTER=<ID=PASS,Description="All filters passed">
            ##FILTER=<ID=FILTER,Description="Genotyping model thinks this site is reference or lower than a threshold score">
            ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
            ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
            ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
            ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Read depth for each allele">
            ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to the closest integer">"""
                      ))

        reference_index_file_path = reference_file_path + ".fai"
        with open(reference_index_file_path, "r") as fai_fp:
            for row in fai_fp:
                columns = row.strip().split("\t")
                contig_name, contig_size = columns[0], columns[1]
                output("##contig=<ID=%s,length=%s>" % (contig_name, contig_size))

        output('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s' % ("\t".join(sample_names)))

    return OutputMethods(
        output,
        output_header,
        close_opened_files,
        gen_output_file
    )
