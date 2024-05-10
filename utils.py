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


def PrintVCFHeader(file, fin_ref, sample, human):
    contigINFO = list() 
    fai_fn = CheckFileExist(fin_ref + ".fai")

    with open(fai_fn, 'r') as fin:
        for line in fin:
            line = line.strip().split("\t")
            if human:
                if line[0] in majorContigs:
                    contigINFO.append([line[0], int(line[1])])
            else:
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
