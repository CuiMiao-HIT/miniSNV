#!/usr/bin/env python
# coding=utf-8

import os
import sys
from subprocess import check_output

def CheckBamExist(fin_bam_list, sfx=""):
    with open(fin_bam_list) as f:
        for row in f:
            fn = row.strip().split('\t')[1]
            if not os.path.isfile(fn+sfx):
                sys.exit("Error: %s not found" %(fn+sfx))

def CheckDirExist(fn):
    if not os.path.exists(fn):
        sys.exit("Error: %s not found" %fn)
    return os.path.abspath(fn)

def CheckFileExist(fn, sfx=""):
    if not os.path.isfile(fn+sfx):
        sys.exit("Error: %s not found" %(fn+sfx))
    return os.path.abspath(fn+sfx)


def CheckCmdExist(cmd):
    if not isinstance(cmd, str):
        return False

    try:
        check_output("which %s" % (cmd), shell=True)
    except:
        sys.exit("Error: %s executable not found" %(cmd))
    return cmd
