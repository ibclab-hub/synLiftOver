#!/usr/bin/env python3

import sys
import re 
from collections import defaultdict

def readQuery(f_query):
    print("read query ...",file=sys.stderr)
    dict_query = {}
    for line in open(f_query):
        line = line.rstrip().split()
        if (line[0] in dict_query) == False:
            dict_query[line[0]] = []

        dict_query[line[0]].append([int(line[1]), int(line[2])])
    return dict_query

def readSize(f_tarsize):
    print("read targetSize ...",file=sys.stderr)
    dict_size = {}
    for line in open(f_tarsize):
        line = line.rstrip().split()
        dict_size[line[0]] = int(line[1])
    return dict_size 

def readSynteny(f_synteny):
    dict_refinfo = defaultdict(dict)
    dict_tarinfo = {}
    print("read synteny ...",file=sys.stderr)
    synteny_id = 0
    for line in open(f_synteny):
        line = line.rstrip().split()
        if (line == ""):
            continue
        synteny_id = line[0]
        refchr = line[1]
        dict_refinfo[refchr][synteny_id] = [int(line[2]), int(line[3])]
        dict_tarinfo[synteny_id] = [line[5], int(line[6]), int(line[7]), "+" if line[4] == line[8] else "-"]

    return dict_refinfo, dict_tarinfo


#def readAxt(chainpath,chromosome):
def readAxt(chromosome, path):
    print("read chain "+chromosome +" ...",file=sys.stderr)
    dict_tarseq = {}
    dict_refseq = {}
    dict_refStart = {}
    dict_refEnd = {}
    dict_tarChrom = {}
    dict_tarStart = {}
    dict_tarEnd = {}
    dict_tarDir = {}
    thisidx = 0
    flag_ref = 1
    for line in open(path+'/'+chromosome+'.axt'):
        line = line.rstrip()
        if (line == "") or (line[0] == "#") :
            continue

        if (line[0]).isdigit() == True:
            flag_ref = 1
            line = line.split()
            thisidx = int(line[0])
            dict_refStart[thisidx] = int(line[2])
            dict_refEnd[thisidx] = int(line[3])
            dict_tarChrom[thisidx] = line[4]
            dict_tarStart[thisidx] = int(line[5])
            dict_tarEnd[thisidx] = int(line[6])
            dict_tarDir[thisidx] = line[7]
            
        else:
            if flag_ref == 1:
                dict_refseq[thisidx] = line
                flag_ref = 0
            else:
                dict_tarseq[thisidx] = line

    return dict_tarseq, dict_refseq, dict_refStart, dict_refEnd, dict_tarChrom, dict_tarStart, dict_tarEnd, dict_tarDir  
