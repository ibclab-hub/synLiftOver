#!/usr/bin/env python3

import sys
import re 
from collections import defaultdict


def printAxt(axtInfo,fo):
    axt_line = str(axtInfo["idx"])+"\t"+axtInfo["refchr"]+"\t"+str(axtInfo["refSt"])+"\t"+str(axtInfo["refEd"])+"\t"+axtInfo["tarchr"]+"\t"+str(axtInfo["tarSt"])+"\t"+str(axtInfo["tarEd"])+"\t"+axtInfo["tardir"]+"\t0\n"
    if (len(axtInfo["refseq"]) != len(axtInfo["tarseq"])):
        print("error !! axtInfo: ", axtInfo)
    fo.write(axt_line)
    fo.write(axtInfo["refseq"]+"\n")
    fo.write(axtInfo["tarseq"]+"\n")
    nulldict = {}
    return nulldict

def printQuery(query,fo):
    for chromosome in query.keys():
        for idx in query[chromosome]:
            fo.write(chromosome+"\t"+str(idx[0])+"\t"+str(idx[1])+"\n")
    return 1