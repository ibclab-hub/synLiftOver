#!/usr/bin/env python3

import sys
import os
import re
import fileReader as fR
import filterFiles as fF
import mergeAxt as mA
import alloc
import subprocess
from collections import defaultdict


filter_flag = 0
merge_flag = 1

if len(sys.argv) != 8:
    sys.exit("Usage: preprocess_axtFiltMerge.py synteny.bed axt_dir output_dir gap filter merge target.size")
f_synteny = sys.argv[1]
path_axt = sys.argv[2]
path_newaxt = sys.argv[3]
min_gapSize = int(sys.argv[4])
filter_flag = int(sys.argv[5])
merge_flag = int(sys.argv[6])
target_sizes = fR.readSize(sys.argv[7])
if any(size <= 0 for size in target_sizes.values()):
    sys.exit("ERROR: chromosome sizes must be positive")

#make directory for axt generation
subprocess.call(['mkdir',path_newaxt])

#read Synteny
dict_refInfo, dict_tarInfo = fR.readSynteny(f_synteny)

axt_files = os.listdir(path_axt)
chromosome_list = [f.split('.')[0] for f in axt_files if f.endswith('.axt')]

#filter and merge 
for chromosome in chromosome_list:
    f = open(path_newaxt+"/"+chromosome+".axt","w")
    d_tarSeq, d_refSeq, d_refSt, d_refEd, d_tarchr, d_tarSt, d_tarEd, d_tarDir = fR.readAxt(chromosome, path_axt)
    filtered_idx, chainRelatedSynId = fF.filterAxt(d_refSt, d_refEd, d_tarSt, d_tarEd, dict_refInfo[chromosome],d_tarDir, dict_tarInfo, filter_flag, d_tarchr, target_sizes)
    
    curAxtInfo = {}
    mergeAxtInfo = {}
    mergeAxtInfo["idx"] = -1
    j = -1
    mergeExistFlag = 0

    for i in filtered_idx:
        print("curidx:",i, file=sys.stderr)
        #set every axt
        curAxtInfo["idx"] = i
        curAxtInfo["refseq"] = d_refSeq[i]
        curAxtInfo["refchr"] = chromosome
        curAxtInfo["tarseq"] = d_tarSeq[i]
        curAxtInfo["tarchr"] = d_tarchr[i]
        curAxtInfo["refSt"] = d_refSt[i]
        curAxtInfo["refEd"] = d_refEd[i]
        curAxtInfo["tarSt"] = d_tarSt[i]
        curAxtInfo["tarEd"] = d_tarEd[i]
        curAxtInfo["tardir"] = d_tarDir[i]
        mergeFailFlag = 0
        
        if (mergeExistFlag == 0):
            mergeAxtInfo = alloc.allocation(curAxtInfo)
            mergeExistFlag=1
            j = i
            continue
        else:
            #merge Exists
            #check if merge is available
            #1 connected idx 
            if(merge_flag != 0):
                if ((i - 1) != j ):
                    print("failed with connected chain", j , i ,file=sys.stderr)
                    mergeFailFlag=1
                #2 sameSynId check
                elif (chainRelatedSynId[i] == -1 or chainRelatedSynId[j] == -1 or chainRelatedSynId[i] != chainRelatedSynId[j]):
                    print("failed with synteny",i,j, file=sys.stderr)
                    mergeFailFlag=1
                #3 chromosome check 
                elif (d_tarchr[i] != d_tarchr[j]):
                    print("failed with chr",i,j, file=sys.stderr)
                    mergeFailFlag=1
                #4 dir check 
                elif (d_tarDir[i] != d_tarDir[j]):
                    print("failed with dir",i,j, file=sys.stderr)
                    mergeFailFlag=1
                #4.1 order check 
                elif (d_tarSt[j] > d_tarEd[i]):
                    print("failed with blockorder",i,j, file=sys.stderr)
                    mergeFailFlag=1

                #5 gapSize and REF/TAR check 
                refDiff = d_refSt[i] - d_refEd[j] - 1
                tarDiff = d_tarSt[i] - d_tarEd[j] - 1
                if not (refDiff == 0  and tarDiff >= 0 and tarDiff < min_gapSize):
                    print("failed with gapSize",i,j,refDiff, tarDiff, file=sys.stderr)
                    mergeFailFlag = 1

            else:
                mergeFailFlag=1
            
            if (mergeFailFlag == 1):
                #merge is not available
                mA.printAxt(mergeAxtInfo,f)
                #set merge with cur  
                mergeAxtInfo = {}
                mergeAxtInfo = alloc.allocation(curAxtInfo)
                mergeExistFlag = 1
                j = i
                continue
            else:
                #merge is available
                #do merge 
                print("can merge",j,i, file=sys.stderr)
                mergeAxtInfo["idx"] = i
                #mergeAxtInfo["refSt"] = min(d_refSt[j], d_refSt[i])
                #mergeAxtInfo["refEd"] = max(d_refEd[j], d_refEd[i])
                #mergeAxtInfo["tarSt"] = min(d_tarSt[j], d_tarSt[i])
                #mergeAxtInfo["tarEd"] = max(d_tarEd[j], d_tarEd[i])

                mergeAxtInfo["refEd"] = d_refEd[i]
                mergeAxtInfo["tarEd"] = d_tarEd[i]
                mergeAxtInfo["refseq"] += "-" * (refDiff+tarDiff) + d_refSeq[i]
                mergeAxtInfo["tarseq"] += "N" * (refDiff+tarDiff) + d_tarSeq[i]
                mergeExistFlag = 1
                #print ("add merge ", i)
                j = i
                continue
    #After for loop process left Axt to print 
    if (mergeAxtInfo["idx"] != -1):
        mA.printAxt(mergeAxtInfo,f)
