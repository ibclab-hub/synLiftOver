import sys
import re 
from collections import defaultdict

def filterQuery(d_query, d_synteny,):
    print("filtering query",file=sys.stderr)
    filtered_query = defaultdict()
    for chromosome in d_query.keys():
        q_list = []
        for this_query in d_query[chromosome] :
            q_start = int(this_query[0])
            q_end = int(this_query[1])
            for synId in d_synteny[chromosome]:
                [s_start,s_end] = d_synteny[chromosome][synId]
                if(s_start <= q_start and q_end <= s_end):
                    #in the synteny
                    q_list.append(this_query)
                    break
        filtered_query[chromosome] = q_list
    return filtered_query

def filterAxt(d_refst, d_refed, d_tarst, d_tared ,d_synteny_chr, d_axt_dir, d_syn_tar, filter_flag, d_tarchr, target_sizes): #for single chromosome 
    print("filtering chain",file=sys.stderr)
    filtered_idxList = []
    #chainIdxMatchSynId = []
    chainIdxMatchSynId = {}
    for idx in d_refst.keys():
        # Compare 1-based inclusive AXT with forward 0-based half-open synteny.
        # Do not modify the original coordinates used for AXT output/merge.
        axt_st = d_refst[idx] - 1
        axt_ed = d_refed[idx]
        axt_t_st = d_tarst[idx] - 1
        axt_t_ed = d_tared[idx]
        dir = d_axt_dir[idx]
        if dir == '-':
            chrom = d_tarchr[idx]
            if chrom not in target_sizes:
                raise ValueError("Missing target chromosome size for reverse AXT: " + chrom)
            size = target_sizes[chrom]
            if not (0 <= axt_t_st < axt_t_ed <= size):
                raise ValueError("Invalid reverse AXT target bounds: " + chrom)
            axt_t_st, axt_t_ed = size - axt_t_ed, size - axt_t_st
        chainIdxMatchSynId[idx] = -1
        if(filter_flag == 0):
            filtered_idxList.append(idx)
        
        for synId in d_synteny_chr:
            [s_start,s_end] = d_synteny_chr[synId]
            s_tar_start = d_syn_tar[synId][1]
            s_tar_end = d_syn_tar[synId][2]
            if d_tarchr[idx] != d_syn_tar[synId][0]:
                continue
            if (dir != d_syn_tar[synId][3]): #axt direction different from synteny direction
                continue
            if (s_start <= axt_st and axt_ed <= s_end):
                if (s_tar_start <= axt_t_st and axt_t_ed <= s_tar_end):
                    #axt in the synteny 
                    if(filter_flag != 0 ):
                        filtered_idxList.append(idx)
                    
                    chainIdxMatchSynId[idx] = synId
                    break
    
    return filtered_idxList, chainIdxMatchSynId
