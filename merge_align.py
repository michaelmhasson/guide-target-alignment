import csv
import numpy as np
from string import maketrans
import itertools
import affalign
def revcom(s):# reverse complement function
    s = s[::-1]
    s = s.translate(maketrans("ATCG","TAGC"))
    return s
gdict = {}
guides = ['GGACCCTTCGTCGACACGAT','GGCATTGTCGATTCGCTCTC','GCGCGAAAATTCTATGAGAA','GATATCGGAGTTGTTGTTCG','GGTTTGACCGTACTTGCAGC','GCGGAATCTCCTGTCCGAGA','GTCGACGATTTTCCACAACC','GGCGATCGACTTACGAATAC','GTCCCTGCACGATACGCATA','GCGGCTACCTATTAACGTAG','GTCGGTTTATTTGACGAACG','GTCGGGTCATTCGTCGATTG','GACCAGGATGGGCACCACCC']
with open('highO1.fasta') as f1, open('highO2.fasta') as f2: # parse through both data files line by line
    for s1, s2 in itertools.izip_longest(f1, f2):
        if s1[0] != '>':
            s1 = s1.strip()
            s2 = s2.strip()
            if len(s1) == 20: # 20 length sequence is the guide
                s2 = revcom(s2)
                if s1 in gdict:# put in nested dicts with
                    if s2 in gdict[s1]:
                        gdict[s1][s2] += 1
                    else:
                        gdict[s1][s2] = 1
                else:
                    gdict[s1] = {}
                    gdict[s1][s2] = 1
            elif len(s2) == 20:
                s1 = revcom(s1)
                if s2 in gdict:
                    if s1 in gdict[s2]:
                        gdict[s2][s1] += 1
                    else:
                        gdict[s2][s1] = 1
                else:
                    gdict[s2] = {}
                    gdict[s2][s1] = 1
pre = 'GCTTTTTTTCTCGAGTACTA'
suff = 'AGGATCCATTAGGCGGCCGC'
for guide in gdict:
    if any(guide == seq for seq in guides):
        w = csv.writer(open(guide+"_aligns.csv", "wb"))
        w.writerow(['guide', 'target', 'alignment', 'count', 'inserts', 'deletes', 'mismatches', 'score','dstart','istart'])
        for target in gdict[guide]:
            if not(guide + "AGG" in target) and not(any(seq in target for seq in guides)):
                [s1a, s2a, adata, i, d, mm, score,dstart,istart] = affalign.alignmat(pre + guide + suff, target)
                if int(score) > 0:
                    w.writerow([s1a, s2a, adata, gdict[guide][target], i, d, mm, score,dstart,istart])