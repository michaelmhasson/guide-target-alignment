import csv
from string import maketrans
import itertools
import affalign


def revcom(s):  # reverse complement function
    s = s[::-1]
    s = s.translate(maketrans("ATCG", "TAGC"))
    return s

gdict = {}
pre = 'GCTTTTTTTCTCGAGTACTA'
suff = 'AGGATCCATTAGGCGGCCGC'
guides = ['GGACCCTTCGTCGACACGAT', 'GGCATTGTCGATTCGCTCTC', 'GCGCGAAAATTCTATGAGAA', 'GATATCGGAGTTGTTGTTCG',
          'GGTTTGACCGTACTTGCAGC', 'GCGGAATCTCCTGTCCGAGA', 'GTCGACGATTTTCCACAACC', 'GGCGATCGACTTACGAATAC',
          'GTCCCTGCACGATACGCATA', 'GCGGCTACCTATTAACGTAG', 'GTCGGTTTATTTGACGAACG', 'GTCGGGTCATTCGTCGATTG',
          'GACCAGGATGGGCACCACCC']


def chooseguide(seq):
    score = '-100000'
    s = []
    for grna in guides:
        [s1a, s2a, adata, i, d, mm, tscore, dstart, istart] = affalign.alignmat(pre + grna + suff, seq)
        if int(tscore) > int(score):
            score = tscore
            s = [s1a, s2a, adata, i, d, mm, tscore, dstart, istart, grna]
    return s

with open('highO1.fasta') as f1, open('highO2.fasta') as f2:  # parse through both data files line by line
    for s1, s2 in itertools.izip_longest(f1, f2):
        if s1[0] != '>':
            s1 = s1.strip()
            s2 = s2.strip()
            if len(s1) == 20:  # 20 length sequence is the guide
                s2 = revcom(s2)
                if s1 in gdict:  # put in nested dicts with
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
w = csv.writer(open("all_nog.csv", "wb"))
w.writerow(['guide', 'guidealign', 'target', 'alignment', 'count', 'inserts', 'deletes', 'mismatches', 'score',
            'dstart', 'istart'])
unedit = 0
chimunedit = 0
chim = 0
nonchim = 0
for guide in gdict:
    for target in gdict[guide]:
        if not(guide + "AGG" in target):
            if not(any(seq + 'AGG' in target for seq in guides)):
                p = chooseguide(target)
                if p[9] == guide:
                    nonchim += 1
                else:
                    chim += 1
                if int(p[6]) > 0:
                    w.writerow([p[9], p[0], p[1], p[2], gdict[guide][target], p[3], p[4], p[5], p[6], p[7], p[8]])
            else:
                chimunedit += 1
        else:
            unedit += 1
