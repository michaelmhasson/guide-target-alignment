import csv
import os
from string import maketrans
import itertools
import affalign


def revcom(s):  # reverse complement function
    s = s[::-1]
    s = s.translate(maketrans("ATCG", "TAGC"))
    return s


def cutends(filenames):  # cutadapt to isolate guide and target sequences from raw data
    os.system('cutadapt -g ACATGCATGGCGGTAATACGGTTATCCACCCTGCAGG...ACCGACTCGGTGCCACTTTTTCAAGTTG '
              '-G ACATGCATGGCGGTAATACGGTTATCCACCCTGCAGG...ACCGACTCGGTGCCACTTTTTCAAGTTG '
              '-o tmpo1.fasta -p tmpo2.fasta ' + filenames)
    os.system('cutadapt -g TAAGTATCCCTTGGAGAACCACCTTGTTG...GTTTAAGAGCTAAGCTGGAAACAGCATAG '
              '-G TAAGTATCCCTTGGAGAACCACCTTGTTG...GTTTAAGAGCTAAGCTGGAAACAGCATAG '
              '-o full1.fasta -p full2.fasta tmpo1.fasta tmpo2.fasta')


def analyze(guides):  # performs analysis from given set of guides
    guidenums = {guides[i]: i+1 for i in range(len(guides))}
    gdict = {}
    with open('full1.fasta') as f1, open('full2.fasta') as f2:  # parse through both data files line by line
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
    pre = 'GCTTTTTTTCTCGAGTACTA'
    suff = 'AGGATCCATTAGGCGGCCGC'
    w = csv.writer(open("full_aligns.csv", "wb"))
    w.writerow(['guide', 'aligned guide', 'target', 'alignment', 'count', 'inserts', 'deletes', 'mismatches', 'score',
                'dstart', 'istart', 'insstr'])
    od = csv.writer(open("efficiency.csv", "wb"))
    od.writerow(['number','guide','count','edits','chimeric','efficiency'])
    otherguides = 0
    for guide in gdict:
        if any(guide == seq for seq in guides):
            nguide = 0
            nedit = 0
            nchim = 0
            for target in gdict[guide]:
                nguide += gdict[guide][target]
                if not(guide + "AGG" in target) and not(any(seq + "AGG" in target for seq in guides)):
                    nedit += gdict[guide][target]
                    [s1a, s2a, adata, i, d, mm, score, dstart, istart, insstr] = affalign.alignmat(pre + guide + suff, target)
                    if int(score) > 0:
                        w.writerow([guide, s1a, s2a, adata, gdict[guide][target], i, d, mm, score, dstart, istart, insstr])
                        if int(mm) > 3 or int(d) == 21 or int(d) == 22:
                            nchim += gdict[guide][target]
                    else:
                        nchim += gdict[guide][target]
            od.writerow([guidenums[guide], guide, str(nguide), str(nedit), str(nchim), str(float(nedit-nchim)/nguide)])
        else:
            for target in gdict[guide]:
                otherguides += gdict[guide][target]
    print(otherguides)


def fullrun(filenames,guides):
    cutends(filenames)
    analyze(guides)

# test run with our set of guides and the day 4 data
'''
default_guides = ['GGACCCTTCGTCGACACGAT', 'GGCATTGTCGATTCGCTCTC', 'GCGCGAAAATTCTATGAGAA', 'GATATCGGAGTTGTTGTTCG',
                  'GGTTTGACCGTACTTGCAGC', 'GCGGAATCTCCTGTCCGAGA', 'GTCGACGATTTTCCACAACC', 'GGCGATCGACTTACGAATAC',
                  'GTCCCTGCACGATACGCATA', 'GCGGCTACCTATTAACGTAG', 'GTCGGTTTATTTGACGAACG', 'GTCGGGTCATTCGTCGATTG',
                  'GACCAGGATGGGCACCACCC']

fullrun('newdata2_1.fastq newdata2_2.fastq',default_guides)
'''
