# 200727 Soichiro Nishiyama

# This takes an output file of PedigreeSim with TEST mode as input
# and calculate ratio of allele transmittance across the genetic map.

# This assumes that the TEST run was on a single chromosome.
# The average transmittance was calculated in each bin-cM.
# This assumes that the id of the founder allele is integer and found in HaploStruct columns.

# gamete-ploidy: integer
# bin-cM: decimal

# HaploStruct column:
# 0	0.6626	1	0.7001	0
# allele id / crossover position / allele id / ///

# usage: python this.py input.pedigresim output gamete-ploidy bin-cM

import sys
from decimal import *
from collections import defaultdict


def create_gamete(inp, ploidy):
    gamete = []
    count = 0
    while True:
        out = []
        line = inp.readline()
        if line[0:1] == "\t" or line[0:1] == "\n" or line == "":
            gamete = ""
            break
        con = line.rstrip().split("\t")
        if con[0] == "mei":
            continue
        count += 1
        c = 0
        for x in con[4:]:
            if c % 2 == 0:
                out.append(int(x))  # genotype
            else:
                out.append(float(x)*100)  # cM
            c += 1
        gamete.append(out)
        if count == ploidy:
            break
    return gamete


def create_phased_line(gametes, end, bins):
    lines = []
    b = int(end / bins)
    for chrom in gametes:
        out = []
        if len(chrom) == 1:
            for i in range(b):
                out.append(chrom[0])
            out.append(chrom[0])  # add at end
        else:
            curr = 0
            curr_idx = 1
            for i in range(b):
                if curr_idx == 0:
                    pass
                elif curr > chrom[curr_idx]:
                    curr_idx += 2
                    if curr_idx >= len(chrom):
                        curr_idx = 0
                out.append(chrom[curr_idx - 1])
                curr = curr + bins
            out.append(chrom[len(chrom) - 1])  # add at end
        lines.append(out)

    return lines  # [[1,1,2,2,2],[2,2,0,0,0]]


def store_in_dict(dic, lines, ploidy, end, bins):
    b = int(end / bins)
    for i in range(ploidy):
        curr = 0
        c = 0
        for j in lines[i]:
            if c == b:
                dic[j][end] += 1
                break
            dic[j][curr] += 1
            curr = curr + bins
            c += 1
    return dic


def check_doubling(check_dics,res_dics,ploidy,bins,end):
    b = int(end/bins)
    for i in range(ploidy):
        curr = 0
        for j in range(b):
            if check_dics[i][curr] == 2:
                res_dics[i][curr] += 1
            curr = curr + bins
        if check_dics[i][end] == 2:
            res_dics[i][end] += 1
    return res_dics


def main():
    inp = open(sys.argv[1])
    o = open(sys.argv[2], 'w')
    ploidy = int(sys.argv[3])
    bins = Decimal(sys.argv[4])

    st = 0

    while True:
        line = inp.readline()
        con = line.rstrip().split("\t")
        if con[0] == "chromosome":
            line = inp.readline()
            con = line.rstrip().split("\t")
            end = Decimal(con[3])
            break

    res_dict_list = []
    for i in range(ploidy):
        res_dict_list.append(defaultdict(int))

    count = 0

    while True:
        gam = create_gamete(inp, ploidy)
        if gam == "":
            break
        lines = create_phased_line(gam, end, bins)
        temp_dict_list = []
        for i in range(ploidy):
            temp_dict_list.append(defaultdict(int))
        d = store_in_dict(temp_dict_list, lines, ploidy, end, bins)
        res_dict_list = check_doubling(d, res_dict_list, ploidy, bins, end)
        count += 1

    inp.close()

    o.write("dist\tfreq\n")

    b = int(end / bins)
    out = []
    pos = []
    curr = 0

    for j in range(b):
        temp = []
        for i in range(ploidy):
            temp.append(str(Decimal(res_dict_list[i][curr])/count))
        out.append(temp)
        pos.append(curr)
        curr = curr + bins
    temp = []
    for i in range(ploidy):
        temp.append(str(Decimal(res_dict_list[i][end]) / count))
    out.append(temp)
    pos.append(end)

    c = 0
    for i in out:
        o.write(str(pos[c])+"\t"+"\t".join(i)+"\n")
        c += 1

    o.close()


main()
