# 200910 SN

# This will search and report the region with transmission rate same as that of centromere.
# You may run this for pedigreeSim out of FDR gamate to check SDR type centromere and vice versa.

# gamete-ploidy: integer
# type: integer (1) to check SDR (2) to check FDR

# python this.py input.pedigresim output gamete-ploidy type

# update: add column for gametes without centromere genotype

import procpedigreesim as prc
import sys
from decimal import *
from collections import defaultdict


def output(out, count, o, ctype, errdic, centro):
    c = max(count.keys()) + 1
    gam = prc.gam_count(count)  # defaultdict cumulative sum of crossovers per gamete
    err = prc.gam_count(errdic)
    cen = prc.gam_count(centro)  # defaultdict cumulative sum of gametes without centromere genotype

    o.write("num-recombination\tgam-count\terror-gam")
    if ctype == 1:  # sdr
        o.write("\tSDRtype")
    if ctype == 2:
        o.write("\tFDRtype")
    o.write("\twithoutCentromere\n")

    for j in range(c):
        o.write(str(j)+"\t"+str(gam[j])+"\t"+str(err[j]))
        o.write("\t" + str(out[j][1])+"\t"+str(cen[j])+"\n")


def main():
    inp = open(sys.argv[1])
    o = open(sys.argv[2], 'w')
    ploidy = int(sys.argv[3])
    ctype = int(sys.argv[4])

    while True:
        line = inp.readline()
        con = line.rstrip().split("\t")
        if con[0] == "chromosome":
            line = inp.readline()
            con = line.rstrip().split("\t")
            end = Decimal(con[3])
            line = inp.readline()
            break

    count = 0
    recom_count = defaultdict(int)
    out = prc.dict_prep(int(end))
    errdic = defaultdict(int)
    valid_pos = defaultdict(int)  # for counting positions without centromere genotype

    while True:
        gam = prc.create_gamete(inp, ploidy)
        if gam == "":
            output(out, recom_count, o, ctype, errdic, valid_pos)
            break
        if gam == "out":
            output(out, recom_count, o, ctype, errdic, valid_pos)
            out = prc.dict_prep(int(end))
            count = 0
            continue
        if count % 250 == 0:
            sys.stderr.write('Checking records: {}\r'.format(count*4))
        a = prc.Gametes(gam, end)
        rpos = a.recomb_pos_eval()
        if len(rpos) >= 1:
            if rpos[0] == "err":
                errdic[len(rpos[1])] += 1
                continue
        rmap = a.recomb_map()
        temp_dict_list = []
        for i in range(ploidy):
            temp_dict_list.append(defaultdict(int))
        d = prc.store_in_dict(temp_dict_list, rmap, ploidy, len(rpos), 1)

        r = prc.GamEval(d, ploidy, 1, len(rpos))
        if ctype == 1:
            res = r.check_sdr()
        elif ctype == 2:
            res = r.check_fdr()
        for i in range(len(res)):
            out[i][res[i]] += 1

        if 1 not in res:
            valid_pos[len(rpos)] += 1
        else:
            valid_pos[res.index(1) - 1] += 1  # centromere genotype position closest to centromere
        recom_count[len(rpos)] += 1
        count += 1

    inp.close()
    o.close()


main()
