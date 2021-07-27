# 200727 Soichiro Nishiyama
#
# This takes an output file of PedigreeSim with TEST mode as input
# and calculate ratio of allele transmittance and allele dosage by the number of recombination event.

# gamete-ploidy: integer
# bin-cM: decimal

# usage: python this.py input.pedigresim output gamete-ploidy

import procpedigreesim as prc
import sys
from decimal import *
from collections import defaultdict
import simprocess as smp


def output(res_list, count, o):
    c = max(count.keys()) + 1
    out = []
    gam = prc.gam_count(count)  # defaultdict cumulative sum of crossovers per gamete

    for i in range(5):
        out.append([])
        # trans, nul, sim, dup, tri

    for j in range(c):
        for i in range(5):
            a = smp.Stat(res_list[i][j])
            res = a.mean_max_min_count()
            out[i].append(res)

    o.write("num-recombination\tgam-count")
    nm1 = ["transmittance", "nuliplex", "simplex", "duplex", "triplex"]
    nm2 = ["-mean", "-min", "-mincount", "-max", "-maxcount"]
    for i in range(5):
        o.write("\t"+"\t".join([nm1[i]+n for n in nm2]))
    o.write("\n")

    for j in range(c):
        o.write(str(j)+"\t"+str(gam[j]))
        for i in range(5):
            o.write("\t" + "\t".join(map(str, out[i][j])))
        o.write("\n")


def main():
    inp = open(sys.argv[1])
    o = open(sys.argv[2], 'w')
    ploidy = int(sys.argv[3])

    while True:
        line = inp.readline()
        con = line.rstrip().split("\t")
        if con[0] == "chromosome":
            line = inp.readline()
            con = line.rstrip().split("\t")
            end = Decimal(con[3])
            line = inp.readline()
            break

    res_list = []
    for i in range(5):
        a = prc.reslist_prep(end, 1)
        res_list.append(a)

    # res lists of transmittance, nul, sim, dup, tri

    # res_dict_list = dict_prep(ploidy)

    count = 0
    recom_count = defaultdict(int)

    while True:
        gam = prc.create_gamete(inp, ploidy)
        if gam == "":
            output(res_list, recom_count, o)
            break
        if gam == "out":
            output(res_list, recom_count, o)

            res_list = []
            for i in range(5):
                a = prc.reslist_prep(end, 1)
                res_list.append(a)
            count = 0
            continue
        if count % 250 == 0:
            sys.stderr.write('Checking records: {}\r'.format(count*4))
        a = prc.Gametes(gam, end)
        rpos = a.recomb_pos()
        rmap = a.recomb_map()
        temp_dict_list = []
        for i in range(ploidy):
            temp_dict_list.append(defaultdict(int))
        d = prc.store_in_dict(temp_dict_list, rmap, ploidy, len(rpos), 1)
        r = prc.GamEval(d, ploidy, 1, len(rpos))
        res = r.check_inheritance()
        res_list = prc.mergeres(res_list, res)
        recom_count[len(rpos)] += 1
        count += 1

    inp.close()
    o.close()


main()
