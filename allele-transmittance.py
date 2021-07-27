# 200727 Soichiro Nishiyama
#
# # This takes an output file of PedigreeSim with TEST mode as input
# # and calculate ratio of allele transmittance and allele dosage across the genetic map.

# This assumes that the TEST run of pedigreeSim was on a single chromosome.
# The average transmittance was calculated in each bin-cM.
# This assumes that the id of the founder allele is integer and found in HaploStruct columns.

# gamete-ploidy: integer
# bin-cM: decimal

# HaploStruct column:
# 0	0.6626	1	0.7001	0
# allele id / crossover position / allele id / ///

# The intervals will be calculated based on the distribution of the transmittance/doubling ratio
# (transmitted/doubled alleles over "ploidy" chromosomes) for each gamete.
# The mean +- 2 * SD was used for this estimation.

# usage: python this.py input.pedigresim output gamete-ploidy bin-cM \

import sys
from decimal import *
from collections import defaultdict
import simprocess
import procpedigreesim as prc


def output(res_list, end, bins, count, o):
    b = int(end / bins)
    out = []
    pos = []
    curr = 0

    for i in range(5):
        out.append([])
        # trans, nul, sim, dup, tri

    for j in range(b):
        for i in range(5):
            a = simprocess.Stat(res_list[i][j])
            res = a.mean_confidence_sd()
            out[i].append(res)
        pos.append(curr)
        curr = curr + bins
    for i in range(5):
        a = simprocess.Stat(res_list[i][j])
        res = a.mean_confidence_sd()
        out[i].append(res)
    pos.append(end)

    c = 0
    o.write(str(count) + " meiosis events generated\ndist")
    nm = ["transmittance", "nuliplex", "simplex", "duplex", "triplex"]
    for i in range(5):
        o.write("\t" + nm[i] + "-mean\t" + nm[i] + "-bottom\t" + nm[i] + "-top")
    o.write("\n")
    for j in range(len(out[0])):
        o.write(str(pos[c]))
        for i in range(5):
            o.write("\t" + "\t".join(map(str, out[i][j])))
        o.write("\n")
        c += 1


def main():
    inp = open(sys.argv[1])
    o = open(sys.argv[2], 'w')
    ploidy = int(sys.argv[3])
    bins = Decimal(sys.argv[4])

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
        a = prc.reslist_prep(end,bins)
        res_list.append(a)
    # res lists of transmittance, nul, sim, dup, tri

    # res_dict_list = dict_prep(ploidy)

    count = 0

    while True:
        gam = prc.create_gamete(inp, ploidy)
        if gam == "":
            output(res_list, end, bins, count, o)
            break
        if gam == "out":
            output(res_list, end, bins, count, o)

            res_list = []
            for i in range(5):
                a = prc.reslist_prep(end, bins)
                res_list.append(a)
            count = 0
            continue
        if count % 250 == 0:
            sys.stderr.write('Checking records: {}\r'.format(count*4))
        a = prc.Gametes(gam, end)
        a.bins = bins
        lines = a.create_phased_line()
        temp_dict_list = []
        for i in range(ploidy):
            temp_dict_list.append(defaultdict(int))
        d = prc.store_in_dict(temp_dict_list, lines, ploidy, end, bins)
        r = prc.GamEval(d, ploidy, bins, end)
        res = r.check_inheritance()
        res_list = prc.mergeres(res_list, res)
        count += 1

    inp.close()
    o.close()


main()
