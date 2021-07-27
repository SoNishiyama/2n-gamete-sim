# 201130 Soichiro Nishiyama
#
# This takes an output file of PedigreeSim with TEST mode as input
# and report average genetic distance by the number of recombination.

# gamete-ploidy: integer

# usage: python this.py input.pedigresim output gamete-ploidy

import procpedigreesim as prc
import sys
from decimal import *
from collections import defaultdict
import simprocess as smp


def output(res_list, o):
    o.write("recomb-count\tnum\tavr\t95%-percentile\t5%-percentile\n")
    for i in range(len(res_list)):
        a = smp.Stat(res_list[i])
        b = a.mean_percentile()
        o.write(str(i+1)+"\t"+str(len(res_list[i]))+"\t"
                +"\t".join([str(c) for c in b])+"\n")


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

    count = 0

    while True:
        gam = prc.create_gamete(inp, ploidy)
        if gam == "":
            output(res_list, o)
            break
        if gam == "out":
            output(res_list, o)

            res_list = []
            count = 0
            continue
        if count % 250 == 0:
            sys.stderr.write('Checking records: {}\r'.format(count*4))
        a = prc.Gametes(gam, end)
        rpos = a.recomb_pos()
        for i in range(len(rpos)):
            if i >= len(res_list):
                res_list.append([rpos[i]])
            else:
                res_list[i].append(rpos[i])
        count += 1

    inp.close()
    o.close()


main()