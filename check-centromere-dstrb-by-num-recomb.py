# 200911 SN

# This will search and report the region with transmission rate same as that of centromere
# at specified numbers of recombination.

# You may run this for pedigreeSim out of FDR gamate to check SDR type centromere and vice versa.

# gamete-ploidy: integer
# type: integer (1) to check SDR (2) to check FDR
# num-recombination: integer

# python this.py input.pedigresim output gamete-ploidy type num-recombination

import procpedigreesim as prc
import sys
from decimal import *
from collections import defaultdict


def main():
    inp = open(sys.argv[1])
    o = open(sys.argv[2], 'w')
    ploidy = int(sys.argv[3])
    ctype = int(sys.argv[4])
    nrec = int(sys.argv[5])

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
    hit = 0
    o.write("gametes\tchromosome\thaplostruct\n")

    while True:
        gam = prc.create_gamete(inp, ploidy)
        if gam == "":
            break
        if gam == "out":
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
        if ctype == 1:
            res = r.check_sdr()
        elif ctype == 2:
            res = r.check_fdr()
        if 1 in res:
            if len(res) >= nrec + 1:
                if res[nrec] == 1:
                    c = 0
                    for i in gam:
                        o.write(str(hit)+"\t"+str(c)+"\t"+"\t".join([str(n) for n in i])+"\n")
                        c += 1
                    hit += 1

        count += 1

    inp.close()
    o.close()


main()
