# 200908 Soichiro Nishiyama

# to process the output of pedigreeSim TEST mode

from collections import defaultdict
import simprocess


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
            gamete = "out"
            break
        count += 1
        c = 0
        for x in con[4:]:
            if c % 2 == 0:
                out.append(int(x))  # genotype
            else:
                out.append(float(x) * 100)  # cM
            c += 1
        gamete.append(out)
        if count == ploidy:
            break
    return gamete


class Gametes:
    def __init__(self, gam, end):
        self.gametes = gam
        self.end = end
        self.bins = 1

    def create_phased_line(self):
        lines = []
        b = int(self.end / self.bins)
        for chrom in self.gametes:
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
                    curr = curr + self.bins
                out.append(chrom[len(chrom) - 1])  # add at end
            lines.append(out)

        return lines  # [[1,1,2,2,2],[2,2,0,0,0]]

    def recomb_pos(self):
        le = []
        for i in self.gametes:
            le.append(int((len(i) - 1) / 2))
        pos = set()
        for i in range(len(self.gametes)):
            for j in range(le[i]):
                pos.add(self.gametes[i][j * 2 + 1])
        pos = sorted(list(pos))
        return pos

    def recomb_map(self):
        rec = self.recomb_pos()
        res = [[] for i in range(len(self.gametes))]
        curr = []
        for i in range(len(self.gametes)):
            curr.append(self.gametes[i][0])
            res[i].append(self.gametes[i][0])
        for j in range(len(rec)):
            for i in range(len(self.gametes)):
                if rec[j] in self.gametes[i][1::2]:
                    idx = self.gametes[i][1::2].index(rec[j])
                    change = self.gametes[i][0::2][idx + 1]
                    res[i].append(change)
                    curr[i] = change
                else:
                    res[i].append(curr[i])
        return res

    def recomb_pos_eval(self):
        le = []
        err = 0
        for i in self.gametes:
            le.append(int((len(i) - 1) / 2))
        pos = defaultdict(float)
        ori = defaultdict(lambda: "def")
        for i in range(len(self.gametes)):
            for j in range(le[i]):
                tar = self.gametes[i][j * 2 + 1]
                before = self.gametes[i][(j + 1) * 2]
                if ori[tar] != "def" and ori[tar] != before:
                    err = 1
                pos[tar] += 1
                ori[tar] = self.gametes[i][j * 2]

        e = list(pos.values())
        if len(e) >= 1:
            if max(e) >= 3 or err == 1:
                pos = ["err", sorted(list(pos.keys()))]
            else:
                pos = sorted(list(pos.keys()))
        else:
            pos = sorted(list(pos.keys()))
        return pos


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


def dict_prep(ploidy):
    dict_list = []
    for i in range(ploidy):
        dict_list.append(defaultdict(int))
    return dict_list


def reslist_prep(end, bins):
    a = [[] for i in range(int(end / bins) + 1)]
    return a


class GamEval:
    def __init__(self, dic, ploidy, bins, end):
        self.dic = dic
        self.ploidy = ploidy
        self.bins = bins
        self.end = end

    def check_inheritance(self):
        a = simprocess.Transmittance(self.dic, self.ploidy, self.bins, self.end)
        trans = a.transmit()
        b = simprocess.Transmittance(self.dic, self.ploidy, self.bins, self.end)
        double = b.doubling()
        res = [trans] + double
        return res

    def check_sdr(self):
        res = self.check_inheritance()  # res lists of transmittance, nul, sim, dup, tri
        out = []
        for i in range(len(res[0])):
            if res[1][i] == 0.5:  # nuliplex
                if res[2][i] == 0:  # simplex
                    out.append(1)
                    continue
            out.append(0)
        return out  # 1 = sdr, 0 = not sdr

    def check_fdr(self):
        res = self.check_inheritance()  # res lists of transmittance, nul, sim, dup, tri
        out = []
        for i in range(len(res[0])):
            if res[1][i] == 0:  # nuliplex
                if res[2][i] == 1:  # simplex
                    out.append(1)
                    continue
            out.append(0)
        return out  # 1 = fdr, 0 = not sdr


def gam_count(count_d):
    res = defaultdict(int)
    key = list(count_d.keys())
    for i in key:
        for j in range(i + 1):
            res[j] += count_d[i]
    return res


def mergeres(res_list, gameval):
    rc = len(gameval[0])
    for i in range(5):  # trans, nul, sim, dup, tri
        for j in range(rc):
            res_list[i][j].append(gameval[i][j])
    return res_list


def cumulative_gametes(defauldic):
    res = defaultdict(int)
    key = list(defauldic.keys())
    for i in key:
        for j in range(i):
            res[j] += 1
    return res
