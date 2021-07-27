import numpy as np
import scipy


class Stat:
    def __init__(self, data):
        self.data = data
        self.conf = 0.95

    def mean_confidence_interval(self):
        a = 1.0 * np.array(self.data)
        n = len(a)
        m, se = np.mean(a), scipy.stats.sem(a)
        h = se * scipy.stats.t.ppf((1 + self.conf) / 2., n - 1)
        return m, m - h, m + h

    def mean_confidence_sd(self):
        a = 1.0 * np.array(self.data)
        m, sd = np.mean(a), np.std(a)
        b = m - sd
        c = m + sd
        if b < 0:
            b = 0
        if c > 1:
            c = 1
        return m, b, c

    def mean_percentile(self):
        dat = np.array(self.data)
        top = np.percentile(dat, self.conf)
        bot = np.percentile(dat, 100 - self.conf)
        avr = np.mean(dat)
        return avr, bot, top

    def mean_max_min_count(self):
        dat = np.array(self.data)
        mx = np.max(dat)
        mx_count = np.count_nonzero(dat == mx)
        mn = np.min(dat)
        mn_count = np.count_nonzero(dat == mn)
        avr = np.mean(dat)
        return avr, mn, mn_count, mx, mx_count


class Transmittance:
    def __init__(self, check_dics, ploidy, bins, end):
        self.check_dics = check_dics
        self.ploidy = ploidy
        self.bins = bins
        self.end = end

    def transmit(self):
        b = int(self.end / self.bins)
        curr = 0
        out = []
        for j in range(b):
            temp = 0
            for i in range(self.ploidy):
                if self.check_dics[i][curr] >= 1:
                    temp += 1
            out.append(temp * 1.0 / self.ploidy)
            curr = curr + self.bins
        temp = 0
        for i in range(self.ploidy):
            if self.check_dics[i][self.end] >= 1:
                temp += 1
        out.append(temp * 1.0 / self.ploidy)
        return out

    def doubling(self):
        b = int(self.end / self.bins)
        curr = 0
        out = [[] for i in range(4)]  # nul, sim, dup, tri
        for j in range(b):
            nul = 0
            sim = 0
            dup = 0
            tri = 0
            for i in range(self.ploidy):
                if self.check_dics[i][curr] == 0:
                    nul += 1
                elif self.check_dics[i][curr] == 1:
                    sim += 1
                elif self.check_dics[i][curr] == 2:
                    dup += 1
                elif self.check_dics[i][curr] == 3:
                    tri += 1
            out[0].append(nul * 1.0 / self.ploidy)
            out[1].append(sim * 1.0 / self.ploidy)
            out[2].append(dup * 1.0 / self.ploidy)
            out[3].append(tri * 1.0 / self.ploidy)
            curr = curr + self.bins
        nul = 0
        sim = 0
        dup = 0
        tri = 0
        for i in range(self.ploidy):
            if self.check_dics[i][curr] == 0:
                nul += 1
            elif self.check_dics[i][curr] == 1:
                sim += 1
            elif self.check_dics[i][curr] == 2:
                dup += 1
            elif self.check_dics[i][curr] == 3:
                tri += 1
        out[0].append(nul * 1.0 / self.ploidy)
        out[1].append(sim * 1.0 / self.ploidy)
        out[2].append(dup * 1.0 / self.ploidy)
        out[3].append(tri * 1.0 / self.ploidy)
        return out

