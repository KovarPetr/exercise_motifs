import itertools
from Bio import SeqIO
import random
from Bio import motifs
from Bio.Seq import Seq

def count_matrix(motifs, prob=False, pseudocount=0):
    length = len(motifs[0])
    count_m = {"A" : [0]*length, "C" : [0]*length, "G" : [0]*length, "T" : [0]*length, }
    for motif in motifs:
        pos = 0
        for letter in motif:
            count_m[letter][pos] = count_m[letter][pos] + 1
            pos += 1
    if prob:
        for letter in count_m.keys():
            count_m[letter] = [(count_m[letter][pos] + pseudocount) / (len(motifs) + 4 * pseudocount)
                               for pos in range(len(count_m[letter]))]
    return count_m

def score(motifs):
    count_m = count_matrix(motifs)
    score = 0
    for position in range(len(motifs[0])):
        maxim = 0
        for letter in count_m.keys():
            if maxim < count_m[letter][position]:
                maxim = count_m[letter][position]
        score += maxim
    return score

def consensus(motifs):
    count_m = count_matrix(motifs)
    final_seq = ["A"] * len(motifs[0])
    for position in range(len(motifs[0])):
        maxim = 0
        for letter in count_m.keys():
            if maxim < count_m[letter][position]:
                maxim = count_m[letter][position]
                final_seq[position] = letter
    return ''.join(final_seq)

def hamming_distance(a, b):
    d = 0
    for position in range(len(a)):
        if a[position] != b[position]:
            d +=1
    return d

def total_distance(pattern, sequences):
    score = 0
    for sequence in sequences:
        scores = [hamming_distance(sequence[i:i+len(pattern)], pattern) for i in range(len(sequence)-len(pattern))]
        score += min(scores)
    return score

class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.l = len(motifs[0])
        self.ppm = count_matrix(motifs, True, pseudocount)

    def lmer_probability(self, lmer):
        prob = 1
        pos = 0
        for letter in lmer:
            prob = round(prob * self.ppm[letter][pos], 4)
            pos += 1
        return prob

    def most_probable_lmer(self, sequence):
        best_seq = ""
        score = 0
        for start in range(len(sequence)-self.l):
            seq = sequence[start:start + self.l]
            score_new = self.lmer_probability(seq)
            if score_new  > score:
                best_seq = seq
                score = score_new
        return best_seq

    def consensus(self):
        final_seq = ["A"] * self.l
        for position in range(self.l):
            maxim = 0
            for letter in self.ppm.keys():
                if maxim < self.ppm[letter][position]:
                    maxim = self.ppm[letter][position]
                    final_seq[position] = letter
        return ''.join(final_seq)

class MotifFinder:

    def __init__(self, sequences, l, seed=None):
        self.sequences = sequences
        self.l = l
        self.rng = random.Random(seed)
        self.windows = []
        for sequence in self.sequences:
            windows = []
            for i in range(len(sequence) - self.l):
                windows.append(sequence[i:i + self.l])
            self.windows.append(windows)

    def total_distance(self, pattern):
        return total_distance(pattern, self.sequences)

    def median_string(self):
        patterns = ["".join(lmer) for lmer in list(itertools.product("ACGT", repeat=self.l))]
        min_dist = 10000
        min_pattern = ""
        for pattern in patterns:
            dist = self.total_distance(pattern)
            if dist < min_dist:
                min_dist = dist
                min_pattern = pattern
        return min_pattern, min_dist

    def randomized_search(self):
        lmers = [self.rng.choice(seq) for seq in self.windows]
        getting_better = True
        while getting_better:
            profile = MotifProfile(lmers, pseudocount=1)
            new_lmers = []
            for seq in self.sequences:
                new_lmers.append(profile.most_probable_lmer(seq))
            score_ = score(lmers)
            new_score = score(new_lmers)
            if new_score > score_:
                lmers = new_lmers
            else:
                getting_better = False
        return lmers, score_

    def best_of(self, runs):
        best_motifs = []
        best_score = 0
        for i in range(runs):
            motifs, score = self.randomized_search()
            if score > best_score:
                best_motifs = motifs
                best_score = score
        return (best_motifs, best_score)

sequences = [str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")]
print(MotifFinder(sequences, 7).median_string())
best_motifs = MotifFinder(sequences, 7).best_of(100)
print(best_motifs)

out = []
for i in range(200):
    med_str = MotifFinder(sequences, 7).median_string()
    cons = consensus(MotifFinder(sequences, 7).randomized_search()[0])
    if med_str == cons:
        out.append(1)
    else:
        out.append(0)
p = sum(out)/200
print("Fraction of succesful: " + p)

out = [["R" , "measured success rate" , "predicted 1 − (1 − p) ^ R"]]
for R in [5, 10, 20, 50]:
    m_s_r = MotifFinder(sequences, 7).best_of(R)
    out.append([str(R), m_s_r, str(1-(1-p)^R)])
print(out)

found = motifs.create([Seq(m) for m in best_motifs])  # the motifs from best_of(100)
found.pseudocounts = 1
pssm = found.pssm
threshold = pssm.distribution().threshold_fpr(0.01)
for sequence in sequences:
    for position, hit_score in pssm.search(Seq(sequence), threshold=threshold):
        ...