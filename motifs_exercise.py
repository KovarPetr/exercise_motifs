from Bio import motifs
from Bio.Seq import Seq

bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
bio.pseudocounts = 1
print(bio.consensus)     # ATGCGTA
print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]

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

profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
print(profile.consensus())                       # ATGCGTA
print(round(profile.lmer_probability("ATGCGTA"), 4))  # 0.0122

two = MotifProfile(["GTAC", "TTAA"])
print(two.most_probable_lmer("ACTGGATGACCC"))    # TGAC
print(round(two.lmer_probability("TGAC"), 4))         # 0.0093