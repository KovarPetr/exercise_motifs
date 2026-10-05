def count_matrix(motifs):
    length = len(motifs[0])
    count_m = {"A" : [0]*length, "C" : [0]*length, "G" : [0]*length, "T" : [0]*length, }
    for motif in motifs:
        pos = 0
        for letter in motif:
            count_m[letter][pos] = count_m[letter][pos] + 1
            pos += 1
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

lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]
red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

print(score(red))                                # 26
print(consensus(best), score(best))              # AGATAG 34
print(hamming_distance("TAAGTT", "TGAATT"))      # 2
print(total_distance("TGCGTT", lecture_dna))     # 13