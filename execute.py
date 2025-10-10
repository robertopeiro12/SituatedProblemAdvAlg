#Part 1 Find the 3 genes (M-gene, S-gene and ORF1AB-gene) in both sequence 
# algorithm: KMP 

def LPS(pattern):
    List = [0] * len(pattern) 
    L = 0  
    R = 1  
    while R < len(pattern):
        if pattern[L] == pattern[R]:
            List[R] = L + 1
            R += 1
            L += 1
        else:
            if L != 0:
                L = List[L - 1]
            else:
                R += 1

    return List

def kmp(text, pattern):
    i = 0  
    j = 0  
    lps = LPS(pattern)  

    while i < len(text):
        if text[i] == pattern[j]:
            i += 1
            j += 1
        else:
            if j > 0:
                j = lps[j - 1]
            else:
                i += 1

        if j == len(pattern):
            return i - j
        
def Z(text, pattern):
    n = len(text)
    m = len(pattern)
    x = n + m + 1

    # Concatenate pattern + "$" + text
    # "$" is a delimiter not present in either string
    C = pattern + "$" + text

    Z = [-1] * x  # Initialize Z array
    L = 0  # Left boundary of the current Z-box
    R = 0  # Right boundary of the current Z-box

    # Build the Z-array
    for i in range(x):
        if i <= R:
            k = i - L
            if Z[k] < R - i + 1:
                # If the previous Z-box can cover this position
                Z[i] = Z[k]
            else:
                # Otherwise, extend the Z-box further
                L = i
                while R < x and C[R - L] == C[R]:
                    R += 1
                Z[i] = R - L
                R -= 1
        else:
            # No Z-box covers i, start a new one
            L = R = i
            while R < x and C[R - L] == C[R]:
                R += 1
            Z[i] = R - L
            R -= 1

    return Z

sars = ["SARS-COV-2-MN908947.3.txt","SARS-COV-2-MT106054.1.txt"]
genes = ["M-gene.txt","S-gene.txt","ORF1AB-gene.txt"]
print("SARS-COV-2-MN908947.3 is Wuhan")
print("SARS-COV-2-MT106054.1 is Texas")
SARS = []
GENES = []

for i in sars:  
    f = open(i, "r")
    lines = f.readlines()
    lines = lines[1:]
    f.close()
    virus  = ""
    for l in lines:
        for c in l: 
            if c != "\n":
                virus = virus + c
    SARS.append(virus)

for i in genes:
    f = open(i, "r")
    lines = f.readlines()
    lines = lines[1:]
    f.close()
    gene  = ""
    for l in lines:
        for c in l: 
            if c != "\n":
                gene = gene + c
    GENES.append(gene)

i = 0
for virus in SARS:
    j=0
    print("Virus: "+sars[i])
    i+=1
    for gene in GENES:
        print("Gene: "+genes[j])
        numIndex = kmp(virus,gene)
        print("Index:", numIndex if numIndex is not None else "None")
        if numIndex is None:
            print("First 12 characters: None")
        else:
            print("First 12 characters:", virus[numIndex : numIndex + 12])
        j+=1
        

# Part 2: Find the longest palindrome in each gene sequence. (i.e., identify the index where it occurs)

def manacher(text):
    text = '#' + '#'.join(text) + '#'
    e = len(text)
    P = [0] * e
    center = limit = 0
    for i in range(1, e - 1):
        if i < limit:
            symmetrical = 2 * center - i
            P[i] = min(limit - i, P[symmetrical])
        gap = P[i] + 1
        while i - gap >= 0 and i + gap < e and text[i - gap] == text[i + gap]:
            P[i] += 1
            gap += 1
        if i + P[i] > limit:
            limit = i + P[i]
            center = i
    max_len = max(P)
    center_index = P.index(max_len)
    start = (center_index - max_len) // 2
    longest_palindrome = text[center_index - max_len:center_index + max_len + 1].replace('#', '')
    return longest_palindrome

print("\n" + "="*20 + "\n")
print("\nPart 2:")
print("\nLongest Palindrome in each gene:")
for i, gene in enumerate(GENES):
    print("Gene: " + genes[i])
    print(manacher(gene))

# Part 3: Find in which sections of the virus each protein is produced (i.e., identify the index of the genome where each of the 24 proteins occurs)
print("\n" + "="*20 + "\n")
print("Part 3:\nProteins in sequence:\n")

#  Build amino dict
amino = {}
with open("codes.txt", "r") as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        if ":" in line:
            key, values = line.split(":", 1)
            key = key.strip()
            codons = [v.strip() for v in values.replace(";", ",").split(",") if v.strip()]
            amino[key] = codons

SARS_codons = []

for virus in SARS:
    virus_codons = [virus[i:i+3] for i in range(0, len(virus), 3)]
    SARS_codons.append(virus_codons)

virustoproteinwuhan = ""
virustoproteintexas = ""
for value in SARS_codons[0]:
    for key, values in amino.items():
        if value in values:
            virustoproteinwuhan+=key
            break
for value in SARS_codons[1]:
    for key, values in amino.items():
        if value in values:
            virustoproteintexas+=key
            break

#  Build protein dict
with open("proteins-seq.txt", "r") as f:
    lines = f.readlines()

proteins = {}
key = None
value = ""

for line in lines:
    line = line.strip()
    if line.startswith(">"): 
        if key: 
            proteins[key] = value
        key = line[1:] 
        value = "" 
    else:
        value += line  

print("Wuhan:\n")
for key,value in proteins.items():
    index= kmp(virustoproteinwuhan,proteins[key][:4])
    if index == None:
        print(f"{key} not found in the genome")
    else:
        print(f"{key} : {index} - 4 amino: {value[:4]} = {SARS[0][index:index+12]}")


print("\nTexas:\n")

for key,value in proteins.items():
    index=kmp(virustoproteintexas,proteins[key][:4])
    if index == None:
        print(f"{key} not found in the genome")
    else:
        print(f"{key} : {index} - 4 amino: {value[:4]} = {SARS[1][index:index+12]}")


# Part 4:
'''
Compare the genomes: Wuhan, 2019, vs. Texas, 2020 Where are they the same? Different? Do the differences result in different amino acids? Specifically

Identify each index ranges where they differ and
for each such range in each genome, specify the codon- (3 letter sequences) and amino acid- (corresponding letter) sequences in those ranges
'''

print("\n" + "="*20 + "\n")
print("Part 4:\nDifferences between Wuhan and Texas:\n")


wuhan = SARS[0]
texas = SARS[1]

diffs = []
for i in range(min(len(wuhan), len(texas))):
    if wuhan[i] != texas[i]:
        diffs.append(i)

ranges = []
if diffs:
    start = diffs[0]
    prev = diffs[0]
    for i in diffs[1:]:
        if i == prev + 1:
            prev = i
        else:
            ranges.append((start, prev))
            start = i
            prev = i
    ranges.append((start, prev))


def translate_codon(codon):
    for aa, codons in amino.items():
        if codon in codons:
            return aa
    return "?"


for start, end in ranges:
    
    codon_index = (start // 3) * 3

    if codon_index + 3 <= len(wuhan) and codon_index + 3 <= len(texas):
        w_codon = wuhan[codon_index:codon_index + 3]
        t_codon = texas[codon_index:codon_index + 3]

        if w_codon != t_codon:
            w_aa = translate_codon(w_codon)
            t_aa = translate_codon(t_codon)
            print(f"Difference in index {codon_index}: Wuhan {w_codon} ({w_aa}) vs Texas {t_codon} ({t_aa})")
