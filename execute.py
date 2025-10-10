# Situated Problem in Advanced Algorithms
# Authors: Roberto Peiro & German Cueto




#Part 1 Find the 3 genes (M-gene, S-gene and ORF1AB-gene) in both sequence 
# algorithm: KMP 

def LPS(pattern):
    """
    Compute the Longest Proper Prefix Suffix (LPS) array for KMP algorithm
    This array is used to skip unnecessary comparisons during pattern matching
    """
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
    """
    Knuth-Morris-Pratt pattern matching algorithm
    Returns the index of first occurrence of pattern in text, or None if not found
    """
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
        

# File names for SARS virus sequences and gene sequences
sars = ["SARS-COV-2-MN908947.3.txt","SARS-COV-2-MT106054.1.txt"]
genes = ["M-gene.txt","S-gene.txt","ORF1AB-gene.txt"]
print("SARS-COV-2-MN908947.3 is Wuhan")
print("SARS-COV-2-MT106054.1 is Texas")

# Lists to store virus and gene sequences
SARS = []
GENES = []

# Read SARS virus sequences from files (skip header line)
for filename in sars:  
    f = open(filename, "r")
    lines = f.readlines()
    lines = lines[1:]
    f.close()
    virus  = ""
    for l in lines:
        for c in l: 
            if c != "\n":
                virus = virus + c
    SARS.append(virus)

# Read gene sequences from files (skip header line)
for filename in genes:
    f = open(filename, "r")
    lines = f.readlines()
    lines = lines[1:] 
    f.close()
    gene  = ""
    for l in lines:
        for c in l: 
            if c != "\n":
                gene = gene + c
    GENES.append(gene)

# Search for each gene in each virus sequence using KMP algorithm
virus_idx = 0
for virus in SARS:
    gene_idx = 0
    print("Virus: "+sars[virus_idx])
    virus_idx += 1
    for gene in GENES:
        print("Gene: "+genes[gene_idx])
        numIndex = kmp(virus,gene)
        print("Index:", numIndex if numIndex is not None else "None")
        if numIndex is None:
            print("First 12 characters: None")
        else:
            print("First 12 characters:", virus[numIndex : numIndex + 12])
        gene_idx += 1
        

# Part 2: Find the longest palindrome in each gene sequence. (i.e., identify the index where it occurs)

def manacher(text):
    """
    Manacher's algorithm to find the longest palindrome in a string
    Returns the longest palindromic substring
    """
    text = '#' + '#'.join(text) + '#'  # Insert separators between characters
    e = len(text)
    P = [0] * e  # Array to store palindrome lengths
    center = limit = 0
    for idx in range(1, e - 1):
        if idx < limit:
            symmetrical = 2 * center - idx
            P[idx] = min(limit - idx, P[symmetrical])
        gap = P[idx] + 1
        while idx - gap >= 0 and idx + gap < e and text[idx - gap] == text[idx + gap]:
            P[idx] += 1
            gap += 1
        if idx + P[idx] > limit:
            limit = idx + P[idx]
            center = idx
    max_len = max(P)
    center_index = P.index(max_len)
    start = (center_index - max_len) // 2
    longest_palindrome = text[center_index - max_len:center_index + max_len + 1].replace('#', '')
    return longest_palindrome

print("\n" + "="*20 + "\n")
print("\nPart 2:")
print("\nLongest Palindrome in each gene:")
# Find the longest palindrome in each gene sequence
for idx, gene in enumerate(GENES):
    print("Gene: " + genes[idx])
    print(manacher(gene))

# Part 3: Find in which sections of the virus each protein is produced (i.e., identify the index of the genome where each of the 24 proteins occurs)
print("\n" + "="*20 + "\n")
print("Part 3:\nProteins in sequence:\n")

# Build amino acid to codon mapping dictionary
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

# Convert virus sequences to codons (groups of 3 nucleotides)
SARS_codons = []

for virus in SARS:
    virus_codons = [virus[codon_idx:codon_idx+3] for codon_idx in range(0, len(virus), 3)]
    SARS_codons.append(virus_codons)

# Translate virus codons to amino acid sequences
virustoproteinwuhan = ""
virustoproteintexas = ""

# Translate Wuhan virus codons to amino acids
for value in SARS_codons[0]:
    for key, values in amino.items():
        if value in values:
            virustoproteinwuhan += key
            break

# Translate Texas virus codons to amino acids
for value in SARS_codons[1]:
    for key, values in amino.items():
        if value in values:
            virustoproteintexas += key
            break

# Build protein sequences dictionary from FASTA file
with open("proteins-seq.txt", "r") as f:
    lines = f.readlines()

proteins = {}
key = None
value = ""

# Parse FASTA format file
for line in lines:
    line = line.strip()
    if line.startswith(">"):  
        if key: 
            proteins[key] = value
        key = line[1:]  
        value = "" 
    else:
        value += line  

# Search for protein sequences in Wuhan virus
print("Wuhan:\n")
for key, value in proteins.items():
    # Search for first 4 amino acids of each protein in translated virus sequence
    index = kmp(virustoproteinwuhan, proteins[key][:4])
    if index == None:
        print(f"{key} not found in the genome")
    else:
        print(f"{key} : {index} - 4 amino: {value[:4]} = {SARS[0][index:index+12]}")

# Search for protein sequences in Texas virus
print("\nTexas:\n")
for key, value in proteins.items():
    index = kmp(virustoproteintexas, proteins[key][:4])
    if index == None:
        print(f"{key} not found in the genome")
    else:
        print(f"{key} : {index} - 4 amino: {value[:4]} = {SARS[1][index:index+12]}")


# Part 4: Compare the genomes: Wuhan, 2019, vs. Texas, 2020 
# Where are they the same? Different? Do the differences result in different amino acids?
# Identify each index ranges where they differ and
# for each such range in each genome, specify the codon- (3 letter sequences) 
# and amino acid- (corresponding letter) sequences in those ranges

print("\n" + "="*20 + "\n")
print("Part 4:\nDifferences between Wuhan and Texas:\n")

wuhan_seq = SARS[0]
texas_seq = SARS[1]

# Build codon-to-amino-acid dictionary for reverse lookup
CODON_TO_AA = {}
for aa, codons in amino.items():
    for codon in codons:
        CODON_TO_AA[codon] = aa

def find_codon_differences(seq1, seq2, label1="Wuhan", label2="Texas"):
    """
    Find differences between two sequences at the codon level
    Returns a list of differences with their positions and amino acid translations
    """
    diffs = []
    n_codons = min(len(seq1), len(seq2)) // 3
    
    for codon_idx in range(n_codons):
        codon1 = seq1[codon_idx*3:codon_idx*3+3]
        codon2 = seq2[codon_idx*3:codon_idx*3+3]
        if codon1 != codon2:
            aa1 = CODON_TO_AA.get(codon1, "?")
            aa2 = CODON_TO_AA.get(codon2, "?")
            diffs.append(
                f"Difference in index {codon_idx*3}: {label1} {codon1} ({aa1}) vs {label2} {codon2} ({aa2})"
            )
    return diffs

# Find and display all codon differences between Wuhan and Texas sequences
differences = find_codon_differences(wuhan_seq, texas_seq)
print("\n".join(differences))


