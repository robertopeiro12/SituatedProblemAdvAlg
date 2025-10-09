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

# Step 1. Build amino dict
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

# Open and read the file
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

if key:
    proteins[key] = value

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


