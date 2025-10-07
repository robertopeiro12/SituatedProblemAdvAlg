
amino = {'A' : ['GCT', 'GCC', 'GCA', 'GCG'], 
'R' : ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
'N' : ['AAT', 'AAC'],
'D' : ['GAT', 'GAC'],
'B' : ['AAT', 'AAC', 'GAT', 'GAC'],
'C' : ['TGT', 'TGC'],
'Q' : ['CAA', 'CAG'],
'E' : ['GAA', 'GAG'],
'Z' : ['CAA', 'CAG' 'GAA', 'GAG'],
'G' : ['GGT', 'GGC', 'GGA', 'GGG'],
'H' : ['CAT', 'CAC'],
'I' : ['ATT', 'ATC', 'ATA'],
'M' : ['ATG' ],
'L' : ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
'K' : ['AAA', 'AAG'],
'F' : ['TTT', 'TTC'],
'P' : ['CCT', 'CCC', 'CCA', 'CCG'],
'S' : ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
'T' : ['ACT', 'ACC', 'ACA', 'ACG'],
'W' : ['TGG'],
'Y' : ['TAT', 'TAC'],
'V' : ['GTT', 'GTC', 'GTA', 'GTG'],
'stop' : ['TAA', 'TGA', 'TAG']
}



f = open("SARS-COV-2-MN908947.3.txt", "r")
lines = f.readlines()
lines = lines[1:]
f.close()

print(lines)

virus  = ""
for l in lines:
	for c in l: 
		if c != "\n":
			virus = virus + c

virus = [virus[i : i + 3] for i in range(0, len(virus), 3)]

print(virus)
print(len(virus))

protein = "MGYINVFAFPFTIYSLLLCRMNSRNYIAQVDVVNFNLT"
orf10 = "ATGGGCTATATAAACGTTTTCGCTTTTCCGTTTACGATATATAGTCTACTCTTGTGCAGAATGAATTCTCGTAACTACATAGCACAAGTAGATGTAGTTAACTTTAATCTCACATAG"

orf10_codons = [ orf10[i:i+3] for i in range(0, len(orf10), 3) ]


for i in range(len(protein)):
	if orf10_codons[i] in amino[ protein[i] ]: 
		print("ok", protein[i], orf10_codons[i], amino[ protein[i] ] )
	else: 
		print("ERROR")
		exit(-1)

