# TD 1 Object Orientated Programming
# 08/12/2015

class changeDNA(object):
	sequence = sequence.lower()
	def __init__(self, sequence):
		sequence.self = sequence
			
# Change DNA in RNA			
	def DNAinRNA(self, sequence):
		sequenceRNA = sequence.replace("t", "u")
	print sequenceRNA

# Detect METcodon in sequence
	def METsequence(self, sequence):
		if "atg" in range(len(sequence)):
			sequenceMET = sequence.replace("atg", "ATG")
	print sequenceMET
	
# Transcript DNA in cDNA
	def DNCincDNA(self, sequence): 
		nvlseq = ""
	for i in range(0, len(sequence)):
		if sequence[i] == "t":
			nvlseq += "a"
		if sequence[i] == "a":
			nvlseq += "t"
		if sequence[i] == "g":
			nvlseq += "c"
		if sequence[i] == "c":
			nvlseq += "g"
		return nvlseq

# Detect if a sequence is valid (only A, T, G, C)
	def validsequence (self, sequence):
		for i in sequence:
			check = 0
			if i != "a"  "t" or "g" and "c":
				check += 1
		if check == 0 :
			print "Sequence valid"
		else:
			print "Sequence not valid"

# Detect STOP codon in a sequence
	def STOPsequence(self, sequence):
		sequence2 = ""	
	for i in range(0, len(sequence), 3):
		if sequence[i:i+3] != "taa" and sequence[i:i+3] != "tag" and sequence[i:i+3] !=  "tga":
			sequence2 += sequence[i:i+3]
#		distinit = i
		else:
		#if distinit % 3 == 0:
			sequence2 += sequence[i:i+3].upper()
		print sequence2
		
# Detect STOP codon any sequence
	def STOPsequence(self, sequence):
		sequence2 = ""
        for j in range(0,3):
			for i in range(0, len(sequence), 3):
				if sequence[i:i+3] != "taa" and sequence[i:i+3] != "tag" and sequence[i:i+3] !=  "tga":
					sequence2 += sequence[i:i+3]
				else:
					sequence2 += sequence[i:i+3].upper()
			print sequence2
		
# Detect START and STOP codons in all 6 possible phases

	def brincodant(self, sequence):
		nvlseq = ""
		for i in range(0, len(sequence)):
			if sequence[i] == "t":
				nvlseq += "a"
			if sequence[i] == "a":
				nvlseq += "t"
			if sequence[i] == "g":
				nvlseq += "c"
			if sequence[i] == "c":
				nvlseq += "g"
		return nvlseq

	def trouverstartstop(self, sequence,nvlseq):
		sequence3 = ""
        nvlseq2 = ""

        for j in range(0,3):
			for i in range(0, len(sequence), 3):
				if sequence[i:i+3] != "taa" and sequence[i:i+3] != "tag" and sequence[i:i+3] !=  "tga" and sequence[i:i+3] != "atg":
					sequence3 += sequence[i:i+3]
				else:
					sequence3 += sequence[i:i+3].upper()
			print "+", j,  sequence3
		
        
        for k in range(0,3):
			for i in range(0, len(nvlseq), 3):
				if nvlseq[i:i+3] != "tta" and nvlseq[i:i+3] != "cta" and nvlseq[i:i+3] !=  "tca" and nvlseq[i:i+3] != "cat":
					nvlseq2 += nvlseq[i:i+3]
				else:
					nvlseq2 += nvlseq[i:i+3].upper()
			print "-", k, nvlseq2
            
            
masequence = changeDNA()

# Data structure

# Hydrophobic score of a peptide
class strdonnee (object):
	def __init__(self, var1):
		pass

	def hydrophobic(peptide):

	peptide = raw_input("Donner votre sequence peptidique")
	peptide = peptide.upper()

	listpep = {"D" : -1, "E" : -1, "K" : -1, "R" : -1, "N" : -0.5, "Q" : -0.5, "S" : -0.5, "I" : 1, "L" : 1, "M" : 1, "V" : 1, "F" : 1, "W" : 1, "C" : 0.5, "Y" : 0.5, "A" : 0, "G" : 0, "H" : 0, "P" : 0, "T" : 0}

	totalcount = 0
	for i in peptide:
		if i in listpep:
			totalcount += listpep[i]
	return totalcount
	print totalcount

# Hydrophobic score of a 5 peptids chain along a protein

peptide2 = raw_input("Donner votre sequence peptidique")
peptide2 = peptide2.upper()
listpep = {"D" : -1, "E" : -1, "K" : -1, "R" : -1, "N" : -0.5, "Q" : -0.5, "S" : -0.5, "I" : 1, "L" : 1, "M" : 1, "V" : 1, "F" : 1, "W" : 1, "C" : 0.5, "Y" : 0.5, "A" : 0, "G" : 0, "H" : 0, "P" : 0, "T" : 0}

signi = ""

def decalagecadre(peptide2)
	signi = ""
	for i in range(len(peptide2):
		peptide = sequence2[i:i+5]
		h = hydrophobic(peptide)
		if h > 1.5 :
			signi += "*"
		else:
			signi += " "
	print sequence
	print signi

# BLAST between 2 peptidic sequences (by match or family)

seq1 = ""
seq2 = ""

positivaminoacid = ["H", "K", "R"]
negativaminoacid = ["D", "E", "Q", "N"]
hydriphobicaminoacid = ["I", "L", "M", "V", "F", "W"]
smallsizeaminoacid = ["S" ,"T", "A" ,"G"]


resultat = ""
def compar(seq1, seq2)
	for i in seq1:
		for i in seq2:
			if i in seq1 == i in seq2:
				resultat += "A"
			elif:
				if i in seq1 = positiveaminoacid and i in seq2 = positiveaminoacid:
					resultat += "+"
				elif i in seq1 = negativaminoacid and i in seq2 = negativaminoacid:
					resultat += "+"			
				elif i in seq1 = hydriphobicaminoacid and i in seq2 = hydriphobicaminoacid:
					resultat += "+"
				elif i in seq1 = smallsizeaminoacid and i in seq2 = smallsizeaminoacid:
					resultat += "+"
			else:
				resultat += " "

	print seq1
	print resultat
	print seq2

# Traduction RNA in proteins

code = {'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R', 'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G'}

def traduction(sequence):
    proteine = ""
    codon = ""
	for i in range(0, len(sequence), 3):
        codon = sequence[i:i+3]
		if codon in code:
			proteine += code[codon]
	print proteine

traduction("AUGUGUAGGGCU")

# Find a concensus sequence among different sequences

def frequence(seq1, seq2, seq3, seq4, seq5):
	consensus = ""
	matrix = [[0 for x in range(len(seq1))] for x in range(4)]  # Or [[],[],[],[]]
	liste = [seq1, seq2, seq3, seq4, seq5]
	for i in range(len(seq1)):
		for j in range(len(liste)):
			if liste[j][i] == "A":
				matrix [0][i] += 1
			elif liste[j][i] == "T":
				matrix [1][i] += 1
			elif liste[j][i] == "G":
				matrix [2][i] += 1
			elif liste[j][i] == "C":
				matrix [3][i] += 1
	print matrix
	
	for i in range(len(seq1)):
		nuclMax=""
		scoreMax=0
		for k in range(0,4):
			if matrix[k][i] >= scoreMax:
				scoreMax = matrix[k][i]
				nuclMax = ("A","T","G","C")[k]
		consensus += nuclMax
	print "Sequence consencus:", consensus

seq1 = "GGTAGCT"
seq2 = "AACGATC"
seq3 = "AACGTTA"
seq4 = "AGCATCG"
seq5 = "ATAGCAA"
frequence(seq1, seq2, seq3, seq4, seq5)




### Read file

myfile = open("file.txt", "r")
contenu = monfichier.read()
monfichier.close()
