# TD 1 Object Orientated Programming
# 08/12/2015

sequence = raw_input("Donner votre sequence")

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
			check +=
	if check = 0:
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
	for j in range(0,3):
	sequence2 = ""
		for i in range(0, len(sequence), 3):
			if sequence[i:i+3] != "taa" and sequence[i:i+3] != "tag" and sequence[i:i+3] !=  "tga":
				sequence2 += sequence[i:i+3]
			else:
			sequence2 += sequence[i:i+3].upper()
		print sequence2
		
# Detect START and STOP codons in all 6 possible phases

def brincodant(sequence):
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

def trouverstartstop(sequence,nvlseq):
	for j in range(0,3):
	sequence2 = ""
		for i in range(0, len(sequence), 3):
			if sequence[i:i+3] != "taa" and sequence[i:i+3] != "tag" and sequence[i:i+3] !=  "tga" and sequence[i:i+3] != "atg":
				sequence2 += sequence[i:i+3]
			else:
			sequence2 += sequence[i:i+3].upper()
		print "+", j,  sequence2
	for k in range(0,3):
	nvlseq2 = ""
		for i in range(0, len(nvlseq), 3):
			if nvlseq[i:i+3] != "tta" and nvlseq[i:i+3] != "cta" and nvlseq[i:i+3] !=  "tca" and nvlseq[i:i+3] != "cat":
				nvlseq2 += nvlseq[i:i+3]
			else:
			nvlseq2 += nvlseq[i:i+3].upper()
		print "-", k, nvlseq2
