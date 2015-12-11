#!/usr/bin/python

# Python course

# -*- coding:utf-8 -*- 

import re

# Recherche d'une expression reguliere telephone
"""
chaine = ""
expression = r"^0[0-9]([ \.\-]?[0-9]{2}){4}$"
while re.search(expression, chaine) is None :
    chaine = input("Tapez votre numero de tel :")
"""

# Changement des expressions AT par TOTO dans une expression
"""
seq = "ATGGTAGATAG"

seq = re.sub("AT.G", "TOTO", seq)
print "Nouvelle sequence :" +seq
"""

# Demande de mot de passe a 6 characters
"""
chn_mdp = r"^[A-Za-z0-9]{6}$"
exp_mdp = re.compile(chn_mdp)
mot_passe = ""
while exp_mdp.search(mot_passe) is None:
	mot_passe = raw_input("Taper votre mot de passe")
"""

# Demande une sequence ATCGn ou atcgn
"""
sequence = raw_input("Taper votre sequence")
motif = r"^[ATCGNatgcn]+$"
motif_comp = re.compile(motif) # On compile le motif dans un objet
m = ""
while motif_comp.search(m) is None:
	m = raw_input("Sequence\n")
"""

# Cherche DEFINITION dans un fichier entier
'''
file  = open("Chr2.gb", "r")

def_motif = r"^DEFINITION\s+[^.]+" #On definit le motif
def_comp = re.compile(def_motif) # On le compile

for ligne in file : # On boucle sur toutes les lignes du fichier
	if (def_comp.search(ligne)): # Si on le trouve, on execute la ligne
		print ligne
		
file.close()
'''

'''
sequence =  r"(A[TG])"
retour = re.match(sequence, "ATGC")
print retour.group(1)


sequence =  r"(A[TG])"
retour = re.findall(sequence, "ATGC")
print retour.group(1)
'''

# Cherche ATGM12000 (ou autre) dans AT.txt puis l'imprime
'''
file = open("AT.txt", "r")
motif = r"AT(?P<NOMDUGROUPE>[\dCM])G([0-9]{5})"
motif_comp = re.compile(motif)

for ligne in file:
	if (motif_comp.search(ligne)):
		print ligne

for ligne in file:
	retour = motif_comp.match(ligne)
	if retour != None:
		print "Chromosome :" + retour.group(1)+ # .group Methode de RE qui va chercher ce qui est dans la n-ieme parenthese de notre motif recherche (0 = tout le motif, 1 = 1ere parenthese, 2 = 2eme parenthese etc... On peut aussi nommer le groupe avec ?P<NOMDUGROUPE> dans le motif, et on l'appelle avec retour.group("NOMDUGROUPE"))
		print "Gene N :" + retour groupe(2)

file.close()
'''

# Recuperer features
'''
file = open("/home/etudiant/Bureau/Chr2.gb", "r")
#motif_comp = re.compile(motif)

for ligne in file:
	retour = re.match(r"^\s{5}(\w+)\s+(\w+\()?(\d+)\.\.(\d+)", ligne)
	if retour != None:
		if retour.group(2) == "complement(":
			print "-",
		else:
			print"+",
		print retour.group(1), retour.group(3), retour.group(4)

file.close()
'''

# Read a fasta
'''
from Bio import SeqIO
for seq_record in SeqIO.parse("sequencecoli.fasta", "fasta"): # Fonction avec 2 parametres, nom et format du fichier. Elle lit le fichier

# Recuperation de l'ID de la sequence
	print seq_record.id
	print seq_record.name
	print seq_record.description

# Recuperation de la sequence
	print repr(seq_record.seq)
	print seq_record.seq
'''

# Read a multifasta
'''
from Bio import SeqIO

for seq_record in SeqIO.parse("multifasta.txt", "fasta"):
	print seq_record.id
	print seq_record.seq
	print len(seq_record.seq)
'''

# Convert a gb in fasta
'''
from Bio import SeqIO

infile = open ("Chr2.gb", "r"):
outfile = open ("Chr2.fasta", "w"):

seq_record = SeqIO.parse(infile, "gb")
SeqIO.write(seq_record, outfile, "fasta")
'''

# Create a record
'''
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

myseq = Seq("ATGGCGGTAG", generic_dna)
myseqrecord = SeqRecord(myseq, id = "MyID", description = "MyDesc")
print myseqrecord.format("gb")
'''
	
# Using a GenBank, write fasta with ID /leng/keyw/date/db_xref in first line
'''
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

for seq_record in SeqIO.parse("Chr2.gb", "gb"):
	print seq_record.id
	print seq_record.description
	print seq_record.annotations
	print seq_record.dbxrefs
	key = ""
	for k in seq_record.annotations["keywords"]:
		key = key + k + " "
	xref = ""
	for ref in seq_record.dbxrefs:
		xref = xref + ref + " "

	myrecord = SeqRecord(seq_record.seq, id = seq_record.id, description = seq_record.description+""+seq_record.annotations["date"]+""+xref+ " "+key)
	print myrecord.format("fasta")
'''

# Genome position using SeqFeature
"""
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import BeforePosition
from Bio.SeqFeature import SeqFeature

start = ExactPosition(250)
stop = BeforePosition(350)

loc = FeatureLocation(start, stop)
feat = SeqFeature(loc, strand = -1, type = "mRNA")
feat.qualifiers["test"] = "Blabla"

print start
print stop
print loc
print feat
"""

# Code to write a seqence with features
'''
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import ExactPosition
from Bio.SeqFeature import BeforePosition
from Bio.SeqFeature import SeqFeature
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


start = ExactPosition(250)
stop = BeforePosition(350)

loc = FeatureLocation(start, stop)
feat = SeqFeature(loc, strand = -1, type = "mRNA")
feat.qualifiers["test"] = "Blabla"

myfeats = []
myfeats.append(feat)

start = ExactPosition(2250)
stop = BeforePosition(2350)

loc = FeatureLocation(start, stop)
feat2 = SeqFeature(loc, strand = 1, type = "CDS")

myfeats.append(feat2)

myseq = Seq("ATCGAGCATCGATC", generic_dna)
myrecord = SeqRecord(myseq, id = "Yes", description = "cool", features = [feat])

print myrecord.format("genbank")
'''

# To send a blast on NCBI website: from Bio.Blast import NCBIWWW


# Les classes
'''
class Impr:
	varclasse = 12
	def imprimer(self, texte):
		print "je suis un langage :", texte
t = Impr()
t.imprimer("Oriente objet")
print t.varclasse
print Impr.varclasse
Impr.varclasse = 14 # On change la variable de classe (qui change la var pour toutes les instances de l'objet)
print t.varclasse
'''

# Les methodes de classes statique
'''
class Test:
	@staticmethod # classe statique : pas besoin d'instancier l'objet pour que la methode marche
	def mamethode(): # Pas de self car statique : si on veut utiliser des variables liees a l'instance, on a besoin du self
		print "ma methode"
Test.mamethode()
t = Test()
t.mamethode()
'''

# Methode __init__ : pour declarer des variables d'instance, le faire dans la methode init. Cree un objet vide
'''
class Test:
	def __init__(self):
		self.mavariable = 42

t = Test()
print t.mavariable
t.mavariable +=1
print t.mavariable

t2 = Test() # On cree une deuxieme instance mais avec la variable de classe et non la variable d'instance
print t2.mavariable
'''
# Methode __new__ : consteucteur d'objet, alloue un espace memoire pour l'objet et retourne un pointeur pour l'objet

# Encapsulation
'''
class Test :
	# Variable de classe
	def __init__(self):
		self.__mavariable = 12
	def affiche(self):
		self.__getv()
		return self.__mavariable
	def __getv(self):
		print "Ici"
	def __str__(self):

t = Test()
print t.affiche
print t.getv() # __getv est appelee via def affiche, imprime ici, puis revient dans def.affiche pour retourner mavariable
'''

'''
Test.__mavariable # J'appelle la classe et sa variable : marche pas. Normal, mavariable est une variable d'instance instanciee dans init

t = Test()
print t.__mavariable # Marche pas : mavariable est utilisable seulement dans l'instance. Pour l'utiliser :
'''

'''
class Test:
	def __init__(self):
		self.__mavariable = 12

	def affiche(self):
		self.__getv()
		return self.__mavariable

	def __getv(self):
		print "Ici"

	def __add__(self,v):
		return "Non je ne veux pas redefinir l'operateur addition"

	def __del__(self):
		print "I'm dead"

	def __repr__(self):
		return "Test!"

	def __str__(self):
		print "Je suis une classe de test"

t = Test()
print t.affiche
print repr(t)
print str(t)
del(t)
'''

# Heritage
'''
class Vehicule(object): # Object juste pour que la classe fille puis heriter de methodes
	def __init__(self, nb): # constructeur avec la variable
		print "je suis un vehicule"
		self.roues =nb
	def nbRoues(self): # Methode publique qui retourne le nb de roues
		return self.roues

class Voiture(Vehicule): # Classe specalise
	def __init__(self): #Cree l'objet voiture de la variable heritee
		super(Voiture, self).__init__(4) # Syntaxe pour appeler des methodes de l'objet dont herite (le parent)
		print "Je suis une une voiture"
		#self.roues = 4 Pas besoin, on a deja declare la variable dans init

	def __del__(self):
		print "Im dead"
		super(Voiture, self).__del__()

t = Voiture()
print t.nbRoues()
'''

# Ad-hoc polyphormism
'''
class Vehicule(object):
	def __init__(self, nb):
		print "je suis un vehicule"
		self.roues =nb
	def nbRoues(self):
		return self.roues

class Voiture(Vehicule):
	def __init__(self):
		super(Voiture, self).__init__(4)
		print "Je suis une une voiture"

	def __del__(self):
		print "Im dead"
		#super(Voiture, self).__del__()

	def __repr__(self): # Methode par defaut propre a python pour representer un objet
		return "Je suis une voiture"

a = 12
b = "toto"

t = Voiture()
print t.nbRoues()
print repr(t)
print repr(a) # Ca devrait representer un objet, mais comme on a redefinit repr, maintenant ca affiche a et b
print repr(b)
'''

# Polymorphisme d'heritage
'''
class Test():
	def __init__(self):
		self.mavariable = 12
	def getvariable(self):
		return self.__mavariable

class Test2(Test):
	def __init__(self):
		self.mavariable = 50
	def getvariable(self):
		print self.__mavariable + 12

t = Test()
print t.getvariable()

t2 = Test2()
print t2.getvariable()
'''

# Classe et methode abstraite
'''
class Vehicule(object):
	def __init__(self):
		print "je suis un vehicule"
		self.roues = 0
	def nbRoues(self):
		return self.roues
	def avance(self):
		raise NotImplementedError # On definit la classe virtuelle

class Velo(Vehicule):
		def avance (self): # Chaque methode est propre a chaque type de vehicule
			print "Je pedale fort"

class Voiture(Vehicule):
	def __init__(self):
		super(Voiture, self).__init__()
		print "Je suis une une voiture"
		self.roues = 4
	def __del__(self):
		print "Im dead"
	def avance(self): # Methode specifique a la voiture
		print "J'accelere"
	def __repr__(self):
		return "Je suis une voiture"

#v = Voiture()
#t = Velo()
#v.avance()
#t.avance()

t1 = Velo()
t2 = Voiture()
t3 = Velo()

print '----------'
tableau = [t1, t2, t3]

for v in tableau:
	v.avance() # Affiche le resultat de la fonction definie dans l'objet
'''

# Exceptions
'''
x = 0

def divise (x,y):
	assert y!=0, "Can't divide by 0"
	return x/y

print divise(10,5)

try:
	x = 10/int(raw_input("Input a number: "))
except ValueError: # Si pas un chiffre, afficher ce qu'il y a dans le value error
	print "Not a number !" 
except ZeroDivisionError:
	print "Can't divide by 0"
print x
'''

# Interface
'''class Test:
	def mamethode(self):
		raise NotImplementedError:
'''




