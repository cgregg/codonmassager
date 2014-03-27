from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
from Bio.Alphabet import generic_protein

my_dna = Seq("ATGGGGAGAAGGCCGTAG", generic_dna)
#print my_dna

#a = my_dna + 'aaa'
#print a

print my_dna.find('AGG')
print my_dna.find('AGA')
print my_dna
print my_dna.count('A')
print len(my_dna)

your_dna = my_dna.complement()
print your_dna
my_rna = my_dna.transcribe()
print my_rna

my_protr = my_rna.translate(table=1, to_stop=True) 
#table = 1 is default std genetic code, http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1
#to_stop=True tells it to stop at stops
print my_protr
my_protd = my_dna.translate(to_stop=True)
print my_protd

#playing with complete CDS'
#yaaX = Seq("GTGAAAAAGATGCAATCTATCGTACTCGCACTTTCCCTGGTTCTGGTCGCTCCCATGGCA" + \
#            "GCACAGGCTGCGGAAATTACGTTAGTCCCGTCAGTAAAATTACAGATAGGCGATCGTGAT" + \
#            "AATCGTGGCTATTACTGGGATGGAGGTCACTGGCGCGACCACGGCTGGTGGAAACAACAT" + \
#            "TATGAATGGCGAGGCAATCGCTGGCACCTACACGGACCGCCGCCACCGCCGCGCCACCAT" + \
#            "AAGAAAGCTCCTCATGATCATCACGGCGGTCATGGTCCAGGCAAACATCACCGCTAA",
#            generic_dna)

#print gene
#YAAX = yaaX.translate(table='Bacterial', cds=True, to_stop=True)
#print YAAX

#playing with codon usage tables
#from Bio.Data import CodonTable
#standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
#mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"]
#print standard_table

#mutable seq objects
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.Alphabet import IUPAC
#my_seq = Seq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)
#mutable_seq = my_seq.tomutable()
#Or just create a mutable seq!
my_seq = MutableSeq("GCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA", IUPAC.unambiguous_dna)
print my_seq
#my_seq_div = my_seq
#my_seq_div[5:8] = 'tag' #how to do insertions????????  only can replace as many characters as indicated.  wait it works now.  
#why 5:8?
#print my_seq #why does this print as my_seq_div with SNP?  
#print my_seq_div
#my_seq_del = my_seq_div.remove("T")
#print my_seq_del
my_seq_rev = my_seq.reverse() #should be able to do my_seq.reverse_complement() as well
print my_seq_rev #this should be working, but it returning None

fin_seq = my_seq_div.toseq() #converts back to immutable Seq Object