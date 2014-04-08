from collections import defaultdict
#import the Sequence class definition     (note: the name of the class is "Seq") 
from Bio.Seq import Seq 
from Bio.Seq import translate #not Translate!!!!!! 
#import the IUPAC alphbet definitions 
from Bio.Alphabet import generic_dna, generic_protein, IUPAC 
from Bio.Data import CodonTable 
import sys 
import optparse
import random

parser = optparse.OptionParser() 
(options,args) = parser.parse_args()

# #setup input/output files, hardcoded
# #output = open('GAACCT_distro.txt','w')
# #input = open('GAACCT_output.txt')

#create a Sequence object and assign it to TolC 
#positive reverse translation control and input for translation
# tolC = Seq("ATGAAGAAATTGCTCCCCATTCTTATCGGCCTGAGCCTTTCTGGGTTCAGTTCGTTGAGC" + \
#         "CAGGCCGAGAACCTGATGCAAGTTTATCAGCAAGCACGCCTTAGTAACCCGGAATTGCGT" + \
#         "AAGTCTGCCGCCGATCGTGATGCTGCCTTTGAAAAAATTAATGAAGCGCGCAGTCCATTA" + \
#         "CTGCCACAGCTAGGTTTAGGTGCAGATTACACCTATAGCAACGGCTACCGCGACGCGAAC" + \
#         "GGCATCAACTCTAACGCGACCAGTGCGTCCTTGCAGTTAACTCAATCCATTTTTGATATG" + \
#         "TCGAAATGGCGTGCGTTAACGCTGCAGGAAAAAGCAGCAGGGATTCAGGACGTCACGTAT" + \
#         "CAGACCGATCAGCAAACCTTGATCCTCAACACCGCGACCGCTTATTTCAACGTGTTGAAT" + \
#         "GCTATTGACGTTCTTTCCTATACACAGGCACAAAAAGAAGCGATCTACCGTCAATTAGAT" + \
#         "CAAACCACCCAACGTTTTAACGTGGGCCTGGTAGCGATCACCGACGTGCAGAACGCCCGC" + \
#         "GCACAGTACGATACCGTGCTGGCGAACGAAGTGACCGCACGTAATAACCTTGATAACGCG" + \
#         "GTAGAGCAGCTGCGCCAGATCACCGGTAACTACTATCCGGAACTGGCTGCGCTGAATGTC" + \
#         "GAAAACTTTAAAACCGACAAACCACAGCCGGTTAACGCGCTGCTGAAAGAAGCCGAAAAA" + \
#         "CGCAACCTGTCGCTGTTACAGGCACGCTTGAGCCAGGACCTGGCGCGCGAGCAAATTCGC" + \
#         "CAGGCGCAGGATGGTCACTTACCGACTCTGGATTTAACGGCTTCTACCGGGATTTCTGAC" + \
#         "ACCTCTTATAGCGGTTCGAAAACCCGTGGTGCCGCTGGTACCCAGTATGACGATAGCAAT" + \
#         "ATGGGCCAGAACAAAGTTGGCCTGAGCTTCTCGCTGCCGATTTATCAGGGCGGAATGGTT" + \
#         "AACTCGCAGGTGAAACAGGCACAGTACAACTTTGTCGGTGCCAGCGAGCAACTGGAAAGT" + \
#         "GCCCATCGTAGCGTCGTGCAGACCGTGCGTTCCTCCTTCAACAACATTAATGCATCTATC" + \
#         "AGTAGCATTAACGCCTACAAACAAGCCGTAGTTTCCGCTCAAAGCTCATTAGACGCGATG" + \
#         "GAAGCGGGCTACTCGGTCGGTACGCGTACCATTGTTGATGTGTTGGATGCGACCACCACG" + \
#         "TTGTACAACGCCAAGCAAGAGCTGGCGAATGCGCGTTATAACTACCTGATTAATCAGCTG" + \
#         "AATATTAAGTCAGCTCTGGGTACGTTGAACGAGCAGGATCTGCTGGCACTGAACAATGCG" + \
#         "CTGAGCAAACCGGTTTCCACTAATCCGGAAAACGTTGCACCGCAAACGCCGGAACAGAAT" + \
#         "GCTATTGCTGATGGTTATGCGCCTGATAGCCCGGCACCAGTCGTTCAGCAAACATCCGCA" + \
#         "CGCACTACCACCAGTAACGGTCATAACCCTTTCCGTAACTGA", generic_dna)

# #transcribe tolC, remember it is necessary for proper codon mapping, but I fugred it out below...
# mRtolC = tolC.transcribe()

# #positive control for translation and input for reverse translation
# TolC = Seq("MKKLLPILIGLSLSGFSSLSQAENLMQVYQQARLSNPELRKSAADRDAAFEKINEARSPL" + \
#         "LPQLGLGADYTYSNGYRDANGINSNATSASLQLTQSIFDMSKWRALTLQEKAAGIQDVTY" + \
#         "QTDQQTLILNTATAYFNVLNAIDVLSYTQAQKEAIYRQLDQTTQRFNVGLVAITDVQNAR" + \
#         "AQYDTVLANEVTARNNLDNAVEQLRQITGNYYPELAALNVENFKTDKPQPVNALLKEAEK" + \
#         "RNLSLLQARLSQDLAREQIRQAQDGHLPTLDLTASTGISDTSYSGSKTRGAAGTQYDDSN" + \
#         "MGQNKVGLSFSLPIYQGGMVNSQVKQAQYNFVGASEQLESAHRSVVQTVRSSFNNINASI" + \
#         "SSINAYKQAVVSAQSSLDAMEAGYSVGTRTIVDVLDATTTLYNAKQELANARYNYLINQL" + \
#         "NIKSALGTLNEQDLLALNNALSKPVSTNPENVAPQTPEQNAIADGYAPDSPAPVVQQTSA" + \
#         "RTTTSNGHNPFRN*", IUPAC.protein) 
 
#create codon table dictionaries(3)
class CodonTable(object):
    """Defines how to make a codon table.
    """

    def __init__(self, source_file):
        #print 'boom'

        # e.g.{'AAA': 'T', 'AAG': 'R'}
        self.codon_to_amino_acid_map = {}
        CUTinput = open('ecoli-codon-usage.txt')
        for line in CUTinput:
        	parts = line.split()
        	#print parts
        	a = 0
           	while a < len(parts):
        		#print parts[a]
        		#print parts[a+1]
        		codon = parts[a]
        		AA = parts[a+1]
        		self.codon_to_amino_acid_map[codon] = AA
        		a = a+5
        
        # Key to list
        #iterate through codon_to_amino_acid_map to find each instance of an amino acid and add the corresponding codon to this new list.  inverse(codon_to_amino_acid_map)...
        # e.g.{'T': ['AAA', 'CCC', 'GGG'], 'A': ['CGC']}
        self.amino_acid_to_codon_map = defaultdict(list) #if the value is going to show up mutliple time in TUCinput, then using defaultdict(set)
        TUCinput = open('ecoli-codon-usage.txt')
        for line in TUCinput:
            parts = line.split()
            a = 0
            while a < len(parts):
                codon = parts[a]
                AA = parts[a+1]
                self.amino_acid_to_codon_map[AA].append(codon) #if using defaultdict(set), then change append(codon) to add(codon)
                a = a+5

        self.amino_acid_to_weight_map = defaultdict(list)
        Freqinput = open('ecoli-codon-usage.txt')
        for line in Freqinput:
            parts = line.split()
            a = 0
            while a < len(parts):
                AA = parts[a+1]
                Freq = parts[a+2]
                self.amino_acid_to_weight_map[AA].append(Freq)
                a = a+5
        
        #path to super map from RC
        self.amino_acid_to_all_map = defaultdict(dict)
        ALLinput = open('ecoli-codon-usage.txt')
        for line in ALLinput:
            #line = line.rstrip('\r\n')
            lineParts = line.split()
            #print lineParts
            z = 0
            # need to do this b/c of error that happens when acccessing a key not there
            while z < len(lineParts):
                if lineParts[z+1] not in self.amino_acid_to_all_map:
                    self.amino_acid_to_all_map[lineParts[z+1]] = {}
                self.amino_acid_to_all_map[lineParts[z+1]][lineParts[z]] = float(lineParts[z+2])
                #print self.amino_acid_to_all_map
                z = z + 5   
        #print self.amino_acid_to_all_map

        self.amino_acid_to_all_cumulativefreq = self.amino_acid_to_all_map
        for aminoAcid in self.amino_acid_to_all_cumulativefreq: #nested for loop for dict in dict
            cumulativeFreq = 0
            for codon in self.amino_acid_to_all_cumulativefreq[aminoAcid]:
                cumulativeFreq += self.amino_acid_to_all_cumulativefreq[aminoAcid][codon]
                self.amino_acid_to_all_cumulativefreq[aminoAcid][codon] = cumulativeFreq
        #print self.amino_acid_to_all_cumulativefreq

    def translate_codon(self, codon):
        AA = self.codon_to_amino_acid_map[codon]
        return AA

    def translate_sequence(self, seq):
        #easy transcriber!  don't need this biopython funciton now...  just call translate_sequence(seq)
        seqtrns = seq.replace('T','U')
        #print seqtrns

        z = 0
        trans = ""
        #print len(seqtrns)
        while z + 2 < len(seqtrns):
            #print z
            codon = seqtrns[z] + seqtrns[z+1] + seqtrns[z+2]
            z = z + 3
            b = self.codon_to_amino_acid_map[codon]
            trans = trans + b
        return trans #returns cause to break out of function currently in, regardles of loops
        
    def reverse_translate_sequence(self, pseq):
        #setup input/output files, as args
        #output = open(args[1],'w')
        #input = open(args[0])
        #for line in input:
        #pseq = line.split()
        print pseq
        #need to iterate through input character by character
        z = 0
        revtrans = ""
        while z < len(pseq):
            AA = pseq[z]
            trp = self.amino_acid_to_all_cumulativefreq[AA]
            #print trp
            trpcodons = trp.keys() #Gleb, did I need to do it this way?  Is this slow?  I couldn't figure out how to get what I wanted out of the dictionary that I got out of the dictionary...  right.
            trpfreqs = trp.values()
            #print trpcodons#[0]
            #print trpfreqs
            choice = random.random()
            print choice
            x = 0
            while choice >= trpfreqs[x]:
                x = x + 1
            else:
                trpchoice = trpcodons[x]
                print trpchoice
                revtrans = revtrans + trpchoice
            z = z + 1
            print z
            print revtrans
        revtrans = revtrans.replace('U','T')
        print revtrans

        output = open('revtransrprtr.txt','w')
        output.write(pseq + '\t' + revtrans)
        output.close()
        # for ID in d:
        #     count = d[ID] #create a variable (integer) that we can call in the write function below.
        #     output.write(ID + '\t' + str(count) + '\n') #write the ID then tab then the count converted to a string followed by a line break to output file.

#output = open()
if __name__ == '__main__':
    my_codon_table = CodonTable('standard_usage.txt')
    #print my_codon_table.translate_codon('AAA')
    print my_codon_table.translate_sequence('ATG')
    print my_codon_table.reverse_translate_sequence('MKKLLIGILSLSTYYFFWMDIPFFFFFFFF')#args[0]')
    #print len(my_codon_table.codon_to_amino_acid_map)
    #print my_codon_table.amino_acid_to_codon_map
    #print len(my_codon_table.amino_acid_to_codon_map)
    #print my_codon_table.amino_acid_to_codon_map['A']
    #print my_codon_table.amino_acid_to_weight_map
    #print my_codon_table.translate_codon('AAA')
    #print my_codon_table.reverse_translate_amino_acid('F')

    #my_second_codon_table = CodonTable('weird_usage.txt')

