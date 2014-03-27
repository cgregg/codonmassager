from collections import defaultdict
#import the Sequence class definition     (note: the name of the class is "Seq") 
from Bio.Seq import Seq 
from Bio.Seq import translate #not Translate!!!!!! 
#import the IUPAC alphbet definitions 
from Bio.Alphabet import generic_dna, generic_rna, generic_protein, IUPAC 
from Bio.Data import CodonTable 

#create a Sequence object and assign it to TolC 
#positive reverse translation control and input for translation
tolC = Seq("ATGAAGAAATTGCTCCCCATTCTTATCGGCCTGAGCCTTTCTGGGTTCAGTTCGTTGAGC" + \
        "CAGGCCGAGAACCTGATGCAAGTTTATCAGCAAGCACGCCTTAGTAACCCGGAATTGCGT" + \
        "AAGTCTGCCGCCGATCGTGATGCTGCCTTTGAAAAAATTAATGAAGCGCGCAGTCCATTA" + \
        "CTGCCACAGCTAGGTTTAGGTGCAGATTACACCTATAGCAACGGCTACCGCGACGCGAAC" + \
        "GGCATCAACTCTAACGCGACCAGTGCGTCCTTGCAGTTAACTCAATCCATTTTTGATATG" + \
        "TCGAAATGGCGTGCGTTAACGCTGCAGGAAAAAGCAGCAGGGATTCAGGACGTCACGTAT" + \
        "CAGACCGATCAGCAAACCTTGATCCTCAACACCGCGACCGCTTATTTCAACGTGTTGAAT" + \
        "GCTATTGACGTTCTTTCCTATACACAGGCACAAAAAGAAGCGATCTACCGTCAATTAGAT" + \
        "CAAACCACCCAACGTTTTAACGTGGGCCTGGTAGCGATCACCGACGTGCAGAACGCCCGC" + \
        "GCACAGTACGATACCGTGCTGGCGAACGAAGTGACCGCACGTAATAACCTTGATAACGCG" + \
        "GTAGAGCAGCTGCGCCAGATCACCGGTAACTACTATCCGGAACTGGCTGCGCTGAATGTC" + \
        "GAAAACTTTAAAACCGACAAACCACAGCCGGTTAACGCGCTGCTGAAAGAAGCCGAAAAA" + \
        "CGCAACCTGTCGCTGTTACAGGCACGCTTGAGCCAGGACCTGGCGCGCGAGCAAATTCGC" + \
        "CAGGCGCAGGATGGTCACTTACCGACTCTGGATTTAACGGCTTCTACCGGGATTTCTGAC" + \
        "ACCTCTTATAGCGGTTCGAAAACCCGTGGTGCCGCTGGTACCCAGTATGACGATAGCAAT" + \
        "ATGGGCCAGAACAAAGTTGGCCTGAGCTTCTCGCTGCCGATTTATCAGGGCGGAATGGTT" + \
        "AACTCGCAGGTGAAACAGGCACAGTACAACTTTGTCGGTGCCAGCGAGCAACTGGAAAGT" + \
        "GCCCATCGTAGCGTCGTGCAGACCGTGCGTTCCTCCTTCAACAACATTAATGCATCTATC" + \
        "AGTAGCATTAACGCCTACAAACAAGCCGTAGTTTCCGCTCAAAGCTCATTAGACGCGATG" + \
        "GAAGCGGGCTACTCGGTCGGTACGCGTACCATTGTTGATGTGTTGGATGCGACCACCACG" + \
        "TTGTACAACGCCAAGCAAGAGCTGGCGAATGCGCGTTATAACTACCTGATTAATCAGCTG" + \
        "AATATTAAGTCAGCTCTGGGTACGTTGAACGAGCAGGATCTGCTGGCACTGAACAATGCG" + \
        "CTGAGCAAACCGGTTTCCACTAATCCGGAAAACGTTGCACCGCAAACGCCGGAACAGAAT" + \
        "GCTATTGCTGATGGTTATGCGCCTGATAGCCCGGCACCAGTCGTTCAGCAAACATCCGCA" + \
        "CGCACTACCACCAGTAACGGTCATAACCCTTTCCGTAACTGA", generic_dna)

#transcribe tolC, remember it is necessary for proper codon mapping
mRtolC = tolC.transcribe()

#positive control for translation and input for reverse translation
TolC = Seq("MKKLLPILIGLSLSGFSSLSQAENLMQVYQQARLSNPELRKSAADRDAAFEKINEARSPL" + \
        "LPQLGLGADYTYSNGYRDANGINSNATSASLQLTQSIFDMSKWRALTLQEKAAGIQDVTY" + \
        "QTDQQTLILNTATAYFNVLNAIDVLSYTQAQKEAIYRQLDQTTQRFNVGLVAITDVQNAR" + \
        "AQYDTVLANEVTARNNLDNAVEQLRQITGNYYPELAALNVENFKTDKPQPVNALLKEAEK" + \
        "RNLSLLQARLSQDLAREQIRQAQDGHLPTLDLTASTGISDTSYSGSKTRGAAGTQYDDSN" + \
        "MGQNKVGLSFSLPIYQGGMVNSQVKQAQYNFVGASEQLESAHRSVVQTVRSSFNNINASI" + \
        "SSINAYKQAVVSAQSSLDAMEAGYSVGTRTIVDVLDATTTLYNAKQELANARYNYLINQL" + \
        "NIKSALGTLNEQDLLALNNALSKPVSTNPENVAPQTPEQNAIADGYAPDSPAPVVQQTSA" + \
        "RTTTSNGHNPFRN*", IUPAC.protein) 
 
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

        #translating on my own, this works, but once I put it into the def translate_codon module, it doesn't
        tolCstr = str(mRtolC)
        trans = ""
        #print trans
        z = 0
        while z + 2 < len(tolCstr):
        #    print tolCstr[0]
            trp = tolCstr[z] + tolCstr[z+1] + tolCstr[z+2]
        #    print trp
            z = z + 3
            if trp in self.codon_to_amino_acid_map:
                #print 'In Here!'
                b = self.codon_to_amino_acid_map[trp]
                trans = trans + b
        print trans

        if str(trans) == str(TolC):
            print 'congrats, you reinvented the wheel'
        else:
            print 'you still suck'

        #print self.amino_acid_to_codon_map
        #if 'M' in self.amino_acid_to_codon_map.keys():
        #    print 'yep, M is in there'

        #reverse translate on my own
        TolCstr = str(TolC)
        rtrns = ""
        z = 0
        while z < len(TolCstr):
            b = TolCstr[z]
            #print b
            z = z + 1
            if b in self.amino_acid_to_codon_map:
                #print 'In Here!'
                cdn = self.amino_acid_to_codon_map[b]
                cdn = cdn[0]
                #print cdn
                rtrns = rtrns + cdn
        print rtrns
        rtrnstolC = Seq(rtrns, generic_rna)
        rtrnstolC = rtrnstolC.back_transcribe()
        print rtrnstolC
        #using biopython to translate
        rtrnsTolCtrns = rtrnstolC.translate(table=1, to_stop=False)
        #print tolCtrns
        if str(TolC) == str(rtrnsTolCtrns):
            print 'you reverse translated the wheel'
        else:
            print 'you still suck'
        print rtrnsTolCtrns

#to do: need to weight the codon picks.  also, it would be nice to re-factor, e.g., change while not minimizing RSCU.

#def translate_codon(self, codon):
    #tolCstr = str(mRtolC)
    #trans = ""
    #print trans
    #z = 0
    #while z + 2 < len(tolCstr):
    #    #print tolCstr[0]
    #    codon = tolCstr[z] + tolCstr[z+1] + tolCstr[z+2]
        #print codon
    #    z = z + 3
    #    if codon in self.codon_to_amino_acid_map:
    #        #print 'In Here!'
    #        b = self.codon_to_amino_acid_map[codon]
    #        trans = trans + b
    #print trans
    #if str(trans) == str(TolC):
    #    print 'congrats, you reinvented the wheel'
    #else:
    #    print 'you still suck'

#def reverse_translate_amino_acid(self, amino_acid):
#    codon_options = self.amino_acid_to_codon_map[amino_acid]

    # Somehow make a choice.  For now, just first one.
    #codon = codon_options[0]

    # Return the one we chose.
#    return codon


if __name__ == '__main__':
    my_codon_table = CodonTable('standard_usage.txt')
    #print my_codon_table.codon_to_amino_acid_map#['AAA']
    #print len(my_codon_table.codon_to_amino_acid_map)
    #print my_codon_table.amino_acid_to_codon_map
    #print len(my_codon_table.amino_acid_to_codon_map)
    #print my_codon_table.amino_acid_to_codon_map['A']
    #print my_codon_table.amino_acid_to_weight_map
    #print my_codon_table.translate_codon('AAA')
    #print my_codon_table.reverse_translate_amino_acid('F')

    #my_second_codon_table = CodonTable('weird_usage.txt')

