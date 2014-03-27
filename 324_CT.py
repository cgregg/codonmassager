from collections import defaultdict
#import the Sequence class definition     (note: the name of the class is "Seq") 
from Bio.Seq import Seq 
from Bio.Seq import translate #not Translate!!!!!! 
#import the IUPAC alphbet definitions 
from Bio.Alphabet import generic_dna, generic_protein, IUPAC 
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

#transcribe tolC, remember it is necessary for proper codon mapping, but I fugred it out below...
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

        self.amino_acid_to_all_map = defaultdict(dict)
        ALLinput = open('ecoli-codon-usage.txt')
        for line in ALLinput:
            line = line.rstrip('\r\n')
            lineParts = line.split('\t')
            print lineParts
            break
            # need to do this b/c of error that happens when acccessing a key not there
            if lineParts[1] not in self.amino_acid_to_all_map:
                self.amino_acid_to_all_map[lineParts[1]] = {}
                self.amino_acid_to_all_map[lineParts[1]][lineParts[0]] = float(lineParts[2])

        # for aminoAcid in self.amino_acid_to_all_map:
        #     cumulativeFreq = 0
        # for codon in self.amino_acid_to_all_map[aminoAcid]:
        #     cumulativeFreq +=` self.amino_acid_to_all_map[aminoAcid][codon]
        #     self.amino_acid_to_all_map[aminoAcid][codon] = cumulativeFreq
        # print self.amino_acid_to_all_map

#Next strategy, is to make a dictionary of lists
# and the 100 items are comprised of your codons
# based on the frequency
# and then you pull out a random number
# and take the mod of 100
# and use it as the index
# and the memory you use to store your table is negligible
# so i'd just need to .append values to the key, based on their freq
# correct...at the point where you have the frequency
# you just x by 100
# and then do a while loop counting down
# and at each step, append said codon to the list
# the only thing to worry about is when multiplying by 100, you need to cast the result into an int
# so you'll get the number as a string
# and in order to multiply it, you have to cast it to a float
# and then you have to take the resulting multiplication and cast THAT to an int
# so something like this:
# numVals = int(float(lineParts[2])*100)
# you might be better off using rounding as opposed to casting
# so rather than int, just use round
# b/c if you use int, 4.6 becomes 4
# as opposed to 5

where for each amino acid, you have a 100 items

    def translate_codon(self, codon):
        AA = self.codon_to_amino_acid_map[codon]
        return AA

    def translate_sequence(self, seq):
        #easy transcriber!  don't need this biopython funciton now..
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
        #this code is now useless, but was my first implementation of the translate algortihm using tolC as sepcific case example
        # tolCstr = str(mRtolC)
        # trans = ""
        # #print trans
        # z = 0
        # while z + 2 < len(tolCstr):
        #     #print tolCstr[0]
        #     codon = tolCstr[z] + tolCstr[z+1] + tolCstr[z+2]
        #     #print codon
        #     z = z + 3
        #     if codon in self.codon_to_amino_acid_map:
        #         #print 'In Here!'
        #         b = self.codon_to_amino_acid_map[codon]
        #         trans = trans + b
        # print trans
        # if str(trans) == str(TolC):
        #     print 'congrats, you reinvented the wheel'
        # else:
        #     print 'you still suck'

    def reverse_translate_amino_acid(self, amino_acid):
        codon_options = self.amino_acid_to_codon_map[amino_acid]

        # Somehow make a choice.  For now, just first one.
        codon = codon_options[0]

        # Return the one we chose.
        return codon


if __name__ == '__main__':
    my_codon_table = CodonTable('standard_usage.txt')
    print my_codon_table.translate_codon('AAA')
    print my_codon_table.translate_sequence('ATGAAGAAATTGCTCCCCATT')
    #print len(my_codon_table.codon_to_amino_acid_map)
    #print my_codon_table.amino_acid_to_codon_map
    #print len(my_codon_table.amino_acid_to_codon_map)
    #print my_codon_table.amino_acid_to_codon_map['A']
    #print my_codon_table.amino_acid_to_weight_map
    #print my_codon_table.translate_codon('AAA')
    #print my_codon_table.reverse_translate_amino_acid('F')

    #my_second_codon_table = CodonTable('weird_usage.txt')

