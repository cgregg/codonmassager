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
import csv

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
    #initialize the function that generates all the codon usage tables in this environment
    def __init__(self, source_file):
        #print 'boom'

        # e.g.{'AAA': 'T', 'AAG': 'R'}
        self.codon_to_amino_acid_map = {}
        CUTinput = open(source_file)
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
        
        #iterate through codon_to_amino_acid_map to find each instance of an amino acid and add the corresponding codon to this new list.  inverse(codon_to_amino_acid_map)...
        # e.g.{'T': ['AAA', 'CCC', 'GGG'], 'A': ['CGC']}
        # self.amino_acid_to_codon_map = defaultdict(list) #if the value is going to show up mutliple time in TUCinput, then using defaultdict(set)
        # TUCinput = open(source_file)
        # for line in TUCinput:
        #     parts = line.split()
        #     a = 0
        #     while a < len(parts):
        #         codon = parts[a]
        #         AA = parts[a+1]
        #         self.amino_acid_to_codon_map[AA].append(codon) #if using defaultdict(set), then change append(codon) to add(codon)
        #         a = a+5

        # self.amino_acid_to_weight_map = defaultdict(list)
        # Freqinput = open(source_file)
        # for line in Freqinput:
        #     parts = line.split()
        #     a = 0
        #     while a < len(parts):
        #         AA = parts[a+1]
        #         Freq = parts[a+2]
        #         self.amino_acid_to_weight_map[AA].append(Freq)
        #         a = a+5
        
        #path to super map from RC
        self.amino_acid_to_all_map = defaultdict(dict)
        ALLinput = open(source_file)
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
        
        # #What follows is the codon table to direct rare codon usage in the beginning of proteins
        # self.amino_acid_to_all1st20_map = defaultdict(dict)
        # ALLinput = open(source_file)
        # for line in ALLinput:
        #     #line = line.rstrip('\r\n')
        #     lineParts = line.split()
        #     #print lineParts
        #     z = 0
        #     # need to do this b/c of error that happens when acccessing a key not there
        #     while z < len(lineParts):
        #         if lineParts[z+1] not in self.amino_acid_to_all1st20_map:
        #             self.amino_acid_to_all1st20_map[lineParts[z+1]] = {}
        #         self.amino_acid_to_all1st20_map[lineParts[z+1]][lineParts[z]] = float(lineParts[z+2])
        #         #print self.amino_acid_to_all1st20_map
        #         z = z + 5   
        # #print self.amino_acid_to_all1st20_map

        # self.amino_acid_to_all1st20_cumulativefreq = self.amino_acid_to_all1st20_map
        # for aminoAcid in self.amino_acid_to_all1st20_cumulativefreq: #nested for loop for dict in dict
        #     cumulativeFreq = 0
        #     for codon in self.amino_acid_to_all1st20_cumulativefreq[aminoAcid]:
        #         cumulativeFreq += self.amino_acid_to_all1st20_cumulativefreq[aminoAcid][codon]
        #         self.amino_acid_to_all1st20_cumulativefreq[aminoAcid][codon] = cumulativeFreq
        #print self.amino_acid_to_all1st20_cumulativefreq

        #This will turn AA into degenerate seqeunce for refactoring
        # self.amino_acid_to_degcodon_map = defaultdict(dict)
        # deginput = open('std-deg-usage.txt')
        # for line in deginput:
        #     #line = line.rstrip('\r\n')
        #     lineParts = line.split()
        #     #print lineParts
        #     z = 0
        #     # need to do this b/c of error that happens when acccessing a key not there
        #     while z < len(lineParts):
        #         if lineParts[z+1] not in self.amino_acid_to_degcodon_map:
        #             self.amino_acid_to_degcodon_map[lineParts[z+1]] = {}
        #         self.amino_acid_to_degcodon_map[lineParts[z+1]][lineParts[z]] = float(lineParts[z+2])
        #         #print self.amino_acid_to_degcodon_map
        #         z = z + 5   
        # print self.amino_acid_to_degcodon_map

        # self.amino_acid_to_degcodon_cumulativefreq = self.amino_acid_to_degcodon_map
        # for aminoAcid in self.amino_acid_to_degcodon_cumulativefreq: #nested for loop for dict in dict
        #     cumulativeFreq = 0
        #     for codon in self.amino_acid_to_degcodon_cumulativefreq[aminoAcid]:
        #         cumulativeFreq += self.amino_acid_to_degcodon_cumulativefreq[aminoAcid][codon]
        #         self.amino_acid_to_degcodon_cumulativefreq[aminoAcid][codon] = cumulativeFreq
        # print self.amino_acid_to_degcodon_cumulativefreq

    #define functions that the script can carry out as needed.  they need the codon/seq/AA/pseq argument passed to them.
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
    
    def reverse_translate_codon(self, AA):
        #print AA
        #need to iterate through input character by character
        revtrans = ""
        trp = self.amino_acid_to_all_cumulativefreq[AA]
        #print trp
        trpcodons = trp.keys() #Gleb, did I need to do it this way?  Is this slow?  I couldn't figure out how to get what I wanted out of the dictionary that I got out of the dictionary...  right.
        trpfreqs = trp.values()
        #print trpcodons#[0]
        #print trpfreqs
        choice = random.random()
        # print choice
        # print AA
        # print trpfreqs
        x = 0
        while choice >= trpfreqs[x]:
            x += 1
        else:
            trpchoice = trpcodons[x]
            #print trpchoice
            revtrans = revtrans + trpchoice
        return revtrans

    def reverse_translate_sequence(self, pseq):
        #setup input/output files, as args
        #output = open(args[1],'w')
        #input = open(args[0])
        #for line in input:
        #pseq = line.split()
        #print pseq
        #need to iterate through input character by character
        z = 0
        revtrans = ""
        while z < len(pseq):#20:
            if z < 20:
                self.table = self.amino_acid_to_all1st20_cumulativefreq
            else:
                self.table = self.amino_acid_to_all_cumulativefreq
            AA = pseq[z]
            trp = self.table[AA]
            #print trp
            trpcodons = trp.keys() #Gleb, did I need to do it this way?  Is this slow?  I couldn't figure out how to get what I wanted out of the dictionary that I got out of the dictionary...  right.
            trpfreqs = trp.values()
            #print trpcodons#[0]
            #print trpfreqs
            choice = random.random()
            #print choice
            x = 0
            while choice >= trpfreqs[x]:
                x = x + 1
            else:
                trpchoice = trpcodons[x]
                #print trpchoice
                revtrans = revtrans + trpchoice
            z = z + 1
            #print z
            #print revtrans
        
        # #print 'table switchover'
        # while z < len(pseq):
        #     AA = pseq[z]
        #     trp = self.amino_acid_to_all_cumulativefreq[AA]
        #     #print trp
        #     trpcodons = trp.keys() #Gleb, did I need to do it this way?  Is this slow?  I couldn't figure out how to get what I wanted out of the dictionary that I got out of the dictionary...  right.
        #     trpfreqs = trp.values()
        #     #print trpcodons#[0]
        #     #print trpfreqs
        #     choice = random.random()
        #     #print choice
        #     x = 0
        #     while choice >= trpfreqs[x]:
        #         x = x + 1
        #     else:
        #         trpchoice = trpcodons[x]
        #         #print trpchoice
        #         revtrans = revtrans + trpchoice
        #     z = z + 1
        #     #print z
        #     #print revtrans
        # print revtrans
        # revtrans = revtrans.replace('U','T')
        # print revtrans

        # output = open('revtransrprtr.txt','w')
        # output.write(pseq + '\t' + revtrans)
        # output.close()
        # for ID in d:
        #     count = d[ID] #create a variable (integer) that we can call in the write function below.
        #     output.write(ID + '\t' + str(count) + '\n') #write the ID then tab then the count converted to a string followed by a line break to output file.

#acutal nuts and bolts of script, exucuting the necessary functions in this environment based on the arguments passed at the command line.
def run_script(args):
    #open the input file, which should be args.scriptinput
    #send to the correct function based on the choice
    # print args.scriptinput
    scriptinput = open(args.scriptinput[0])
    outputfh = open(args.output[0],'w')
    output = csv.writer(outputfh)#output.txt','w')
    my_codon_table = CodonTable('rEcoli-codon-usage.txt')
    my_rare_codon_table = CodonTable('rEcoli-codon-usage-1st20.txt')
    scriptinputreader = csv.reader(scriptinput)
    if args.mode[0] == 'translate':
        for row in scriptinputreader:
            line = row[1]
            # print line
            codonseq = ''
            ti = 0
            line = line.replace('T','U')
            while ti + 2 < len(line):
                codon = line[ti] + line[ti+1] + line[ti+2]
                codonseq += my_codon_table.translate_codon(codon)
                ti += 3
            output.writerow([row[0]] + [codonseq]) 
    elif args.mode[0] == 'revtranslate':
        for row in scriptinputreader:
            line = row[1]
            # print line
            codonseq = ''
            ti = 0
            for character in line:
                if ti < 20:
                    table = my_rare_codon_table
                else:
                    table = my_codon_table
                codonseq += table.reverse_translate_codon(character)
                ti += 1
                #print ti
            codonseq = codonseq.replace('U','T')
            output.writerow([row[0]] + [codonseq])
    elif args.mode[0] == 'transcribe':
        for row in scriptinputreader:
            line = row[1]
            # print line
            codonseq = ''
            codonseq = line.replace('T','U')
            output.writerow([row[0]] + [codonseq])
            # else:
            #     print args.mode
    elif args.mode[0] == 'calCUT':
        calCUT = {}
        for row in scriptinputreader:
            line = row[1]
            ti = 0
            while ti + 2 < len(line):
                codon = line[ti] + line[ti+1] + line[ti+2]
                #print codon
                if not codon in calCUT:
                    calCUT[codon] = 0
                calCUT[codon] += 1
                ti += 3
                #print ti
        for codon in calCUT:
            count = calCUT[codon]
            cudun = codon.replace('T','U')
            AA = my_codon_table.translate_codon(cudun)
            output.writerow([codon] + [AA] + [count])
                
    # print codonseq
    # if str(codonseq) == str(TolC):
    #     print "nice"
    # else:
    #     print "you suck"
    outputfh.close()    

#output = open()
if __name__ == '__main__':
    print 'hello'
    import argparse
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('scriptinput', metavar='N',type=str, nargs=1, help='input file')
    parser.add_argument('mode', type=str, nargs=1, choices = ['translate', 'revtranslate', 'transcribe', 'calCUT'], help='what should I do?')
    parser.add_argument('output', type=str, nargs=1, help='output file')
    args = parser.parse_args()
    # print args#(args.accumulate(args.scriptinput))
    run_script(args)


# for line in input:
#     #print line #GK really likes print lines as a means to debug
#     #print line.split()
#     parts = line.split() #split input string whose fields are separated by white space into separate strings, demarcated by that white space.
#     ID = parts[1] #create variable and make it the second column ([1], a bucket)
#     if not ID in d: #if the new bucket is not in d, add it, add make its info=0.  d['key'] = x, makes x the string associated with key 'key'.
#         d[ID] = 0
#     d[ID] = d[ID] + 1        #counter outside of if statement
#     #sys.exit() #needs (), sys is function and functions expect (), then will look inside for arguments or nothing.
# #print d
# #sys.exit()

# for ID in d:
#     count = d[ID] #create a variable (integer) that we can call in the write function below.
#     output.write(ID + '\t' + str(count) + '\n') #write the ID then tab then the count converted to a string followed by a line break to output file.

# output.close() #always close output file.  function's require ().