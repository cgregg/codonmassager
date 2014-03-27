# Example 3-1 http://www1.chapman.edu/~fahy/bioinformatics/chapter3.html

#import the Sequence class definition     (note: the name of the class is "Seq")
from Bio.Seq import Seq
from Bio.Seq import translate #not Translate!!!!!!
#import the IUPAC alphbet definitions
from Bio.Alphabet import generic_dna, generic_protein, IUPAC
from Bio.Data import CodonTable

input = open('')

print CodonTable.standard_dna_table
sys.exit()
my_CodonTable = CodonTable
#def back_translate (prot, CUT)


#Fetch translators with codon usage determined (NCBI)
#standard_translator = translate.unambiguous_dna_by_id[1]
#mito_translator = translate.unambiguous_dna_by_id[2]

#create a Sequence object and assign it to TolC
TolC = Seq("MKKLLPILIGLSLSGFSSLSQAENLMQVYQQARLSNPELRKSAADRDAAFEKINEARSPL" + \
		"LPQLGLGADYTYSNGYRDANGINSNATSASLQLTQSIFDMSKWRALTLQEKAAGIQDVTY" + \
		"QTDQQTLILNTATAYFNVLNAIDVLSYTQAQKEAIYRQLDQTTQRFNVGLVAITDVQNAR" + \
		"AQYDTVLANEVTARNNLDNAVEQLRQITGNYYPELAALNVENFKTDKPQPVNALLKEAEK" + \
		"RNLSLLQARLSQDLAREQIRQAQDGHLPTLDLTASTGISDTSYSGSKTRGAAGTQYDDSN" + \
		"MGQNKVGLSFSLPIYQGGMVNSQVKQAQYNFVGASEQLESAHRSVVQTVRSSFNNINASI" + \
		"SSINAYKQAVVSAQSSLDAMEAGYSVGTRTIVDVLDATTTLYNAKQELANARYNYLINQL" + \
		"NIKSALGTLNEQDLLALNNALSKPVSTNPENVAPQTPEQNAIADGYAPDSPAPVVQQTSA" + \
		"RTTTSNGHNPFRN*", IUPAC.protein)

#tolC = TolC.back_translate[1]#(table='Standard')
#tolC = back_translate(TolC, table='Standard', stop_symbol='*', to_stop=True, cds=True)
#tolC = standard_translator.back_translate(TolC)
print tolC

output = open('20140318revtrans.txt','w')
output.write(tolC + '\n')
output.close()
