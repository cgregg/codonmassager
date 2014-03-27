self.amino_acid_to_all_map = defaultdict(dict)
for line in ALLinput:
	line = line.rstrip('\r\n')
	lineParts = line.split('\t')
	# need to do this b/c of error that happens when acccessing a key not there
	if lineParts[1] not in self.amino_acid_to_all_map:
		self.amino_acid_to_all_map[lineParts[1]] = {}
		self.amino_acid_to_all_map[lineParts[1]][lineParts[0]] = float(lineParts[2])

for aminoAcid in self.amino_acid_to_all_map:
	cumulativeFreq = 0
for codon in self.amino_acid_to_all_map[aminoAcid]:
	cumulativeFreq += self.amino_acid_to_all_map[aminoAcid][codon]
	self.amino_acid_to_all_map[aminoAcid][codon] = cumulativeFreq