from genbank.locus import Locus


for header,sequence in kseq(filename):
	yield Locus(header,sequence)

