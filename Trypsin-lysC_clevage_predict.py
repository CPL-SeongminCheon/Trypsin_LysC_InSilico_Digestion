import sys,os,argparse
from Bio import SeqIO as seq
from itertools import combinations

def tryptic_cut(sequence,allow_miss_clevage,min_length):

	sequence = sequence.rstrip("*").upper()
	cut_locus=[0]

	for i in range(0, len(sequence)-1):    #Trypsin/lys-C clevage locus identify
		if sequence[i] == "K" or sequence[i] == "R":
			cut_locus.append(i+1)
	if cut_locus[-1] != len(sequence):      #if  not K in R in end of peptide, append end point in cut_locus 
		cut_locus.append(len(sequence))

	all_predict_peptide = []     # All predicted clevate peptide in list
	allow_tryptic_peptide = []   # only allow miss clevage and min_length peptide in list

	if len(cut_locus) > 2:     #Trypsin clevage site in peptide sequence

		for combi in combinations(cut_locus , 2): 	 #return combination by tuple in list
			all_predict_peptide.append(sequence[combi[0]:combi[1]])

		for line in all_predict_peptide:   # Check and Filtering by allow miss clevage sites
			if line[-1] =="K" or line[-1] =="R":
				line_checked = line[:-1]
			else:
				line_checked = line
			if (line_checked.count("R") + line_checked.count("K")) <= allow_miss_clevage:
				if len(line) >= min_length:
					allow_tryptic_peptide.append(line)
				else:
					continue

	else: 		#No trypsin clevage site in peptide sequence
		if len(sequence) >= min_length:
			allow_tryptic_peptide.append(sequence)

	return allow_tryptic_peptide


parser = argparse.ArgumentParser(description = "process trypsin/lys-c clevate site prediction from peptide sequence")
parser.add_argument("-i","--input",help="Input peptide fasta file name")
parser.add_argument("-o","--output",help="Output peptide fasta file name")
parser.add_argument("-miss","--missing",type=int,help="Allow miss clevage number  0,1,2 or more, default is 2", default=2)
parser.add_argument("-min", "--minimum_length", type=int, help="Minimum length of trypsin-lysC digested peptide length, default is 6", default=6)
args = parser.parse_args()

input = args.input
output = args.output
allow_miss_clevage = args.missing
min_length = args.minimum_length

if input is None or output is None:
	os.system("python Trypsin-lysC_clevage_predict.py -h") 
	sys.exit()

if input == output:
	print("Warning: input file name is same to output file name!!  See usage!")
	os.system("python Trypsin-lysC_clevage_predict.py -h")
	sys.exit()


out = open(output,"w")

for recod in seq.parse(input,"fasta"):
	seqid = str(recod.id.split(" ")[0])
	sequence = str(recod.seq)
	tryptic_peptide = tryptic_cut(sequence, allow_miss_clevage, min_length)

	flag = 0
	for line in tryptic_peptide:
		flag +=1
		
		if line[-1] =="K" or line[-1] =="R":
			check = line[:-1]
		else:
			check = line
		count_miss = check.count("R") + check.count("K")
		
		out.write(">" + seqid+"@TrypticLysC"+str(flag).zfill(4)+"MissClevage"+str(count_miss)+'\n'+line+'\n')

