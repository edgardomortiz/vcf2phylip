#!/usr/bin/env python2
# -*- coding: utf-8 -*-



'''
The script converts a collection of SNPs in VCF format into a PHYLIP file for phylogenetic analysis
'''



__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "1.0"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2017-10-10"



import sys

filename = sys.argv[1]                    # First argument is the name of the VCF file
try: min_sample_locus = int(sys.argv[2])  # Second argument is minimum of samples per SNP, default is 4 for phylogenetics
except (ValueError, IndexError): min_sample_locus = 4

# Output filename will be the same as input file, indicating the minimum of samples specified
# e.g. VCF: sample.vcf, PHY: sample_min4.phy
outfile = filename.split(".")[0]+"_min"+str(min_sample_locus)+"."+filename.split(".")[-1].replace("vcf","phy").replace("VCF","phy")

# Dictionary of IUPAC ambiguities for nucleotides
amb = {("A","A"):"A",
	   ("A","C"):"M",
	   ("A","G"):"R",
	   ("A","N"):"A",
	   ("A","T"):"W",
	   ("C","A"):"M",
	   ("C","C"):"C",
	   ("C","G"):"S",
	   ("C","N"):"C",
	   ("C","T"):"Y",
	   ("G","A"):"R",
	   ("G","C"):"S",
	   ("G","G"):"G",
	   ("G","N"):"G",
	   ("G","T"):"K",
	   ("N","A"):"A",
	   ("N","C"):"C",
	   ("N","G"):"G",
	   ("N","N"):"N",
	   ("N","T"):"T",
	   ("T","A"):"W",
	   ("T","C"):"Y",
	   ("T","G"):"K",
	   ("T","N"):"T",
	   ("T","T"):"T"}

# Start processing of VCF file
with open(filename) as vcf:

	# Initialize line counter
	num_lines = 0

	# Dictionary that holds as keys the column index of the sample and as value a list
	# of two elements, element [0] is the name of the sample, element[1] is the DNA sequence
	seq = {}

	# Keep track of longest sequence name for padding with spaces in the output file
	longest_name = 0

	# Start parsing lines in VCF
	for line in vcf:

		# Split line into fields, 
		broken = line.strip("\n").split("\t")

		# Get sample names to create empty dictionary from the header of the file
		if broken[0] == "#CHROM":
			for i in range(9, len(broken)):
				seq[i] = [broken[i], ""]
				longest_name = max(longest_name, len(broken[i]))

		# Now process the SNPs one by one
		elif line[0] != "#":

			# Print progress every 10,000 lines
			num_lines += 1
			if num_lines % 10000 == 0:
				print str(num_lines)+" SNPs processed"

			# If SNP meets minimum of samples requirement
			if int(broken[7].split(";")[0].replace("NS=","")) >= min_sample_locus:

				# Create a dictionary for genotype to nucleotide translation
				# each SNP may code the nucleotides in a different manner
				nuc = {str(0):broken[3], ".":"N"}
				for n in range(0, len(broken[4].split(","))):
					nuc[str(n+1)] = broken[4].split(",")[n]

				# Translate genotype into nucleotides and the obtain the IUPAC ambiguity
				# for heterozygous SNPs, and append to DNA sequence of each sample
				for i in range(9, len(broken)):
					gt = broken[i].split(":")[0].replace("|","/")
					seq[i][1] += amb[(nuc[gt.split("/")[0]], nuc[gt.split("/")[1]])]
vcf.close()

# Write PHYLIP file
output = open(outfile, "w")

# PHYLIP header is number of samples and number of sites in alignment
phylip = str(len(seq))+" "+str(len(seq[9][1]))+"\n"
output.write(phylip)

# Initialize sample counter
sample_num = 0

# Loop the dictionary to obtain the DNA sequence while writing to output file
for sample in seq.keys():

	# Add padding with spaces after sequence name so every nucleotide starts
	# at the same position (only aestethic)
	padding = (longest_name + 3 - len(seq[sample][0])) * " "
	
	# Write sequence to file
	output.write(seq[sample][0]+padding+seq[sample][1]+"\n")
	
	# Print progress of writing the PHYLIP file
	sample_num += 1
	print "Sample "+str(sample_num)+" of "+str(len(seq))+" written to phylip file"
output.close()
print "Done!\n"

# END
