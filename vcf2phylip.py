#!/usr/bin/env python2
# -*- coding: utf-8 -*-



'''
The script converts a collection of SNPs in VCF format into a PHYLIP file for
phylogenetic analysis. The code is optmizied to process VCF files with sizes 
>1GB. For small VCF files the algorithm slows down as the number of taxa increases
(but is still fast).
'''



__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "1.0"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2017-10-10"



import sys

# First argument is the name of the VCF file
filename = sys.argv[1]

# Second argument is minimum of samples per SNP, default is 4 for phylogenetics
try: min_sample_locus = int(sys.argv[2])
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



# Process header of VCF file
with open(filename) as vcf:

	# List to store sample names
	sample_names = []

	# Keep track of longest sequence name for padding with spaces in the output file
	len_longest_name = 0

	# Look for the line in the VCF header with the sample names
	for line in vcf:
		if line.startswith("#CHROM"):

			# Split line into fields
			broken = line.strip("\n").split("\t")

			# Create a list of sample names and the length of the longest name
			for i in range(9, len(broken)):
				sample_names.append(broken[i])
				len_longest_name = max(len_longest_name, len(broken[i]))
			break
vcf.close()



# We need to create an intermediate file to hold the sequence data vertically
# and then transpose it to create the PHYLIP file
temporal = open(outfile+".tmp", "w")

# Store the index of the last column of VCF file
index_last_sample = len(sample_names)+9

# Start processing SNPs of VCF file
with open(filename) as vcf:

	# Initialize line counter
	snp_num = 0
	while True:

		# Load large chunks of file into memory
		vcf_chunk = vcf.readlines(100000)
		if not vcf_chunk:
			break

		# Now process the SNPs one by one
		for line in vcf_chunk:
			if not line.startswith("#"):

				# Split line into columns
				broken = line.strip("\n").split("\t")

				# Print progress every 100,000 lines
				snp_num += 1
				if snp_num % 100000 == 0:
					print str(snp_num)+" SNPs processed"

				# If SNP meets minimum of samples requirement
				if int(broken[7].split(";")[0].replace("NS=","")) >= min_sample_locus:

					# Create a dictionary for genotype to nucleotide translation
					# each SNP may code the nucleotides in a different manner
					nuc = {str(0):broken[3], ".":"N"}
					for n in range(len(broken[4].split(","))):
						nuc[str(n+1)] = broken[4].split(",")[n]

					# Translate genotypes into nucleotides and the obtain the IUPAC ambiguity
					# for heterozygous SNPs, and append to DNA sequence of each sample
					site_tmp = [(amb[(nuc[broken[i][0]], nuc[broken[i][2]])]) for i in range(9, index_last_sample)]

					# Write entire row of single nucleotide genotypes to temporary file
					temporal.write("\t".join(site_tmp)+"\n")
vcf.close()
temporal.close()



# Write PHYLIP file
output = open(outfile, "w")

# PHYLIP header is number of samples and number of sites in alignment
header = str(len(sample_names))+" "+str(snp_num)+"\n"
output.write(header)
for s in range(0, len(sample_names)):
	with open(outfile+".tmp") as tmp_seq:
		seqout = ""
		for line in tmp_seq:
			seqout += line.strip("\n").split("\t")[s]

		# Add padding with spaces after sequence name so every nucleotide starts
		# at the same position (only aestethic)
		padding = (len_longest_name + 3 - len(sample_names[s])) * " "

		# Write sample name and its corresponding sequence
		output.write(sample_names[s]+padding+seqout+"\n")

		# Print current progress
		print "Sample "+str(s+1)+" of "+str(len(sample_names))+" added to PHYLIP file."

output.close()

print "Done!\n"

# END
