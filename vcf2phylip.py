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
__version__     = "1.2"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2017-10-10"


import sys
import os
import argparse


def main():
	parser = argparse.ArgumentParser(description="Converts SNPs in VCF format into a PHYLIP matrix for phylogenetic analysis")
	parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
		help="Name of the input VCF file")
	parser.add_argument("-s", "--min-samples-locus", action="store", dest="min_samples_locus", type=int, default=4,
		help="Minimum of samples required to be present at a locus, default=4 since is the minimum for phylogenetics")
	parser.add_argument("-o", "--outgroup", action="store", dest="outgroup", default="",
		help="Name of the outgroup in the matrix. Sequence will be written as first taxon in the alignment")
	args = parser.parse_args()

	filename = args.filename
	min_samples_locus = args.min_samples_locus
	outgroup = args.outgroup

	# Output filename will be the same as input file, indicating the minimum of samples specified
	# e.g. VCF: sample.vcf, PHY: sample_min4.phy
	outfile = filename+".min"+str(min_samples_locus)+".vcf"

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
		snp_filtered = 0
		while True:

			# Load large chunks of file into memory
			vcf_chunk = vcf.readlines(100000)
			if not vcf_chunk:
				break

			# Now process the SNPs one by one
			for line in vcf_chunk:
				if not line.startswith("#") or line == "":

					# Split line into columns
					broken = line.strip("\n").split("\t")

					# Print progress every 100,000 lines
					snp_num += 1
					if snp_num % 100000 == 0:
						print str(snp_num)+" SNPs processed"

					# If SNP meets minimum of samples requirement
					if int(broken[7].split(";")[0].replace("NS=","")) >= min_samples_locus:

						# Add to running sum of filtered SNPs
						snp_filtered += 1

						# Create a dictionary for genotype to nucleotide translation
						# each SNP may code the nucleotides in a different manner
						nuc = {str(0):broken[3], ".":"N"}
						for n in range(len(broken[4].split(","))):
							nuc[str(n+1)] = broken[4].split(",")[n]

						# Translate genotypes into nucleotides and the obtain the IUPAC ambiguity
						# for heterozygous SNPs, and append to DNA sequence of each sample
						site_tmp = ''.join([(amb[(nuc[broken[i][0]], nuc[broken[i][2]])]) for i in range(9, index_last_sample)])

						# Write entire row of single nucleotide genotypes to temporary file
						temporal.write(site_tmp+"\n")
	vcf.close()
	temporal.close()



	# Write PHYLIP file
	output = open(outfile, "w")

	# PHYLIP header is number of samples and number of sites in alignment
	header = str(len(sample_names))+" "+str(snp_filtered)+"\n"
	output.write(header)

	# Index of outgroup in list of sample names
	idx_outgroup = "NA"

	# Write outgroup as first sequence in alignment if the name is specified
	if outgroup in sample_names:
		idx_outgroup = sample_names.index(outgroup)
		with open(outfile+".tmp") as tmp_seq:
			seqout = ""
			for line in tmp_seq:
				seqout += line[idx_outgroup]

			# Add padding with spaces after sequence name so every nucleotide starts
			# at the same position (only aestethic)
			padding = (len_longest_name + 3 - len(sample_names[idx_outgroup])) * " "

			# Write sample name and its corresponding sequence
			output.write(sample_names[idx_outgroup]+padding+seqout+"\n")

			# Print current progress
			print "Outgroup "+outgroup+" added to PHYLIP file."

	# Write sequences of the ingroup
	for s in range(0, len(sample_names)):
		if s != idx_outgroup:
			with open(outfile+".tmp") as tmp_seq:
				seqout = ""
				for line in tmp_seq:
					seqout += line[s]

				# Add padding with spaces after sequence name so every nucleotide starts
				# at the same position (only aestethic)
				padding = (len_longest_name + 3 - len(sample_names[s])) * " "

				# Write sample name and its corresponding sequence
				output.write(sample_names[s]+padding+seqout+"\n")

				# Print current progress
				print "Sample "+str(s+1)+" of "+str(len(sample_names))+" added to PHYLIP file."

	output.close()


	# Remove temporary file
	os.remove(outfile+".tmp")


	print "Done!\n"


if __name__ == "__main__":
    main()

