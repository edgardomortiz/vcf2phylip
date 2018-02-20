#!/usr/bin/env python2
# -*- coding: utf-8 -*-



'''
The script converts a collection of SNPs in VCF format into a NEXUS file for
phylogenetic analysis. The code is optmizied to process VCF files with sizes 
>1GB. For small VCF files the algorithm slows down as the number of taxa increases
(but is still fast). The NEXUS fiile can be also codified as binary characters
in the case of SNPs for analysis with SNAPP.
'''



__author__      = "Edgardo M. Ortiz"
__version__     = "1.0"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2018-02-14"


import sys
import os
import argparse


def main():
	parser = argparse.ArgumentParser(description="Converts SNPs in VCF format into a NEXUS matrix for phylogenetic analysis")
	parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
		help="Name of the input VCF file")
	parser.add_argument("-s", "--min-samples-locus", action="store", dest="min_samples_locus", type=int, default=4,
		help="Minimum of samples required to be present at a locus, default=4 since is the minimum for phylogenetics")
	parser.add_argument("-o", "--outgroup", action="store", dest="outgroup", default="",
		help="Name of the outgroup in the matrix. Sequence will be written as first taxon in the alignment")
	parser.add_argument("-b", "--binary", action="store_const", dest="nexusbin", const=1, default=0,
		help="Write binary format for SNPs analysis in SNAPP, default=nucleotide alignment if ommitted")
	args = parser.parse_args()

	filename = args.filename
	min_samples_locus = args.min_samples_locus
	outgroup = args.outgroup
	nexusbin = args.nexusbin

	# Output filename will be the same as input file, indicating the minimum of samples specified
	# e.g. VCF: sample.vcf, NEX: sample_min4.nex
	outfile = filename.replace(".vcf",".min"+str(min_samples_locus)+".nex")

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

	gen_bin = {"./.":"?",
			   ".|.":"?",
			   "0/0":"0",
			   "0|0":"0",
			   "0/1":"1",
			   "0|1":"1",
			   "1/0":"1",
			   "1|0":"1",
			   "1/1":"2",
			   "1/1":"2"}

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
				if not line.startswith("#") and line.strip("\n") != "": # pyrad sometimes produces an empty line after the #CHROM line

					# Split line into columns
					broken = line.strip("\n").split("\t")

					# Print progress every 100,000 lines
					snp_num += 1
					if snp_num % 100000 == 0:
						print str(snp_num)+" SNPs processed"

					# When nexusbin is 0 the program writes nucleotide alignment
					if nexusbin == 0: 

						# If SNP meets minimum of samples requirement
						if "NS=" in broken[7]: # for stacks, pyrad, ipyrad NS means Number of Samples in locus
							
							if int(broken[7].split("NS=")[1].split(";")[0]) >= min_samples_locus:

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

						# Optimize by reducing repetition of code when possible, otherwise seems to work fine
						elif "WT=" in broken[7] and "HET=" in broken[7] and "HOM=" in broken[7] and "NC=" in broken[7]:  # Regina's format from samtools
							if len(sample_names) - int(broken[7].split(";")[4].replace("NC=","")) >= min_samples_locus:

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

					# If nexusbin is 1 then write binary format for SNPs
					if nexusbin == 1: 

						# If SNP meets minimum of samples requirement
						if "NS=" in broken[7]: # for stacks, pyrad, ipyrad NS means Number of Samples in locus
							# If SNP meets the minimum samples requirement AND is biallelic
							if int(broken[7].split("NS=")[1].split(";")[0]) >= min_samples_locus and len(broken[4]) == 1:

								# Add to running sum of filtered SNPs
								snp_filtered += 1

								# Translate genotypes into nucleotides and the obtain the IUPAC ambiguity
								# for heterozygous SNPs, and append to DNA sequence of each sample
								site_tmp = ''.join([(gen_bin[broken[i][0:3]]) for i in range(9, index_last_sample)])

								# Write entire row of single nucleotide genotypes to temporary file
								temporal.write(site_tmp+"\n")

						# Optimize by reducing repetition of code when possible, otherwise seems to work fine
						elif "WT=" in broken[7] and "HET=" in broken[7] and "HOM=" in broken[7] and "NC=" in broken[7]:
							# If we use Regina's VCF from samtools  AND the SNP is biallelic:
							if len(sample_names) - int(broken[7].split(";")[4].replace("NC=","")) >= min_samples_locus and len(broken[4]) ==1:

								# Add to running sum of filtered SNPs
								snp_filtered += 1

								# Translate genotypes into nucleotides and the obtain the IUPAC ambiguity
								# for heterozygous SNPs, and append to DNA sequence of each sample
								site_tmp = ''.join([(gen_bin[broken[i][0:3]]) for i in range(9, index_last_sample)])

								# Write entire row of single nucleotide genotypes to temporary file
								temporal.write(site_tmp+"\n")

	vcf.close()
	temporal.close()



	# Write NEXUS file
	output = open(outfile, "w")

	# NEXUS header is number of samples and number of sites in alignment
	nexustype = "DNA"
	missing_data = "N"
	if nexusbin == 1:
		nexustype = "SNP"
		missing_data = "?"

	header = "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=" + str(len(sample_names)) + " NCHAR=" + str(snp_filtered) + ";\n\tFORMAT DATATYPE=" + nexustype + " MISSING=" + missing_data + " GAP=- ;\nMATRIX\n\n"
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
			print "Outgroup "+outgroup+" added to NEXUS file."

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
				print "Sample "+str(s+1)+" of "+str(len(sample_names))+" added to NEXUS file."

	output.write(";\nEND;\n")

	output.close()


	# Remove temporary file
	os.remove(outfile+".tmp")


	print "Done!\n"


if __name__ == "__main__":
    main()

