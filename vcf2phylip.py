#!/usr/bin/env python2
# -*- coding: utf-8 -*-



'''
The script converts a collection of SNPs in VCF format into a PHYLIP, FASTA, 
NEXUS, or binary NEXUS file for phylogenetic analysis. The code is optimized
to process VCF files with sizes >1GB. For small VCF files the algorithm slows
down as the number of taxa increases (but is still fast).
'''


__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "1.5"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2018-04-24"


import sys
import os
import argparse


def main():
	parser = argparse.ArgumentParser(description="Converts SNPs in VCF format into an alignment for phylogenetic analysis")
	parser.add_argument("-i", "--input", action="store", dest="filename", required=True,
		help="Name of the input VCF file")
	parser.add_argument("-m", "--min-samples-locus", action="store", dest="min_samples_locus", type=int, default=4,
		help="Minimum of samples required to be present at a locus, default=4 since is the minimum for phylogenetics.")
	parser.add_argument("-o", "--outgroup", action="store", dest="outgroup", default="",
		help="Name of the outgroup in the matrix. Sequence will be written as first taxon in the alignment.")
	parser.add_argument("-p", "--phylip-disable", action="store_true", dest="phylipdisable", default=False,
		help="A PHYLIP matrix is written by default unless you enable this flag")
	parser.add_argument("-f", "--fasta", action="store_true", dest="fasta", default=False,
		help="Write a FASTA matrix, disabled by default")
	parser.add_argument("-n", "--nexus", action="store_true", dest="nexus", default=False,
		help="Write a NEXUS matrix, disabled by default")
	parser.add_argument("-b", "--nexus-binary", action="store_true", dest="nexusbin", default=False,
		help="Write a binary NEXUS matrix for analysis of biallelic SNPs in SNAPP, disabled by default")
	args = parser.parse_args()


	filename = args.filename
	min_samples_locus = args.min_samples_locus
	outgroup = args.outgroup
	phylipdisable = args.phylipdisable
	fasta = args.fasta
	nexus = args.nexus
	nexusbin = args.nexusbin
	ploidy = 2

	# Dictionary of IUPAC ambiguities for nucleotides
	# '*' means deletion for GATK (and other software?)
	# Deletions are ignored when making the consensus
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
		   ("T","T"):"T",
		   ("*","*"):"-",
		   ("A","*"):"A",
		   ("*","A"):"A",
		   ("C","*"):"C",
		   ("*","C"):"C",
		   ("G","*"):"G",
		   ("*","G"):"G",
		   ("T","*"):"T",
		   ("*","T"):"T",
		   ("N","*"):"N",
		   ("*","N"):"N"}


	# Dictionary for translating biallelic SNPs into SNAPP
	# 0 is homozygous reference
	# 1 is heterozygous
	# 2 is homozygous alternative
	gen_bin = {"./.":"?",
			   ".|.":"?",
			   "0/0":"0",
			   "0|0":"0",
			   "0/1":"1",
			   "0|1":"1",
			   "1/0":"1",
			   "1|0":"1",
			   "1/1":"2",
			   "1|1":"2"}


	# Process header of VCF file
	with open(filename) as vcf:

		# Create a list to store sample names
		sample_names = []

		# Keep track of longest sequence name for padding with spaces in the output file
		len_longest_name = 0

		# Look for the line in the VCF header with the sample names
		for line in vcf:
			if line.startswith("#CHROM"):

				# Split line into fields
				broken = line.strip("\n").split("\t")

				# If the minimum-samples-per-locus parameter is larger than the number of
				# species in the alignment make it the same as the number of species
				if min_samples_locus > len(broken[9:]):
					min_samples_locus = len(broken[9:])

				# Create a list of sample names and the keep track of the longest name length
				for i in range(9, len(broken)):
					name_sample = broken[i].replace("./","") # GATK adds "./" to sample names sometimes
					sample_names.append(name_sample)
					len_longest_name = max(len_longest_name, len(name_sample))

			# Find out the ploidy of the genotypes, just distinguishes if sample is not haploid vs n-ploid
			elif not line.startswith("#"):
				broken = line.strip("\n").split("\t")
				for j in range(9, len(broken[9:])):
					if broken[j].split(":")[0] not in [".",".|.","./."]:
						ploidy = len(broken[j].split(":")[0])
						# print broken[j]
					break
	vcf.close()


	# Output filename will be the same as input file, indicating the minimum of samples specified
	outfile = filename.replace(".vcf",".min"+str(min_samples_locus))

	# We need to create an intermediate file to hold the sequence data 
	# vertically and then transpose it to create the matrices
	if fasta or nexus or not phylipdisable:
		temporal = open(outfile+".tmp", "w")
	
	# if binary NEXUS is selected also create a separate temporal
	if nexusbin:
		temporalbin = open(outfile+".bin.tmp", "w")


	##################
	# PROCESS VCF FILE
	index_last_sample = len(sample_names)+9

	# Start processing SNPs of VCF file
	with open(filename) as vcf:

		# Initialize line counter
		snp_num = 0
		snp_accepted = 0
		snp_shallow = 0
		snp_multinuc = 0
		snp_biallelic = 0
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
					for g in range(9,len(broken)):
						if broken[g].split(":")[0] in [".", ".|."]:
							broken[g] = "./."

					# Keep track of number of genotypes processed
					snp_num += 1

					# Print progress every 500000 lines
					if snp_num % 500000 == 0:
						print str(snp_num)+" genotypes processed"

					# Check if the SNP has the minimum of samples required
					if (len(broken[9:]) - ''.join(broken[9:]).count("./.")) >= min_samples_locus:
						
						# Check that ref genotype is a single nucleotide and alternative genotypes are single nucleotides
						if len(broken[3]) == 1 and (len(broken[4])-broken[4].count(",")) == (broken[4].count(",")+1):

							# Add to running sum of accepted SNPs
							snp_accepted += 1

							# If nucleotide matrices are requested
							if fasta or nexus or not phylipdisable:

								# Create a dictionary for genotype to nucleotide translation
								# each SNP may code the nucleotides in a different manner
								nuc = {str(0):broken[3], ".":"N"}
								for n in range(len(broken[4].split(","))):
									nuc[str(n+1)] = broken[4].split(",")[n]

								# Translate genotypes into nucleotides and the obtain the IUPAC ambiguity
								# for heterozygous SNPs, and append to DNA sequence of each sample
								if ploidy > 1:
									site_tmp = ''.join([(amb[(nuc[broken[i][0]], nuc[broken[i][2]])]) for i in range(9, index_last_sample)])
								else:
									site_tmp = ''.join([(amb[(nuc[broken[i][0]], nuc[broken[i][0]])]) for i in range(9, index_last_sample)])

								# Write entire row of single nucleotide genotypes to temporary file
								temporal.write(site_tmp+"\n")

							# Write binary NEXUS for SNAPP if requested
							if nexusbin:

								# Check taht the SNP only has two alleles
								if len(broken[4]) == 1:
									
									# Add to running sum of biallelic SNPs
									snp_biallelic += 1

									# Translate genotype into 0 for homozygous ref, 1 for heterozygous, and 2 for homozygous alt
									binsite_tmp = ''.join([(gen_bin[broken[i][0:3]]) for i in range(9, index_last_sample)])

									# Write entire row to temporary file
									temporalbin.write(binsite_tmp+"\n")

						else:
							# Keep track of loci rejected due to multinucleotide genotypes
							snp_multinuc += 1
							# Keep track of loci rejected due to exceeded missing data
							snp_shallow += 1

					else:
						# Keep track of loci rejected due to exceeded missing data
						snp_shallow += 1

		# Print useful information about filtering of SNPs
		print str(snp_num) + " genotypes processed in total"
		print "\n"
		print str(snp_shallow) + " genotypes were excluded because they exceeded the amount of missing data allowed"
		print str(snp_multinuc) + " genotypes passed missing data filter but were excluded for not being SNPs"
		print str(snp_accepted) + " SNPs passed the filters"
		if nexusbin:
			print str(snp_biallelic) + " SNPs were biallelic and selected for binary NEXUS"
		print "\n"

	vcf.close()
	if fasta or nexus or not phylipdisable:
		temporal.close()
	if nexusbin:
		temporalbin.close()


	#######################
	# WRITE OUTPUT MATRICES

	if not phylipdisable:
		output_phy = open(outfile+".phy", "w")
		header_phy = str(len(sample_names))+" "+str(snp_accepted)+"\n"
		output_phy.write(header_phy)

	if fasta:
		output_fas = open(outfile+".fasta", "w")

	if nexus:
		output_nex = open(outfile+".nexus", "w")
		header_nex = "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=" + str(len(sample_names)) + " NCHAR=" + str(snp_accepted) + ";\n\tFORMAT DATATYPE=DNA" + " MISSING=N" + " GAP=- ;\nMATRIX\n"
		output_nex.write(header_nex)

	if nexusbin:
		output_nexbin = open(outfile+".bin.nexus", "w")
		header_nexbin = "#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX=" + str(len(sample_names)) + " NCHAR=" + str(snp_biallelic) + ";\n\tFORMAT DATATYPE=SNP" + " MISSING=?" + " GAP=- ;\nMATRIX\n"
		output_nexbin.write(header_nexbin)


	# Store index of outgroup in list of sample names
	idx_outgroup = "NA"

	# Write outgroup as first sequence in alignment if the name is specified
	if outgroup in sample_names:
		idx_outgroup = sample_names.index(outgroup)

		if fasta or nexus or not phylipdisable:
			with open(outfile+".tmp") as tmp_seq:
				seqout = ""

				# This is where the transposing happens
				for line in tmp_seq:
					seqout += line[idx_outgroup]

				# Write FASTA line
				if fasta:
					output_fas.write(">"+sample_names[idx_outgroup]+"\n"+seqout+"\n")
				
				# Pad sequences names and write PHYLIP or NEXUS lines
				padding = (len_longest_name + 3 - len(sample_names[idx_outgroup])) * " "
				if not phylipdisable:
					output_phy.write(sample_names[idx_outgroup]+padding+seqout+"\n")
				if nexus:
					output_nex.write(sample_names[idx_outgroup]+padding+seqout+"\n")

				# Print current progress
				print "Outgroup, "+outgroup+", added to the matrix(ces)."

		if nexusbin:
			with open(outfile+".bin.tmp") as bin_tmp_seq:
				seqout = ""

				# This is where the transposing happens
				for line in bin_tmp_seq:
					seqout += line[idx_outgroup]

				# Write line of binary SNPs to NEXUS
				padding = (len_longest_name + 3 - len(sample_names[idx_outgroup])) * " "
				output_nexbin.write(sample_names[idx_outgroup]+padding+seqout+"\n")

				# Print current progress
				print "Outgroup, "+outgroup+", added to the binary matrix."


	# Write sequences of the ingroup
	for s in range(0, len(sample_names)):
		if s != idx_outgroup:
			if fasta or nexus or not phylipdisable:
				with open(outfile+".tmp") as tmp_seq:
					seqout = ""

					# This is where the transposing happens
					for line in tmp_seq:
						seqout += line[s]

					# Write FASTA line
					if fasta:
						output_fas.write(">"+sample_names[s]+"\n"+seqout+"\n")
					
					# Pad sequences names and write PHYLIP or NEXUS lines
					padding = (len_longest_name + 3 - len(sample_names[s])) * " "
					if not phylipdisable:
						output_phy.write(sample_names[s]+padding+seqout+"\n")
					if nexus:
						output_nex.write(sample_names[s]+padding+seqout+"\n")

					# Print current progress
					print "Sample "+str(s+1)+" of "+str(len(sample_names))+", "+sample_names[s]+", added to the nucleotide matrix(ces)."

			if nexusbin:
				with open(outfile+".bin.tmp") as bin_tmp_seq:
					seqout = ""

					# This is where the transposing happens
					for line in bin_tmp_seq:
						seqout += line[s]

					# Write line of binary SNPs to NEXUS
					padding = (len_longest_name + 3 - len(sample_names[s])) * " "
					output_nexbin.write(sample_names[s]+padding+seqout+"\n")

					# Print current progress
					print "Sample "+str(s+1)+" of "+str(len(sample_names))+", "+sample_names[s]+", added to the binary matrix."

	if not phylipdisable:
		output_phy.close()
	if fasta:
		output_fas.close()
	if nexus:
		output_nex.write(";\nEND;\n")
		output_nex.close()
	if nexusbin:
		output_nexbin.write(";\nEND;\n")
		output_nexbin.close()

	if fasta or nexus or not phylipdisable:
		os.remove(outfile+".tmp")
	if nexusbin:
		os.remove(outfile+".bin.tmp")


	print "\nDone!\n"


if __name__ == "__main__":
    main()

