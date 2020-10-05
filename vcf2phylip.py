#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""
The script converts a collection of SNPs in VCF format into a PHYLIP, FASTA, 
NEXUS, or binary NEXUS file for phylogenetic analysis. The code is optimized
to process VCF files with sizes >1GB. For small VCF files the algorithm slows
down as the number of taxa increases (but is still fast).

Any ploidy is allowed, but binary NEXUS is produced only for diploid VCFs.
"""


__author__      = "Edgardo M. Ortiz"
__credits__     = "Juan D. Palacio-Mejía"
__version__     = "2.4"
__email__       = "e.ortiz.v@gmail.com"
__date__        = "2020-10-04"


import argparse
import gzip
import os
import random
import sys


# Dictionary of IUPAC ambiguities for nucleotides
# '*' is a deletion in GATK, deletions are ignored in consensus, lowercase consensus is udes when an
# 'N' or '*' is part of the genotype. Capitalization is used by some software but ignored by Geneious
# for example
ambiguities = {"*"    :"-", "A"    :"A", "C"    :"C", "G"    :"G", "N"    :"N", "T"     :"T",
               "*A"   :"a", "*C"   :"c", "*G"   :"g", "*N"   :"n", "*T"   :"t",
               "AC"   :"M", "AG"   :"R", "AN"   :"a", "AT"   :"W", "CG"   :"S",
               "CN"   :"c", "CT"   :"Y", "GN"   :"g", "GT"   :"K", "NT"   :"t",
               "*AC"  :"m", "*AG"  :"r", "*AN"  :"a", "*AT"  :"w", "*CG"  :"s",
               "*CN"  :"c", "*CT"  :"y", "*GN"  :"g", "*GT"  :"k", "*NT"  :"t",
               "ACG"  :"V", "ACN"  :"m", "ACT"  :"H", "AGN"  :"r", "AGT"  :"D",
               "ANT"  :"w", "CGN"  :"s", "CGT"  :"B", "CNT"  :"y", "GNT"  :"k",
               "*ACG" :"v", "*ACN" :"m", "*ACT" :"h", "*AGN" :"r", "*AGT" :"d",
               "*ANT" :"w", "*CGN" :"s", "*CGT" :"b", "*CNT" :"y", "*GNT" :"k",
               "ACGN" :"v", "ACGT" :"N", "ACNT" :"h", "AGNT" :"d", "CGNT" :"b",
               "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b", "*ACGNT":"N"}


# Dictionary for translating biallelic SNPs into SNAPP, only for diploid VCF
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


def extract_sample_names(vcf_file):
    """
    Extract sample names from VCF file
    """
    if vcf_file.endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    sample_names = []
    with opener(vcf_file, "rt") as vcf:
        for line in vcf:
            line = line.strip("\n")
            if line.startswith("#CHROM"):
                record = line.split("\t")
                sample_names = [record[i].replace("./", "") for i in range(9, len(record))]
                break
    return sample_names


def is_anomalous(record, num_samples):
    """
    Determine if the number of samples in current record corresponds to number of samples described
    in the line '#CHROM'
    """
    return bool(len(record) != num_samples + 9)


def is_snp(record):
    """
    Determine if current VCF record is a SNP (single nucleotide polymorphism) as opposed to MNP 
    (multinucleotide polymorphism)
    """
    return bool(len(record[3]) == 1 
                and len(record[4]) - record[4].count(",") == record[4].count(",") + 1)


def num_genotypes(record, num_samples):
    """
    Get number of genotypes in VCF record, total number of samples - missing genotypes
    """
    missing = 0
    for i in range(9, num_samples + 9):
        if record[i].startswith("."):
            missing += 1
    return num_samples - missing


def get_matrix_column(record, num_samples, resolve_IUPAC):
    """
    Transform a VCF record into a phylogenetic matrix column with nucleotides instead of numbers
    """
    nt_dict = {str(0): record[3].replace("-","*"), ".": "N"}
    alt = record[4].replace("-", "*")
    alt = alt.split(",")
    for n in range(len(alt)):
        nt_dict[str(n+1)] = alt[n]
    column = ""
    for i in range(9, num_samples + 9):
        genotype = record[i].split(":")[0].replace("/", "").replace("|", "")
        if resolve_IUPAC:
            column += nt_dict[random.choice(genotype)]
        else:
            column += ambiguities["".join(sorted(set([nt_dict[j] for j in genotype])))]
    return column


def get_matrix_column_bin(record, num_samples):
    """
    If VCF is diploid, return an alignment column in NEXUS binary from a VCF record
    """
    column = ""
    for i in range(9, num_samples + 9):
        genotype = record[i].split(":")[0]
        if len(genotype) == 3:
            column += gen_bin[genotype]
        else:
            column += "?"
    return column


def main():
    parser = argparse.ArgumentParser(description=__doc__, 
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--input",
        action = "store",
        dest = "filename",
        required = True,
        help = "Name of the input VCF file, can be gzipped")
    parser.add_argument("-m", "--min-samples-locus",
        action = "store",
        dest = "min_samples_locus",
        type = int,
        default = 4,
        help = "Minimum of samples required to be present at a locus (default=4)")
    parser.add_argument("-o", "--outgroup",
        action = "store",
        dest = "outgroup",
        default = "",
        help = "Name of the outgroup in the matrix. Sequence will be written as first taxon in the "
               "alignment.")
    parser.add_argument("-p", "--phylip-disable",
        action = "store_true",
        dest = "phylipdisable",
        help = "A PHYLIP matrix is written by default unless you enable this flag")
    parser.add_argument("-f", "--fasta",
        action = "store_true",
        dest = "fasta",
        help = "Write a FASTA matrix, disabled by default")
    parser.add_argument("-n", "--nexus",
        action = "store_true",
        dest = "nexus",
        help = "Write a NEXUS matrix, disabled by default")
    parser.add_argument("-b", "--nexus-binary",
        action = "store_true",
        dest = "nexusbin",
        help = "Write a binary NEXUS matrix for analysis of biallelic SNPs in SNAPP, only diploid "
               "genotypes will be processed, disabled by default.")
    parser.add_argument("-r", "--resolve-IUPAC",
        action = "store_true",
        dest = "resolve_IUPAC",
        help = "Randomly resolve heterozygous genotypes to avoid IUPAC ambiguities in the matrices")
    parser.add_argument("-v", "--version",
        action = "version",
        version = "%(prog)s {version}".format(version=__version__))
    args = parser.parse_args()


    filename = args.filename
    min_samples_locus = args.min_samples_locus
    outgroup = args.outgroup.split(",")[0].split(";")[0]
    phylipdisable = args.phylipdisable
    fasta = args.fasta
    nexus = args.nexus
    nexusbin = args.nexusbin
    resolve_IUPAC = args.resolve_IUPAC


    # Get samples names and number of samples in VCF
    sample_names = extract_sample_names(filename)
    num_samples = len(sample_names)
    if len(sample_names) == 0:
        print("\nSample names not found in VCF, your file may be corrupt or missing the header.\n")
        sys.exit()
    print("\nConverting file '{}':\n".format(filename))
    print("Number of samples in VCF: {:d}".format(len(sample_names)))

    # If the 'min_samples_locus' is larger than the actual number of samples in VCF readjust it
    min_samples_locus = min(num_samples, min_samples_locus)

    # Output filename will be the same as input file, indicating the minimum of samples specified
    if filename.endswith(".gz"):
        outfile = filename.replace(".vcf.gz",".min"+str(min_samples_locus))
    else:
        outfile = filename.replace(".vcf",".min"+str(min_samples_locus))
    # We need to create an intermediate file to hold the sequence data vertically and then transpose 
    # it to create the matrices
    if fasta or nexus or not phylipdisable:
        temporal = open(outfile+".tmp", "w")
    # If binary NEXUS is selected also create a separate temporal
    if nexusbin:
        temporalbin = open(outfile+".bin.tmp", "w")


    ##########################
    # PROCESS GENOTYPES IN VCF

    if filename.endswith(".gz"):
        opener = gzip.open
    else:
        opener = open

    with opener(filename, "rt") as vcf:
        # Initialize line counter
        snp_num = 0
        snp_accepted = 0
        snp_shallow = 0
        mnp_num = 0
        snp_biallelic = 0

        while 1:
            # Load large chunks of file into memory
            vcf_chunk = vcf.readlines(50000)
            if not vcf_chunk:
                break

            for line in vcf_chunk:
                line = line.strip()

                if line and not line.startswith("#"): # skip empty and commented lines
                    # Split line into columns
                    record = line.split("\t")
                    # Keep track of number of genotypes processed
                    snp_num += 1
                    # Print progress every 500000 lines
                    if snp_num % 500000 == 0:
                        print("{:d} genotypes processed.".format(snp_num))
                    if is_anomalous(record, num_samples):
                        print("Skipped potentially malformed line: {}".format(line))
                        continue
                    else:
                        # Check if the SNP has the minimum number of samples required
                        if num_genotypes(record, num_samples) < min_samples_locus:
                            # Keep track of loci rejected due to exceeded missing data
                            snp_shallow += 1
                            continue
                        else:
                            # Check that neither REF nor ALT contain MNPs
                            if is_snp(record):
                                # Add to running sum of accepted SNPs
                                snp_accepted += 1
                                # If nucleotide matrices are requested
                                if fasta or nexus or not phylipdisable:
                                    # Transform VCF record into an alignment column
                                    site_tmp = get_matrix_column(record, num_samples, resolve_IUPAC)
                                    # Uncomment for debugging
                                    # print(site_tmp)
                                    # Write entire row of single nucleotide genotypes to temp file
                                    temporal.write(site_tmp+"\n")
                                # Write binary NEXUS for SNAPP if requested
                                if nexusbin:
                                    # Check that the SNP only has two alleles
                                    if len(record[4]) == 1:
                                        # Add to running sum of biallelic SNPs
                                        snp_biallelic += 1
                                        # Translate genotype into 0 for homozygous REF, 1 for 
                                        # heterozygous, and 2 for homozygous ALT
                                        binsite_tmp = get_matrix_column_bin(record, num_samples)
                                        # Write entire row to temporary file
                                        temporalbin.write(binsite_tmp+"\n")
                            else:
                                # Keep track of loci rejected due to multinucleotide genotypes
                                mnp_num += 1

        # Print useful information about filtering of SNPs
        print("Total of genotypes processed: {:d}".format(snp_num))
        print("Genotypes excluded because they exceeded the amount "
              "of missing data allowed: {:d}".format(snp_shallow))
        print("Genotypes that passed missing data filter but were "
              "excluded for being MNPs: {:d}".format(mnp_num))
        print("SNPs that passed the filters: {:d}".format(snp_accepted))
        if nexusbin:
            print("Biallelic SNPs selected for binary NEXUS: {:d}".format(snp_biallelic))
        print("")

    if fasta or nexus or not phylipdisable:
        temporal.close()
    if nexusbin:
        temporalbin.close()


    #######################
    # WRITE OUTPUT MATRICES

    if not phylipdisable:
        output_phy = open(outfile+".phy", "w")
        output_phy.write("{:d} {:d}\n".format(len(sample_names), snp_accepted))

    if fasta:
        output_fas = open(outfile+".fasta", "w")

    if nexus:
        output_nex = open(outfile+".nexus", "w")
        output_nex.write("#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n\tFORMAT "
                         "DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n".format(len(sample_names),
                                                                                      snp_accepted))

    if nexusbin:
        output_nexbin = open(outfile+".bin.nexus", "w")
        output_nexbin.write("#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={:d} NCHAR={:d};\n\tFORMAT "
                            "DATATYPE=SNP MISSING=? GAP=- ;\nMATRIX\n".format(len(sample_names),
                                                                                     snp_biallelic))

    # Get length of longest sequence name
    len_longest_name = 0
    for name in sample_names:
        if len(name) > len_longest_name:
            len_longest_name = len(name)

    # Write outgroup as first sequence in alignment if the name is specified
    idx_outgroup = None
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
                print("Outgroup, '{}', added to the matrix(ces).".format(outgroup))

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
                print("Outgroup, '{}', added to the binary matrix.".format(outgroup))

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
                    print("Sample {:d} of {:d}, '{}', added to the nucleotide matrix(ces).".format(
                                                           s+1, len(sample_names), sample_names[s]))

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
                    print("Sample {:d} of {:d}, '{}', added to the binary matrix.".format(
                                                           s+1, len(sample_names), sample_names[s]))

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

    print( "\nDone!\n")

if __name__ == "__main__":
    main()
