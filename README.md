# vcf2phylip
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2540861.svg)](https://doi.org/10.5281/zenodo.2540861)  
Convert SNPs in VCF format to PHYLIP, NEXUS, binary NEXUS, or FASTA alignments for phylogenetic analysis

## _Brief description_
This script takes as input a VCF file and will use the SNP genotypes to create a matrix for phylogenetic analysis in the PHYLIP (relaxed version), FASTA, NEXUS, or binary NEXUS formats. For heterozygous SNPs the consensus is made and the IUPAC nucleotide ambiguity codes are written to the final matrix(ces), any ploidy level is allowed and automatically detected. The code is optimized for large VCF matrices (hundreds of samples and millions of genotypes), for example, in our tests it processed a 20GB VCF (~3 million SNPs x 650 individuals) in ~27 minutes. The initial version of the script just produced a PHYLIP matrix but now we have added other popular formats, including the binary NEXUS file to run SNPs analysis with the SNAPP plugin in BEAST (only for diploid samples).

Additionally, you can choose a minimum number of samples per SNP to control the final amount of missing data. Since phylogenetic software usually root the trees at the first sequence in the alignment (e.g. RAxML, IQTREE, and MrBayes), the script also allows you to specify an OUTGROUP sequence that will be written in the first place in the alignment.

Compressed VCF files can be directly analyzed but the extension must be `.vcf.gz`.

The script has been tested with VCF files produced by [*pyrad v.3.0.66*](https://github.com/dereneaton/pyrad), [*ipyrad v.0.7.x*](http://ipyrad.readthedocs.io/), [*Stacks v.1.47*](http://catchenlab.life.illinois.edu/stacks/), [*dDocent*](http://ddocent.com/), [*GATK*](https://software.broadinstitute.org/gatk/), and [*freebayes*](https://github.com/ekg/freebayes).

Please don't hesitate to open an [`Issue`](https://github.com/edgardomortiz/vcf2phylip/issues) if you find any problem or suggestions for a new feature.

## _Usage_
Just type `python vcf2phylip.py -h` to show the help of the program:

```
usage: vcf2phylip.py [-h] -i FILENAME [-m MIN_SAMPLES_LOCUS] [-o OUTGROUP]
                     [-p] [-f] [-n] [-b]

The script converts a collection of SNPs in VCF format into a PHYLIP, FASTA,
NEXUS, or binary NEXUS file for phylogenetic analysis. The code is optimized
to process VCF files with sizes >1GB. For small VCF files the algorithm slows
down as the number of taxa increases (but is still fast).

Any ploidy is allowed, but binary NEXUS is produced only for diploid VCFs.

optional arguments:
  -h, --help            show this help message and exit
  -i FILENAME, --input FILENAME
                        Name of the input VCF file, can be gzipped
  -m MIN_SAMPLES_LOCUS, --min-samples-locus MIN_SAMPLES_LOCUS
                        Minimum of samples required to be present at a locus,
                        default=4 since is the minimum for phylogenetics.
  -o OUTGROUP, --outgroup OUTGROUP
                        Name of the outgroup in the matrix. Sequence will be
                        written as first taxon in the alignment.
  -p, --phylip-disable  A PHYLIP matrix is written by default unless you
                        enable this flag
  -f, --fasta           Write a FASTA matrix, disabled by default
  -n, --nexus           Write a NEXUS matrix, disabled by default
  -b, --nexus-binary    Write a binary NEXUS matrix for analysis of biallelic
                        SNPs in SNAPP, disabled by default
```

## _Examples_
In the following examples you can omit `python` if you change the permissions of `vcf2phylip.py` to executable.

_Example 1:_ Use default parameters to create a PHYLIP matrix with a minimum of 4 samples per  SNP:
```bash
python vcf2phylip.py --input myfile.vcf
# Which i equivalent to:
python vcf2phylip.py -i myfile.vcf
# This command will create a PHYLIP called myfile_min4.phy
```

_Example 2:_ Create a PHYLIP and a FASTA matrix using a minimum of 60 samples per SNP:
```bash
python vcf2phylip.py --input myfile.vcf --fasta --min-samples-locus 60
# Which is equivalent to:
python vcf2phylip.py -i myfile.vcf -f -m 60
# This command will create a PHYLIP called myfile_min60.phy and a FASTA called myfile_min60.fasta
```

_Example 3:_ Create all output formats, and select "sample1" as outgroup:
```bash
python vcf2phylip.py --input myfile.vcf --outgroup sample1 --fasta --nexus --nexus-binary
# Which is equivalent to:
python vcf2phylip.py -i myfile.vcf -o sample1 -f -n -b
# This command will create a PHYLIP called myfile_min4.phy, a FASTA called myfile_min4.fasta, a NEXUS called myfile_min4.nexus, and a binary NEXUS called myfile_min4.bin.nexus
```

_Example 4:_ If, for example, you wish to disable the creation of the PHYLIP matrix and only create a NEXUS matrix:
```bash
python vcf2phylip.py --input myfile.vcf --phylip-disable --nexus
# Which is equivalent to:
python vcf2phylip.py -i myfile.vcf -p -n
# This command will create only a NEXUS matrix called myfile_min4.nexus
```

## _Credits_
- Code: [Edgardo M. Ortiz](mailto:e.ortiz.v@gmail.com)
- Data and testing: [Juan D. Palacio-Mejía](mailto:jdpalacio@gmail.com)

## _Citation_
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2540861.svg)](https://doi.org/10.5281/zenodo.2540861)  
**Ortiz, E.M. 2019.** vcf2phylip v2.0: convert a VCF matrix into several matrix formats for phylogenetic analysis. DOI:10.5281/zenodo.2540861

