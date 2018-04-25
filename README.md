# vcf2phylip
Convert SNPs in VCF format to PHYLIP for phylogenetic analysis

## _Brief description_
This script takes as input a VCF file and will use the SNP genotypes to create a matrix for phylogenetic analysis in the PHYLIP, FASTA, NEXUS, or binary NEXUS formats. For heterozygous SNPs the consensus is made and the IUPAC nucleotide ambiguity codes are written to the final matrix(ces). The code is optimized for large VCF matrices (hundreds of samples and millions of genotypes), for example, in our tests it processed a 20GB VCF (~3 million SNPs x 650 individuals) in ~27 minutes. The initial version of the script just produced a PHYLIP matrix but now we have added other popular formats, including the binary NEXUS file to run SNPs analysis with the SNAPP plugin in BEAST.

Additionally, you can choose a minimum number of samples per SNP to control the final amount of missing data. Since phylogenetic software usually root the trees at the first sequence in the alignment (e.g. RAxML, IQTREE, and MrBayes), the script also allows you to specify an OUTGROUP sequence that will be written in the first place in the alignment.

The script has been tested with VCF files produced by [*pyrad v.3.0.66*](https://github.com/dereneaton/pyrad), [*ipyrad v.0.7.x*](http://ipyrad.readthedocs.io/), [*Stacks v.1.47*](http://catchenlab.life.illinois.edu/stacks/), [*dDocent*](http://ddocent.com/), and [*GATK*](https://software.broadinstitute.org/gatk/).

Please don't hesitate to open an "Issue" if you find any problem.

## _Usage_
Just type `python vcf2phylip.py -h` to show the help of the program:

```
usage: vcf2phylip.py [-h] -i FILENAME [-m MIN_SAMPLES_LOCUS] [-o OUTGROUP]
                     [-p] [-f] [-n] [-b]

Converts SNPs in VCF format into an alignment for phylogenetic analysis

optional arguments:
  -h, --help            show this help message and exit
  -i FILENAME, --input FILENAME
                        Name of the input VCF file
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
_Example 1:_ Use default parameters to create a PHYLIP matrix with a minimum of 4 samples per SNP:
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
