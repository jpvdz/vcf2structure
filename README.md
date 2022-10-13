# vcf2structure
Convert variant data stored in VCF format to the STRUCTURE format. Requires a VCF file with bi-allelic SNPs and a population map, which is a CSV file containing sample IDs (first column) and population labels (second column). The columns need to be named `sample_id` and `population`.   

## Dependencies
Requires the `scikit-allel`, `numpy` and `pandas` packages. 

## How to run
The following example reads in a VCF file called `example.vcf` and a population map `example_popmap.csv`, converts the VCF into a STRUCTURE file called `example.str` with samples encoded over two rows (i.e., one row per sample is set to `False`):

`./vcf2structure.py example.vcf example_popmap.csv example.str False`

If you instead want samples encoded on a single row (and each locus encoded over two columns, with one allele per column) set the fourth argument to `True`. 