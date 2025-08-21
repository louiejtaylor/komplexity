# Komplexity

A k-mer-based approach for masking low-complexity sequences. See also the sister implemetation in Rust: https://github.com/eclarke/komplexity

## Install

    git clone https://github.com/louiejtaylor/komplexity

or just [grab the binary](https://github.com/louiejtaylor/komplexity/blob/master/komplexity).

## Usage

    ./komplexity [-fa] [-fq] [-out] [-k] [-win] [-h]

    -fa: Input filename (.fasta format)
    -fq: Input filename (.fastq format)
    -k=4: k-mer size
    -out: Output filename (default: appends "_masked" to the input filename) 
    -win=100: window length

## Known issues
 - Quality scores not preserved when masking [#5](https://github.com/louiejtaylor/komplexity/issues/5)
