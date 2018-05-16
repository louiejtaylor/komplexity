# Komplexity

A k-mer-based approach for masking low-complexity sequences. See also the [currently faster] sister implemetation in Rust: https://github.com/eclarke/komplexity

## Install

    git clone https://github.com/louiejtaylor/komplexity

or just [grab the binary](https://github.com/louiejtaylor/komplexity/blob/master/komplexity).

## Usage

    ./komplexity -in [-out] [-k] [-win] [-h]

    -in="test.fa": Input filename
    -k=4: k-mer size
    -out="test.fa_filtered": Output filename 
    -win=100: window length

Still under construction..for now works only on fasta
