#!/bin/bash

#echo fasta
#../komplexity -fa test2.fa ; cat test2.fa_filtered
echo fastq
../komplexity -fq test.fq ; cat test.fq_filtered
