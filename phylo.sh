[ -e $2 ] && rm $2
/project/haplotyping/code_NP/halotyper_upload/code/megacc -a param_UPGMA_nucleotide.mao -d $1 -o $2 -n
rm $1
