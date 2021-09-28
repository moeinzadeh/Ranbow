samtools view -H ../../mapping/allfinal.rmdup.bam > sel.sam
rm sel.fasta
grep "^#" ../../mapping/final.vcf > sel.vcf
while read p; do
	samtools faidx ../../refs/geneX10.merged.fasta $p >> sel.fasta
	samtools view ../../mapping/allfinal.rmdup.bam $p >> sel.sam
	grep $p ../../mapping/final.vcf >> sel.vcf
	echo $p
done <scaffolds.list
samtools view -Sb sel.sam > sel.bam
samtools index sel.bam
rm  sel.sam
