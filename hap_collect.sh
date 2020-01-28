cd $1


cat  <(awk '{print "@SQ\tSN:"$1"\tLN:"$2}' $4) <( cat  ./*/ranbow.single.hap.sam | sort -k3,3 -k4,4n) > 1.sam
#python2.7 correct_readname.py  > 1.sam
samtools view -Sb 1.sam > 1.bam
samtools sort 1.bam -o ranbow.single.hap.bam
samtools index ranbow.single.hap.bam
#rm 1.bam 1.sam

##cat  <(awk '{print "@SQ\tSN:"$1"\tLN:"$2}' $4) <( cat  ./*/ranbow.pair.hap.sam | sort -k3,3 -k4,4n) > 1.sam
##samtools view -Sb 1.sam > 1.bam
##samtools sort 1.bam ranbow.pair.hap
##samtools index ranbow.pair.hap.bam
##rm 1.bam 1.sam

cat  ./*/ranbow.single.hap > ranbow.single.hap
##cat  ./*/ranbow.pair.hap > ranbow.pair.hap

cat  ./*/ranbow.single.hap.fasta > ranbow.single.hap.fasta
##cat  ./*/ranbow.pair.hap.fasta   > ranbow.pair.hap.fasta


