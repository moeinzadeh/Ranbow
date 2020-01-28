cd $1
#for ((i = 1; i <= $2; i++)); do paste -d"\t" haplen.$i.txt  <( cat tmp$i.nwk_all | tr -d -c [='('=][=')']'\n' ) ; done |  awk -v x=$3 '{if($1>=x)print $2}' | sort | uniq -c
cat haplen.*.txt | cut -f 2 | sort | grep -E '^\(' | uniq -c 
#cat haplen.*.txt | awk '{if($6!="[]" && $2 ~ /^\(/ )print $2,$3,$6,$7,$8}'  > collect.txt
#cat haplen.*.txt  > collect_repeat.txt
cat haplen.*.txt  > collect.txt
