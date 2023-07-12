#! /bin/bash
runID=$2

#cd ${runID}/results/picard/hs_metrics/


grep -A 1 '^BAIT' $1 | awk 'BEGIN {{FS="\t"}} {{for (i=1;i<=NF;i++) {{a[NR"_"i]=$i}}}} END {{for (i=1;i<=NF;i++) {{print a[1"_"i]"\t"a[2"_"i]}}}}' #> $3



#for i in clean/*.fq.gz; do 
#  zcat $i | egrep -c '^@' >> clean/cleanReads.tsv; 
#done 
#
#
#for i in contam/*.fq.gz; do 
#  zcat $i | egrep -c '^@' >> contam/contamReads.tsv; 
#done 
#
#
#paste <(ls clean/*.fq.gz) clean/cleanReads.tsv contam/contamReads.tsv > readsCleanContam.tsv
#
#paste readsCleanContam.tsv <(awk '{ tmp= $3/($2+$2); printf"%0.5f\n", tmp }' readsCleanContam.tsv) > readsCleanContamPerc.tsv
#
#cat readsCleanContamPerc.tsv | tr ',' '.' > percReadsContam.tsv
#
#rm readsCleanContamPerc.tsv
#rm clean/cleanReads.tsv
#rm contam/contamReads.tsv
#rm readsCleanContam.tsv

#cd ../../../..
#readarray -t tar < <(cat ../../../../config/fusions.tsv)

#for i in ${tar[@]}; do 
#	x=$(echo $i | cut -f1 -d','); 
#	y=$(echo $i | cut -f2 -d','); 
#	echo $x $y; 
#	egrep $x fusions.tsv | egrep $y; 
#done

