#! /bin/bash
runID=$1

#cd /home/schimar/smk/workflow/${runid}
cd $runID

id=${runID##*/}                 #$(echo $runID | cut -f1 -d'/')
# ------------------


mkdir -p stats/tmp/

rm stats/tmp/*

# ---------



egrep 'Pairs:' logs/bbmerge/*.log | cut -f2 >> stats/tmp/nRds.txt

egrep 'Joined:' logs/bbmerge/*.log >> stats/tmp/bbmerge_joined.txt

egrep 'Avg Insert:' logs/bbmerge/*.log >> stats/tmp/avg_insert_size.txt

egrep 'Standard Deviation' logs/bbmerge/*.log >> stats/tmp/avg_insert_size_sd.txt

egrep 'Insert range:' logs/bbmerge/*.log >> stats/tmp/insert_range.txt



# ------------



## picard hs metrics  ---  non-dedupe
cat results/fq2cons/picard/metrics/cons/*.metricsFrmt.tsv | egrep 'ON_BAIT_VS_SELECTED' | awk '{ print $2 }' > stats/tmp/hs_on_bait_vs_sel.txt

cat results/fq2cons/picard/metrics/cons/*.metricsFrmt.tsv | egrep 'PCT_TARGET_BASES_2X' | awk '{ print $2 }' > stats/tmp/hs_target2x.txt
cat results/fq2cons/picard/metrics/cons/*.metricsFrmt.tsv | egrep 'PCT_TARGET_BASES_10X' | awk '{ print $2 }' > stats/tmp/hs_target10x.txt
cat results/fq2cons/picard/metrics/cons/*.metricsFrmt.tsv | egrep 'PCT_TARGET_BASES_20X' | awk '{ print $2 }' > stats/tmp/hs_target20x.txt
cat results/fq2cons/picard/metrics/cons/*.metricsFrmt.tsv | egrep 'PCT_TARGET_BASES_30X' | awk '{ print $2 }' > stats/tmp/hs_target30x.txt

cat results/fq2cons/picard/metrics/cons/*.metricsFrmt.tsv | egrep 'MEDIAN_TARGET_COVERAGE' | awk '{ print $2 }' > stats/tmp/hs_med_target_cov.txt
cat results/fq2cons/picard/metrics/cons/*.metricsFrmt.tsv | egrep 'MEAN_TARGET_COVERAGE' | awk '{ print $2 }' > stats/tmp/hs_mean_target_cov.txt



# n raw reads  
#cat results/reads/trimmed/*.stats.txt | egrep '#Total' | awk '{ print $2/2 }' > stats/tmp/nTotRds.txt
#cat results/reads/trimmed/*.stats.txt | egrep '#Total' | awk '{ print $2 }' > stats/tmp/nTotRdsX2.txt

# n reads after trimming 
#cat logs/bbduk/*.log | egrep 'Result:' | cut -f2 -d$'\t' | cut -f1 -d' ' > stats/tmp/nTrmdRds.txt

# percent reads after trimming 
#cat logs/bbduk/*.log | egrep 'Result:' | cut -f2 -d$'\t' | egrep -o '[0-9.]+\%' | cut -f1 -d'%' > stats/tmp/percTrmdRds.txt

# percent mapped reads
#egrep '0 mapped \([0-9]+' results/mapped/*.fs.txt | egrep -o '[0-9.]+\%' | cut -f1 -d'%' > stats/tmp/percMpdRds.txt


#egrep '0 mapped \([0-9]+' results/mapped/dedupe/*.fs.txt | egrep -o '[0-9.]+\%' | cut -f1 -d'%' > stats/tmp/percMpdRdsDedupe.txt

# n mapped reads 
#cat results/mapped/*.fs.txt | egrep '0 mapped \([0-9]+' | awk '{ print int($1) }' > stats/tmp/nMpdRds.txt
#cat results/mapped/dedupe/*.fs.txt | egrep '0 mapped \([0-9]+' | awk '{ print int($1) }' > stats/tmp/nMpdRdsDedupe.txt
#
#
#
#
#cat logs/sambamba-markdup/*.log | egrep 'found [0-9]+ duplicates' | awk '{ print $2 }' > stats/tmp/nDups.txt
#
## median coverage global
#cat results/cov/dedupe/*.med.cov.txt > stats/tmp/med_cov.txt
#
## median target coverage 
#cat results/cov/dedupe/*.med.gcov.txt > stats/tmp/med_gcov.txt
#
## median coverage global (for reads before dedupe, NOTE: only different when markdup -r was used) 
#cat results/cov/*.med.cov.txt > stats/tmp/med_cov_nd.txt
#
## median target coverage (for reads before dedupe, NOTE: only different when markdup -r was used) 
#cat results/cov/*.med.gcov.txt > stats/tmp/med_nd.gcov.txt
#
#
#paste stats/tmp/nDups.txt stats/tmp/nTotRdsX2.txt | awk '{OFMT="%f"; print $1/$2 }' | tr ',' '.' > stats/tmp/percDup.txt
#
#awk '{OFMT="%.2f"; print $4 }' results/reads/percReadsContam.tsv | uniq > stats/tmp/percContam.txt
#
## picard hs metrics  ---  non-dedupe
#cat results/picard/hs_metrics/*.metrics.tsv | egrep 'ON_BAIT_VS_SELECTED' | awk '{ print $2 }' > stats/tmp/hs_on_bait_vs_sel.txt
#
#cat results/picard/hs_metrics/*.metrics.tsv | egrep 'PCT_TARGET_BASES_2X' | awk '{ print $2 }' > stats/tmp/hs_target2x.txt
#cat results/picard/hs_metrics/*.metrics.tsv | egrep 'PCT_TARGET_BASES_10X' | awk '{ print $2 }' > stats/tmp/hs_target10x.txt
#cat results/picard/hs_metrics/*.metrics.tsv | egrep 'PCT_TARGET_BASES_20X' | awk '{ print $2 }' > stats/tmp/hs_target20x.txt
#cat results/picard/hs_metrics/*.metrics.tsv | egrep 'PCT_TARGET_BASES_30X' | awk '{ print $2 }' > stats/tmp/hs_target30x.txt
#
#cat results/picard/hs_metrics/*.metrics.tsv | egrep 'MEDIAN_TARGET_COVERAGE' | awk '{ print $2 }' > stats/tmp/hs_med_target_cov.txt
#cat results/picard/hs_metrics/*.metrics.tsv | egrep 'MEAN_TARGET_COVERAGE' | awk '{ print $2 }' > stats/tmp/hs_mean_target_cov.txt
#
## picard hs metrics  ---  dedupe
#cat results/picard/hs_metrics/dedupe/*.metrics.tsv | egrep 'ON_BAIT_VS_SELECTED' | awk '{ print $2 }' > stats/tmp/hs_on_bait_vs_sel.ddp.txt
#
#cat results/picard/hs_metrics/dedupe/*.metrics.tsv | egrep 'PCT_TARGET_BASES_2X' | awk '{ print $2 }' > stats/tmp/hs_target2x.ddp.txt
#cat results/picard/hs_metrics/dedupe/*.metrics.tsv | egrep 'PCT_TARGET_BASES_10X' | awk '{ print $2 }' > stats/tmp/hs_target10x.ddp.txt
#cat results/picard/hs_metrics/dedupe/*.metrics.tsv | egrep 'PCT_TARGET_BASES_20X' | awk '{ print $2 }' > stats/tmp/hs_target20x.ddp.txt
#cat results/picard/hs_metrics/dedupe/*.metrics.tsv | egrep 'PCT_TARGET_BASES_30X' | awk '{ print $2 }' > stats/tmp/hs_target30x.ddp.txt
#
#cat results/picard/hs_metrics/dedupe/*.metrics.tsv | egrep 'MEDIAN_TARGET_COVERAGE' | awk '{ print $2 }' > stats/tmp/hs_med_target_cov.ddp.txt
#cat results/picard/hs_metrics/dedupe/*.metrics.tsv | egrep 'MEAN_TARGET_COVERAGE' | awk '{ print $2 }' > stats/tmp/hs_mean_target_cov.ddp.txt




# ----------------
#ls results/mapped/dedupe/*.bam | cut -f1 -d'.' | cut -f4 -d'/' > stats/tmp/idkeys.txt
ls results/bbmerge/*.fq.gz | cut -f3 -d'/' | cut -f1 -d'.' > stats/tmp/idkeys.txt 
 
cd stats/tmp/

#echo -e "id\tnTotalReads\tnTotRdsX2\tnRdsTrm\tpercRdsTrm\tpercContam\tnDups\tpercDups\tnMpdRds\tpercMpdRds\tnMpdRdsDedupe\tpercMappedRdsDedupe\tmed_nodedupe_gcov\tmed_nodedupe_covBed\tmed_gcov\tmed_covBed\ths_on_bait_vs_sel\ths_med_target_cov\ths_mean_target_cov\ths_target2x\ths_target10x\ths_target20x\ths_target30x\ths.dedupe_on_bait_vs_sel\ths.dedupe_med_target_cov\ths.dedupe_mean_target_cov\ths.dedupe_target2x\ths.dedupe_target10x\ths.dedupe_target20x\ths.dedupe_target30x" > ../runStats.${id}.tsv

echo -e "id\tnTotRds.txt\tmerge_joined\tavg_insert_size\tavg_insert_size_sd\tinsert_range" > ../runStats.${id}.tsv

#paste idkeys.txt nTotRds.txt nTotRdsX2.txt nTrmdRds.txt percTrmdRds.txt percContam.txt nDups.txt percDup.txt nMpdRds.txt percMpdRds.txt nMpdRdsDedupe.txt percMpdRdsDedupe.txt med_nd.gcov.txt med_cov_nd.txt med_gcov.txt med_cov.txt hs_on_bait_vs_sel.txt hs_med_target_cov.txt hs_mean_target_cov.txt hs_target2x.txt hs_target10x.txt hs_target20x.txt hs_target30x.txt hs_on_bait_vs_sel.ddp.txt hs_med_target_cov.ddp.txt hs_mean_target_cov.ddp.txt hs_target2x.ddp.txt hs_target10x.ddp.txt hs_target20x.ddp.txt hs_target30x.ddp.txt >> ../runStats.${id}.tsv
paste idkeys.txt nRds.txt <(cat bbmerge_joined.txt | cut -f3) <(cat avg_insert_size.txt | cut -f2) <(cat avg_insert_size_sd.txt | cut -f2) <(cat insert_range.txt | cut -f2) >> ../runStats.${id}.tsv
cd ../../../





## -------------------------------------------------
#for i in results/reads/clean/*.fq.gz; do 
#  zcat $i | egrep -c '^@' >> stats/tmp/cleanReads.tsv; 
#done 
#
#
#for i in results/reads/contam/*.fq.gz; do 
#  zcat $i | egrep -c '^@' >> contam/contamReads.tsv; 
#done 
#
#
#paste <(ls clean/*.fq.gz) clean/cleanReads.tsv contam/contamReads.tsv > readsCleanContam.tsv
#
#paste readsCleanContam.tsv <(awk '{ tmp= $3/($2+$3); printf"%0.5f\n", tmp }' readsCleanContam.tsv) > readsCleanContamPerc.tsv
#
#cat readsCleanContamPerc.tsv | tr ',' '.' > percReadsContam.tsv
#
#rm readsCleanContamPerc.tsv
#rm clean/cleanReads.tsv
#rm contam/contamReads.tsv
#rm readsCleanContam.tsv
#
#cd ../..
##readarray -t tar < <(cat ../../../../config/fusions.tsv)
#
##for i in ${tar[@]}; do 
##	x=$(echo $i | cut -f1 -d','); 
##	y=$(echo $i | cut -f2 -d','); 
##	echo $x $y; 
##	egrep $x fusions.tsv | egrep $y; 
##done
#
