# https://kasperdanielhansen.github.io/genbioconductor/html/biomaRt.html


library(biomaRt) 
library(rhdf5)

mart <- useMart("ensembl")

nsmbl <- useDataset("hsapiens_gene_ensembl", mart)



# -------------------------------

s2 <- read.table("s2_genes.txt", header= F)$V1

exls <- list()
for(i in 1:length(s2)){
    ensg <- getBM(attributes= c("ensembl_gene_id", "external_gene_name", "ensembl_transcript_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"), filters= "external_gene_name", values= s2[i], mart= nsmbl)
	exls[[i]] <- split(ensg, ensg$ensembl_transcript_id)
	names(exls)[i] <- s2[i]
}

# ---------------

h5createFile("exons_ensembl_winters2018_Aug30_23.h5")


for(i in 1:length(exls)){
	#h5createGroup("tst.h5", names(exls)[i])
	if (length(exls[[i]]) != 0){
		for(j in 1:length(exls[[i]])){
			h5write(exls[[i]][[j]], "exons_ensembl_winters2018_Aug30_23.h5", names(exls[[i]])[j])
		}
	}
}

# print sessionInfo() to log file (and time stamnp for ensembl?) 

# --------------------------------------------------------------------

#listAttributes(nsmbl)
#listFilters(nsmbl)

# MET - ENSG00000105976.16
# EGFR - ENSG00000146648.21


# from https://support.bioconductor.org/p/56955/
#met <- getSequence(id="ENSG00000105976.16", type= "ensembl_gene_id", seqType="gene_exon", mart= nsmbl)
# or: 

#met <- getBM(attributes=c('ensembl_exon_id', "exon_chrom_start","exon_chrom_end","gene_exon"), filters = "ensembl_gene_id", values="ENSG00000105976", mart=ensembl, bmHeader=TRUE)
#names(met) <- c("end", "seq", "id", "start")
#egfr <- getBM(attributes=c('ensembl_exon_id', "exon_chrom_start","exon_chrom_end","gene_exon"), filters = 90, values="ENSG00000146648", mart=ensembl, bmHeader=TRUE)
#names(egfr) <- c("end", "seq", "id", "start")


# gene names as filters 
#"external_gene_name"

#abl1 <- getBM(attributes= c("ensembl_gene_id", "ensembl_exon_id", "exon_chrom_start", "exon_chrom_end"), filters= "external_gene_name", values= "ABL1", mart= nsmbl)



# COSMIC db  
#

#cosmicMart = useEnsembl(biomart="snp", dataset = hsapiens_snp_som")



