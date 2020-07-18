
####
#import libraries
library(here)
library(tidyverse)
library(DESeq2)
library(tximport)
library(tximportData)
library(readr)
library(biomaRt)
library(tidymodels)
library(baguette)


###
file_path <- here("data/encode/k562")
anno_file<-here('data/Homo_sapiens.GRCh38.90.merge_transcriptome.fa.processed.txt')
files = Sys.glob(paste(file_path, "*.tsv.gz", sep="/"))
sample_types = c("cyto", "cyto", "nuc", "nuc")
output_file =  here("data/encode/k562/deseq_k562.txt") # to output K562 kallisto output
output_file_annotated = here("data/encode/k562/deseq_k562_annotated.txt") # the name of the output file

########################
# process files
########################

# import the annotation file
setwd(file_path)
gene_anno = read.table(anno_file, header=TRUE, colClasses="character")
tx2g = gene_anno[c('transcript', 'gene')]

# import the kallisto files and aggregate to gene level
# tx <- tximport(files, type='kallisto', txOut=TRUE, tx2gene=tx2g)
# $Abundance = TPM
# $counts = expected count
tx <- tximport(files, type='rsem', txOut=TRUE, tx2gene=tx2g)
gene <- summarizeToGene(tx, tx2g)

# prepare sample table for DESeq2
sample_table <- data.frame(condition=factor(sample_types))
rownames(sample_table) <- colnames(gene$counts)

# run DESeq2
dds <- DESeqDataSetFromTximport(gene, colData=sample_table, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

# save the output
counts_adj <- counts(dds, normalized=TRUE)
colnames(counts_adj) = paste(sample_table$condition, "counts_adj", sep="_")
final <- cbind(as.data.frame(res), counts_adj)
final <-rownames_to_column(final, var = "ensembl_gene_id_version")
write.table(final, output_file, sep="\t")


#deep biomart annotation
ensembl <- useEnsembl(biomart = "ensembl", dataset = 'hsapiens_gene_ensembl')
listAttributes(mart = ensembl)
annotation <- getBM(attributes = c("ensembl_gene_id_version", "percentage_gene_gc_content","gene_biotype","tmhmm"), mart=ensembl)
annotation<-annotation %>% as.tibble() %>% mutate_if(is.character,as.factor)
annotation %>% dplyr::count(gene_biotype)


final_annotated<-left_join(final[,1:7],annotation,by="ensembl_gene_id_version")

final_annotated<-final_annotated %>% 
  mutate(which_localization = case_when(
  log2FoldChange >0.75  ~ "cytosol", 
  log2FoldChange <0.75  ~ "nucleus", 
  TRUE ~ "non-localized")) %>% 
  mutate_if(is.character,as.factor) %>% 
  dplyr::select(which_localization,ensembl_gene_id_version,percentage_gene_gc_content,gene_biotype,tmhmm)

write.table(final_annotated, output_file_annotated, sep="\t")


final_annotated %>% 
  add_count(which_localization,gene_biotype) %>% 
  dplyr::select(which_localization,gene_biotype,n) %>% unique() %>% 
  drop_na %>% 
  ggplot(aes(which_localization,  n,color = gene_biotype
  )) +
  geom_boxplot(alpha = 0.2, size = 1.5, show.legend = TRUE) +
  labs(x = NULL, y = "localization")
ggsave(here("plots", "Gene_biotype_distribution.jpg"))



