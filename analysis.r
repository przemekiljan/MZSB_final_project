library(phyloseq)
library(ggplot2)
library(vegan)

# Importing data
pre_otu<-read.table("mini_emp.tsv",sep="\t",header=T,row.names=1)
colnames(pre_otu) <- gsub("X","",colnames(pre_otu))
# otu <- as.matrix((pre_otu)+2.225074e-308) ## command that adds miniml value to each count
otu <- as.matrix(pre_otu)
abundance = otu_table(otu,taxa_are_rows = TRUE)
taxonomy <- as.matrix(read.table("Tax.csv",sep=",", header = T, row.names = 1))
tax_final = tax_table(taxonomy)
meta_data <- read.table('CAMDA_2019_EMP_metainformation.tsv',sep='\t', header=T,row.names=1)
row.names(meta_data) = tolower(row.names(meta_data))
meta_final <- sample_data(meta_data)
tree <- read_tree('97_otus.tree')
ds <- merge_phyloseq(abundance,tax_final,meta_final,tree)

#Rarefication step can be ommited, especially due to potential significance of low-count taxa
#ds.rarefied = rarefy_even_depth(ds, rngseed=2137, sample.size=0.1*min(sample_sums(ds)), replace=F)

#Merging samples by environmental biome
#ds.phylum = tax_glom(ds, taxrank = "Phylum", NArm = FALSE)
#ds.merged = merge_samples(ds.phylum, "envo_biome_2")
#ds.prop = transform_sample_counts(ds.merged, function(x) x/sum(x))

# plot_bar(ds.prop, fill="Phylum") + facet_wrap(~envo_biome_2, scales= "free_x", nrow=1)

#Merging samples by environmental biome in same studies
#Initialization
subset_1= subset_samples(ds, sample_data(ds)$"study_id"== sort(unique(sample_data(ds)$study_id))[1])
merged_1= merge_samples(subset_1, "envo_biome_2")
merged_1@sam_data$"envo_biome_2" = sample_names(merged_1)
sample_names(merged_1) = paste(sample_names(merged_1),sample_data(merged_1)$"study_id"[1], sep = "_")
phylo_fill = merged_1
#Propagation
for(i in sort(unique(sample_data(ds)$study_id))[2:length(sort(unique(sample_data(ds)$study_id)))]) {
  subset_i= subset_samples(ds, sample_data(ds)$"study_id"== i)
  merged_i= merge_samples(subset_i, "envo_biome_2")
  merged_i@sam_data$"envo_biome_2" = sample_names(merged_i)
  sample_names(merged_i) = paste(sample_names(merged_i),sample_data(merged_i)$"study_id"[1], sep = "_")
  phylo_fill = merge_phyloseq(phylo_fill, merged_i)
  
}
#Normalizng counts to table of proportions
phylo_fill_prop = transform_sample_counts(phylo_fill, function(x) x/sum(x))

#Measuring alpha diversity of samples
alph_rich = estimate_richness(ds, measures = c("Chao1","Shannon"))
pairwise.wilcox.test(alph_rich$Shannon, sample_data(ds)$envo_biome_2, p.adjust.method = 'hochberg')
pairwise.wilcox.test(alph_rich$Chao1, sample_data(ds)$envo_biome_2, p.adjust.method = 'hochberg')

#Beta diversity
bray_dist = phyloseq::distance(phylo_fill_prop, method="bray")
ordination = ordinate(phylo_fill_prop, method="PCoA", distance=bray_dist)
plot_ordination(phylo_fill_prop,
                ordination,
                color="envo_biome_2") 
                + theme(aspect.ratio=1) + geom_point(size=1) 
                + geom_text(aes(label=study_id),size=3,hjust=0, vjust=0)

#Multivariate differential abundance testing  
anosim(otu_table(phylo_fill_prop), sample_data(phylo_fill_prop)$envo_biome_2)
anosim(otu_table(phylo_fill_prop), sample_data(phylo_fill_prop)$study_id)

#It is also possible to subset samples by the sample feature, e.g. study_id
# canada = subset_samples(ds, sample_data(ds)$"study_id"==632)
# cd.rarefied = rarefy_even_depth(canada, rngseed=2137, sample.size=0.9*min(sample_sums(ds)), replace=F)
# cd.phylum = tax_glom(cd.rarefied, taxrank="Phylum", NArm=FALSE)
# 
# plot_bar(cd.phylum, fill="Phylum") + facet_wrap(~envo_biome_2, scales= "free_x", nrow=1)
#   
#Different subsetting
# schenzen = subset_samples(ds, sample_data(ds)$"study_id"==1521)
# schenzen.rarefied = rarefy_even_depth(schenzen, rngseed=2137, sample.size=0.9*min(sample_sums(ds)), replace=F)
# schenzen.phylum = tax_glom(cd.rarefied, taxrank="Phylum", NArm=FALSE)
# schenzen.merged = merge_samples(schenzen.phylum, "envo_biome_2")

# Plotting general abundance per environmental biome
# merged <- merge_samples(ds, 'envo_biome_2')
# plot_heatmap(merged, taxa.order = "Phylum")




