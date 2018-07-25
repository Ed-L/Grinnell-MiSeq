library(ape)
library(ggplot2)
library(grid)
library(plyr)
library(reshape)
library(reshape2)
library(ggthemes)
library(gridExtra)
library(corrplot)
library(dplyr)
library(vegan)
library(SpadeR)
library(limma)

library(RColorBrewer)
library(phyloseq)

library(DESeq2)


darte_ed_16s <- import_mothur(mothur_shared_file = "~/Desktop/stability.opti_mcc.shared", mothur_constaxonomy_file = "~/Desktop/stability.cons.taxonomy")


metadata <- read.csv("~/Desktop/16slabels.csv", row.names=1)
sample <- sample_data(metadata)
sample_data(darte_ed_16s)<- sample

# set rarefying depth
after_remove_low_depth <- prune_samples(sample_sums(darte_ed_16s) >= 400, darte_ed_16s)
sample_sums(after_remove_low_depth)
set.seed(1)
rare <- rarefy_even_depth(after_remove_low_depth, sample.size = 400,rngseed=TRUE)
sample_sums(rare)

#remove NTC
to_remove <- c("NTC")
pruned <- prune_samples(!(rownames(sample_data(rare)) %in% to_remove), rare)

 # filter out OTUs less than 10
#darte_ed_16s_filter <- filter_taxa(pruned, function(x) sum(x) > 10, TRUE) 
 
 # relative abundance
darte_ed_16s_filter_re <- transform_sample_counts(pruned, function(x) x /sum(x))

#Get rid of small taxa
darte_ed_16s_filter2 <- filter_taxa(darte_ed_16s_filter_re, function(x) sum(x) > .001, TRUE) 

darte_ed_16s_filter_re_g = tax_glom(darte_ed_16s_filter2, "Rank2")

#stacked bar by sample
plot_bar(darte_ed_16s_filter_re_g, fill="Rank2", facet_grid= ".~Type") + scale_fill_hue() +scale_x_discrete(limits=c("AH110MT2"))

#Ed attempt at boxplot by phylum
tax <- tax_glom(darte_ed_16s_filter_re_g, taxrank="Rank2")
phylum10 = names(sort(taxa_sums(tax), TRUE)[1:10])
top10=prune_taxa(phylum10, tax)
dat <- psmelt(top10)
ggplot(dat, aes(Type, Abundance))  + facet_wrap(~Rank2, ncol=5, scales= "free") + geom_boxplot(aes(x = Type, y = Abundance, fill= Type))+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


''#Jin's stacked bar
# You need to change "Type" into your name of treatment
group_bar_chart <- function(physeq){
  fin_table  = data.frame()
  for (group in unique(sample_data(physeq)$DT) ){
    temp_physeq = prune_samples( (sample_data(physeq)$DT == group), physeq)
    
    new_table <- data.frame(taxa_sums(temp_physeq), tax_table(temp_physeq)[,2])
    colnames(new_table) = c("abundance", "phylum")
    sum_table = data.frame()
    for (x in unique(new_table$phylum) ){
      temp = subset(new_table, phylum==x)
      su = sum(temp$abundance)
      sum_table = rbind(sum_table, data.frame(su, temp[1,2]))
    }
    
    colnames(sum_table) = c("abundance", "phylum")
    perc_table = sum_table
    for (i in 1:nrow(sum_table)){
      perc_table[i,1] = sum_table[i,1] / sum(sum_table[,1])
    }
    
    other_table = data.frame()
    other = c()
    for (i in 1:nrow(perc_table)){
      if(perc_table[i,1] > 0.01) {
        other_table = rbind(other_table, perc_table[i,])
      }else{
        other = c(other, perc_table[i,1])
      }
    }
    
    sum(other)
    tep = data.frame(sum(other), "other")
    colnames(tep) = c("abundance", "phylum")
    tfin = rbind(other_table, tep)
    ttfin = cbind(tfin,group)
    fin_table = rbind(fin_table, ttfin)
    phy = unique(fin_table$phylum)
  }
  return(fin_table)
}


#stacked bar
fin_table <- group_bar_chart(pruned)
fin_table$group = factor(fin_table$group, levels = c("11/8/2016 Soil","11/12/2016 Manure" ,"11/15/2016 Manure Line","11/15/2016 Soil"))
(p <- ggplot(fin_table, aes(x=group,y=abundance, fill=phylum))+geom_bar(stat="identity",color="black")+labs(y="relative abundance")+ scale_fill_hue()+theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)))
tax_table(darte_ed_16s_filter_re)[1:10]
otu_table(darte_ed_16s_filter_re)[1:10]
otu_table(darte_ed_16s_filter)
ggsave(p, filename = "stackedbar.pdf", units = "in", width = 6, height = 4)


#NMDS
data.selected <- t(otu_table(pruned))
mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$DT <- as.factor(sample_data(rare) $DT)
pairwise.adonis(data.selected, as.factor(sample_data(rare) $DT))


ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=DT)) + geom_point(size=5) +
geom_hline(yintercept=0.0, colour="grey", lty=2)+
geom_vline(xintercept=0.0, colour="grey",lty=2) +
scale_color_brewer(palette="Set1") + theme(legend.position = c(0.11, 0.9)) +
theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))+stat_ellipse()+labs(color="Day")


# show the NMDS coordination num 
data.sm


#statistic method
pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='bonferroni')
{

co = combn(unique(as.character(factors)),2)
pairs = c()
F.Model =c()
R2 = c()
p.value = c()

for(elem in 1:ncol(co)){
if(sim.function == 'daisy'){
library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
} else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}

ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] , permutations=9999);
pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
F.Model =c(F.Model,ad$aov.tab[1,4]);
R2 = c(R2,ad$aov.tab[1,5]);
p.value = c(p.value,ad$aov.tab[1,6])
}
p.adjusted = p.adjust(p.value,method=p.adjust.m)
sig = c(rep('',length(p.adjusted)))
sig[p.adjusted <= 0.05] <-'.'
sig[p.adjusted <= 0.01] <-'*'
sig[p.adjusted <= 0.001] <-'**'
sig[p.adjusted <= 0.0001] <-'***'

pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
return(pairw.res)
}


#Diversity
p <- plot_richness(pruned, x="DT", measures=c("Shannon"))
ggplot(p$data, aes(x = DT, y = value, color = NULL))+geom_boxplot(alpha=0)+geom_jitter(width=0.2)

divtest <- aov(value ~ DT,p$data)
summary(divtest)
TukeyHSD(divtest)


#network
ig = make_network(pruned, type = "samples", distance = "bray", max.dist = 0.85)
plot_network(ig, pruned, color = "DT", line_weight = 0.4, label = NULL)


# prune to just the top 100 most abundant OTUs across all samples (crude).
GP100 = prune_taxa(names(sort(taxa_sums(pruned), TRUE))[1:100], pruned)
jg = make_network(GP100, "taxa", "jaccard", 0.3)
plot_network(jg, pruned, "taxa",color="Rank2", line_weight = 0.4, label = NULL)

#heatmap
gpt <- subset_taxa(pruned, Rank1=="Bacteria")
gpt <- prune_taxa(names(sort(taxa_sums(pruned),TRUE)[1:300]), gpt)
plot_heatmap(gpt, "NMDS", "bray", "Type", "Rank3", low="#66CCFF", high="#000033")

#Another heatmap way
heatmap(otu_table(gpt))


# filter out OTUs less than 1
filtered <- filter_taxa(darte_ed_16s, function(x) sum(x) > 1, TRUE) 
#Get rid of small taxa and NTC
NoSmallTax <- filter_taxa(filtered, function(x) sum(x) > .001, TRUE) 
NoNTC <- prune_samples(!(rownames(sample_data(NoSmallTax)) %in% to_remove), NoSmallTax)


diagdds = phyloseq_to_deseq2(NoNTC, ~ DT)

#manually calculate gm mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

#results table
res = results(diagdds, contrast=c("DT", "11/8/2016 Soil", "11/12/2016 Manure"))
alpha = 0.01
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab), ], "matrix"))
sigtab

sigman <- subset(sigtab, log2FoldChange > 0)
sigman <- sigman[c("log2FoldChange")]
sigman <- setNames(cbind(rownames(sigman), sigman, row.names = NULL), 
         c("OTU"))


res2 = results(diagdds, contrast=c("DT", "11/8/2016 Soil", "11/15/2016 Soil"))
sigman1 <- res[c("log2FoldChange")]
sigman1$Date <- "11/15/16"
sigman1 <- setNames(cbind(rownames(sigman1), sigman1, row.names = NULL), 
         c("OTU", "log2foldchange", "Date"))
sigman2 <- subset(sigman1, OTU %in% sigman$OTU)
sigman2 <- setNames(cbind(rownames(sigman2), sigman2, row.names = NULL), 
                    c("OTU", "log2foldchange", "Date"))
sigman2 <- as.data.frame(sigman2)
#may want to p.adjust()
#line graph
ggplot(sigman2, aes(x=Date, y=log2foldchange, color=OTU))+geom_line()+scale_color_hue()

#Bar graph summary of sig changes
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Rank2, function(x) max(x))
x = sort(x, TRUE)
sigtab$Rank2 = factor(as.character(sigtab$Rank2), levels=names(x))

# Genus order
#x = tapply(sigtab$log2FoldChange, sigtab$Genus, function(x) max(x))
#x = sort(x, TRUE)

y <- sigtab[c("log2FoldChange","Rank2")]
#Plot Phylum diff

ggplot(y, aes(x=Rank2, y=log2FoldChange))+geom_boxplot()

#export CSV
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(sigtab), file="sig_results.csv")
resOrdered

#RDA

#import gene abundance data
RDAdata.gene <- read.csv("~/Desktop/RDA r .csv", row.names=1)
env <- RDAdata.gene[, c('ermB.Abundance','ErmC.Abundance')]  # select only two explanatory variables


spe <- t(otu_table(darte_ed_16s_filter_re_g))


spe.log <- log1p (spe)  # species data are in percentage scale which is strongly rightskewed, better to transform them
spe.hell <- decostand (spe.log, 'hell')  # we are planning to do tb-RDA, this is Hellinger pre-transformation
tbRDA <- rda (spe.hell ~ ermB.Abundance+ErmC.Abundance, data = env)  # calculate tb-RDA with two explanatory variables
tbRDA


plot(tbRDA,display=c("species","bp"),type="text",scaling=-1)

anova.cca(tbRDA, by="term")
anova.cca(tbRDA, by="margin")
anova.cca(tbRDA, by="axis")
permutest(tbRDA)

#RDA for indv species
otu <- spe.hell[, c('Otu00001')]
spRDA <- rda (otu ~ ermb + ermc, data = env)
spRDA
anova.cca(spRDA, by="term")
permutest(spRDA)

