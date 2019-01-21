library(ape)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(vegan)
library(RColorBrewer)
library(phyloseq)
library(DESeq2)



darte_ed_16s <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared", mothur_constaxonomy_file = "stability.cons.taxonomy")

metadata <- read.csv("GrinnellMiseqv2R.csv", row.names=1)
sample <- sample_data(metadata)
sample_data(darte_ed_16s)<- sample

#average number of reads
average <- mean(sample_sums(darte_ed_16s))

#Rarefaction curve is taken from online tutorial
rf <- read.table("stability.opti_mcc.groups.rarefaction.txt", header=T)
# Clean data file (no high or low confidence intervals)

# grep command may need to be changed. The idea is to pick out columns that are not hci or lci
rareDataColumnsPrimary <- rf[grepl("X0", names(rf))]

# this takes numsampled column from original data
numSampled <- rf[ ,c("numsampled")]

# this combines numsampled column with grepped columns
rareData <- cbind(numSampled,rareDataColumnsPrimary)

## Plot data with ggplot2

# transform data to long form
longRareData <- melt(rareData, id.vars = "numSampled")
head(longRareData)

# plotting
options(scipen = 100000)
(rarefaction <- ggplot(data = longRareData,
       aes(x = numSampled, y = value, group=variable)) +
  geom_line() +
  labs(x = "Number of Samples", y = "Number of OTUs") +
  theme_bw())
ggsave(rarefaction, filename = "rarefaction.pdf", units = "in", width = 6, height = 6, dpi=600)

# set rarefying depth
after_remove_low_depth <- prune_samples(sample_sums(darte_ed_16s) >= 19239, darte_ed_16s)
sample_sums(after_remove_low_depth)
set.seed(1)
rare <- rarefy_even_depth(after_remove_low_depth, sample.size = 19239,rngseed=TRUE)
sample_sums(rare)

#remove NTC
to_remove <- c("NTC")
#pruned <- prune_samples(!(rownames(sample_data(rare)) %in% to_remove), rare)

# filter out OTUs less than 10
#darte_ed_16s_filter <- filter_taxa(pruned, function(x) sum(x) > 10, TRUE) 

# relative abundance
darte_ed_16s_filter_re <- transform_sample_counts(rare, function(x) x /sum(x))

#Get rid of small taxa
#darte_ed_16s_filter2 <- filter_taxa(darte_ed_16s_filter_re, function(x) sum(x) > .001, TRUE) 

darte_ed_16s_filter_re_g = tax_glom(darte_ed_16s_filter_re, "Rank2")

#stacked bar by sample
plot_bar(darte_ed_16s_filter_re_g, fill="Rank2", facet_grid= ".~Type") + scale_fill_hue() +scale_x_discrete(limits=c("OT3S1R5G197"))

#Ed attempt at boxplot by phylum
tax <- tax_glom(darte_ed_16s_filter_re_g, taxrank="Rank2")
phylum10 = names(sort(taxa_sums(tax), TRUE)[1:10])
top10=prune_taxa(phylum10, tax)
dat <- psmelt(top10)

#test for each phyla
Acido <- subset(dat, dat$Rank2 == "Acidobacteria")
Actino <- subset(dat, dat$Rank2 == "Actinobacteria")
unclass <- subset(dat, dat$Rank2 == "Bacteria_unclassified")
Bactero <- subset(dat, dat$Rank2 == "Bacteroidetes")
cand <- subset(dat, dat$Rank2 == "candidate_division_WPS-1")
chloro <- subset(dat, dat$Rank2 == "Chloroflexi")
Firm <- subset(dat, dat$Rank2 == "Firmicutes")
Planc <- subset(dat, dat$Rank2 == "Planctomycetes")
Prote <- subset(dat, dat$Rank2 == "Proteobacteria")
Verr <- subset(dat, dat$Rank2 == "Verrucomicrobia")


summary(aov(Abundance~Date_Type, Verr))
TukeyHSD(aov(Abundance~Date_Type, unclass))
capture.output(TukeyHSD(aov(Abundance~Date_Type, Verr)),file="Reltukey.csv")

dat$Date_Type <- factor(dat$Date_Type, levels = c("manure","fall pre-manure soil","manure line","fall post-manure soil","spring 1 soil","spring 2 soil"))
#boxplot
(boxplotbyphylum <- ggplot(dat, aes(Date_Type, Abundance))+
  facet_wrap(~Rank2, ncol=5, scales= "free") +
  geom_boxplot(aes(x = Date_Type, y = Abundance, fill= Type), outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title.x = element_blank(),legend.position = "bottom", strip.text = element_text(size=8))+
  ylab(label = "Relative Abundance"))
ggsave(boxplotbyphylum, filename = "BoxPlotbyPhylum.eps", units = "in", width = 6.8, height = 6, dpi=600)




''#Jin's stacked bar
# You need to change "Type" into your name of treatment
group_bar_chart <- function(physeq){
  fin_table  = data.frame()
  for (group in unique(sample_data(physeq)$Date_Type) ){
    temp_physeq = prune_samples( (sample_data(physeq)$Date_Type == group), physeq)
    
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

bar.colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#aa6e28", "#f032e6", "#999999", "#8DD3C7", "#fffac8", "#BEBADA", "#fabebe", "#800000")
mycolors = c(brewer.pal(name="Set1", n = 9),brewer.pal(name="Set3", n = 5))

fin_table <- group_bar_chart(rare)
fin_table$group = factor(fin_table$group, levels = c("manure","fall pre-manure soil","manure line","fall post-manure soil","spring 1 soil","spring 2 soil"))
(p <- ggplot(fin_table, aes(x=group,y=abundance, fill=phylum))+
    geom_bar(stat="identity",color="black")+labs(y="relative abundance")+ 
    scale_fill_manual(values=bar.colors)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank())+
    ylab("Relative Abundance"))
ggsave(p, filename = "RelativeStackedBar.eps", units = "in", width = 6, height = 6, dpi=600)

tax_table(darte_ed_16s_filter_re)[1:10]
otu_table(darte_ed_16s_filter_re)[1:10]
otu_table(darte_ed_16s_filter)

storedtax <- tax_glom(darte_ed_16s_filter_re, taxrank="Rank2")
tax_table(storedtax)

#NMDS with manure
data.selected <- t(otu_table(rare))
mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$Date <- as.factor(sample_data(rare)$Date)

adonis(data.selected ~ as.factor(sample_data(rare)$Date) * as.factor(sample_data(rare)$Type))
adonis(data.selected ~ as.factor(sample_data(rare)$Date_Type))
pairwise.adonis(data.selected, as.factor(sample_data(rare)$Date_Type))
capture.output(pairwise.adonis(data.selected, as.factor(sample_data(rare) $Date_Type)), file="pairwiseadonis.csv")


(NMDS <- ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=Date_Type)) + geom_point(size=2) +
  geom_hline(yintercept=0.0, colour="grey", lty=2)+
  geom_vline(xintercept=0.0, colour="grey",lty=2) +
  scale_color_brewer(palette="Set1") +
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))
  +labs(color= "Sample"))
  #+stat_ellipse())
ggsave(NMDS, filename = "NMDSeverything.eps", units = "in", width = 6, height = 6, dpi=300)

# show the NMDS coordination num 
data.sm

#NMDS without manure
remove_manure <- rare
sample_data(remove_manure) <- subset(sample_data(rare), sample_data(rare)$Date_Type != "manure")

data.selected <- t(otu_table(remove_manure))
mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$Date_Type <- as.factor(sample_data(remove_manure)$Date_Type)

adonis(data.selected ~ as.factor(sample_data(remove_manure)$Date)*as.factor(sample_data(remove_manure)$Type))
pairwise.adonis(data.selected, as.factor(sample_data(remove_manure)$Date_Type))
scapture.output(pairwise.adonis(data.selected, as.factor(sample_data(rare) $Date_Type)), file="pairwiseadonis.csv")




(NMDSNoMan<- ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=Date_Type)) + geom_point(size=2) +
  geom_hline(yintercept=0.0, colour="grey", lty=2)+
  geom_vline(xintercept=0.0, colour="grey",lty=2) +
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))+
  stat_ellipse()+
  labs(color="Sample"))
ggsave(NMDSNoMan, filename = "NMDSNoManure(needsresizing).eps", units = "in", width = 8, height = 8, dpi =300)

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
p <- plot_richness(rare, x="Date_Type", measures=c("Shannon", "Chao1"))
p$Date_Type <- factor(p$Date_Type, limits = c("manure","fall pre-manure soil","manure line","fall post-manure soil","spring 1 soil","spring 2 soil"))
(diversity <- ggplot(p$data, aes(x = Date_Type, y = value, fill = Type))+
  scale_x_discrete(limits=c("manure","fall pre-manure soil","manure line","fall post-manure soil","spring 1 soil","spring 2 soil"))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2)+
  ylab("Diversity Metric")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank())+
  labs(fill="Sample Type")+
  facet_wrap(~variable, scales="free"))

chao1 <- subset(p$data, variable == "Chao1")
diversity + ylim(750, 4250)

shannon <- subset(p$data, variable == "Shannon")
diversity +  ylim(3.75, 7.25)

ggsave(diversity, filename = "Diversity.eps", units = "in", width = 6.8, height = 6, dpi =300)

pdata <- p$data
Chao1 <- subset(pdata, pdata$variable == "Chao1")
divtest <- aov(value ~ Date_Type, Chao1)
summary(divtest)
TukeyHSD(divtest)
capture.output(TukeyHSD(divtest),file="Chao1tukey.csv")

Shannon <- subset(pdata, pdata$variable == "Shannon")
divtest <- aov(value ~ Date_Type, Shannon)
summary(divtest)
TukeyHSD(divtest)
capture.output(TukeyHSD(divtest),file="Shannontukey.csv")

divprof<- renyi(t(otu_table(rare)), hill=TRUE)
plot(divprof)


# filter out OTUs less than 1
filtered <- filter_taxa(rare, function(x) sum(x) > 1, TRUE) 
#Get rid of small taxa and NTC
NoSmallTax <- filter_taxa(filtered, function(x) sum(x) > .001, TRUE) 
NoNTC <- prune_samples(!(rownames(sample_data(NoSmallTax)) %in% NoSmallTax), NoSmallTax)


diagdds = phyloseq_to_deseq2(NoNTC, ~ Date_Type)

#manually calculate gm mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

#results table
res0 = results(diagdds, contrast=c("Date_Type", "manure", "fall pre-manure soil"), altHypothesis = "greater")
alpha = 0.05
sigtab0 = res0
sigtab0 = res0[which(res0$padj < alpha), ]
#sigtab0 = res0[which(res0$log2FoldChange > 0), ]
sigtab0 = cbind(as(sigtab0, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab0), ], "matrix"))
sigtab0 = sigtab0[c("log2FoldChange","padj","Rank2","Rank3","Rank4")]
sigtab0$Date_Type<-"manure"

sigtabL = results(diagdds, contrast=c("Date_Type", "manure line", "fall pre-manure soil"), altHypothesis = "greater")
#sigtabL = sigtabL[which(sigtabL$padj < alpha), ]
sigtabL = cbind(as(sigtabL, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtabL), ], "matrix"))
#sigtabL= subset(sigtabL, rownames(sigtabL) %in% c(rownames(sigtab0)))
sigtabL = sigtabL[c("log2FoldChange","padj","Rank2","Rank3","Rank4")]
sigtabL$Date_Type<-"manure line"

sigtab1 = results(diagdds, contrast=c("Date_Type", "fall post-manure soil", "fall pre-manure soil"), altHypothesis = "greater")
#sigtab1 = sigtab1[which(sigtab1$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab1), ], "matrix"))
#sigtab1 = subset(sigtab1, rownames(sigtab1) %in% c(rownames(sigtab0)))
sigtab1 = sigtab1[c("log2FoldChange","padj","Rank2","Rank3","Rank4")]
sigtab1$Date_Type<-"fall post-manure soil"

sigtab2 = results(diagdds, contrast=c("Date_Type", "spring 1 soil", "fall pre-manure soil"), altHypothesis = "greater")
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab2), ], "matrix"))
#sigtab2 = subset(sigtab2, rownames(sigtab2) %in% c(rownames(sigtab0)))
sigtab2 = sigtab2[c("log2FoldChange","padj","Rank2","Rank3","Rank4")]
sigtab2$Date_Type<-"spring 1 soil"

sigtab3 = results(diagdds, contrast=c("Date_Type", "spring 2 soil", "fall pre-manure soil"), altHypothesis = "greater")
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab3), ], "matrix"))
#sigtab3 = subset(sigtab3, rownames(sigtab3) %in% c(rownames(sigtab0)))
sigtab3 = sigtab3[c("log2FoldChange","padj","Rank2", "Rank3","Rank4")]
sigtab3$Date_Type<-"spring 2 soil"

sigtab0 <- setNames(cbind(rownames(sigtab0), sigtab0, row.names = NULL), 
                    c("OTU","log2FoldChange","pval","Phyla","Class", "Order","Date_Type"))
sigtabL <- setNames(cbind(rownames(sigtabL), sigtabL, row.names = NULL), 
                    c("OTU","log2FoldChange","pval","Phyla","Class","Order", "Date_Type"))
sigtab1 <- setNames(cbind(rownames(sigtab1),  sigtab1, row.names = NULL),  
                    c("OTU","log2FoldChange","pval","Phyla","Class","Order", "Date_Type"))
sigtab2 <- setNames(cbind(rownames(sigtab2),  sigtab2, row.names = NULL), 
                    c("OTU","log2FoldChange","pval","Phyla","Class","Order", "Date_Type"))
sigtab3 <- setNames(cbind(rownames(sigtab3),  sigtab3, row.names = NULL), 
                    c("OTU","log2FoldChange","pval","Phyla","Class","Order", "Date_Type"))

rbind <- rbind(sigtabL, sigtab1, sigtab2, sigtab3)
noOTUbind <- rbind[c("log2FoldChange","Phyla", "Date_Type", "pval")]

padj <- p.adjust(rbind$pval, "bonferroni", length(rbind$pval))

Mtrack<- subset(rbind, rbind$OTU %in% sigtab0$OTU)

Mtrack$Date_Type <- factor(Mtrack$Date_Type, levels = c("manure line","fall post-manure soil","spring 1 soil","spring 2 soil"))

#condense by phyla
(Mantrack <- (ggplot(Mtrack, aes(x=Date_Type, y=log2FoldChange))+
  theme_bw()+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  facet_wrap(~Order, scales = "free", nrow=2)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1, size=5), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  ylab("Log2 Fold Change from Fall Pre-Manure Soil Abudance")+
  xlab("Sample")))
ggsave(plot = last_plot(), filename = "ManureONLYDispersionPlot.eps", units = "in", width = 12, height = 8)

#lanying edit
mtrackrelabu <- prune_taxa(taxa_names(darte_ed_16s_filter_re) %in% sigtab0$OTU, darte_ed_16s_filter_re)
dat <- psmelt(mtrackrelabu)
dat$Date_Type <- factor(dat$Date_Type, levels = c("manure","fall pre-manure soil","manure line","fall post-manure soil","spring 1 soil","spring 2 soil"))

(boxplotbyphylumMAN <- ggplot(dat, aes(Date_Type, Abundance))+
    facet_wrap(~Rank2, nrow=2, scales= "free") +
    geom_boxplot(aes(x = Date_Type, y = Abundance, fill= Rank2), outlier.shape = NA)+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1, size=5), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    ylab(label = "Relative Abundance"))
ggsave(boxplotbyphylum, filename = "BoxPlotbyPhylum.eps", units = "in", width = 6.8, height = 6, dpi=600)

datnoman <- subset(dat, dat$Date_Type!="manure")

(boxplotbyphylumMAN <- ggplot(datnoman, aes(Date_Type, Abundance))+
    facet_wrap(~Rank4, scales= "free") +
    geom_jitter()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1, size=5), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
    ylab(label = "Relative Abundance"))
ggsave(boxplotbyphylum, filename = "BoxPlotbyPhylum.eps", units = "in", width = 6.8, height = 6, dpi=600)



#stat test
manres<- aov(log2FoldChange~Phyla+Date_Type,Mtrack)
summary(manres)
TukeyHSD(manres)

co <- Mtrack[order(Mtrack$Class),]
po <- Mtrack[order(Mtrack$Phyla),]
Mtrack1 = Mtrack[which(Mtrack$Date_Type == "manure line" & Mtrack$pval < 0.05),]
sub <- Mtrack[c("OTU", "log2FoldChange", "pval","Date_Type")]
sub <- cbind(as(sub, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab0), ], "matrix"))
write.csv(sub, file = "MtrackOTU.csv")

#Unique Manure OTUs
manuresamples <- subset(sample_data(rare), sample_data(rare)$Date_Type=="manure")
premansamples <- subset(sample_data(rare), sample_data(rare)$Date_Type=="fall pre-manure soil")
test <- subset(sample_data(rare), sample_data(rare)$Date_Type=="fall post-manure soil")
testOTUs <- otu_table(rare)[,rownames(test)]

manureOTUs <- otu_table(rare)[,rownames(manuresamples)]
premanOTUs <- otu_table(rare)[,rownames(premansamples)]
filteredman <- filter_taxa(manureOTUs, function(x) sum(x) > 1, TRUE)
filteredpreman <- filter_taxa(premanOTUs, function(x) sum(x) > 1, TRUE)

manureonly <- subset(filteredman, !(rownames(filteredman) %in% row.names(filteredpreman)))
Mtab = subset(sigtab2, sigtab2$OTU %in% c(rownames(manureonly)))
Msigtab = subset(sigtab0, sigtab0$OTU %in% c(rownames(manureonly)))
shared <- subset(testOTUs, rownames(testOTUs) %in% rownames(manureonly))
filteredtest <- filter_taxa(shared, function(x) sum(x) > 1, TRUE)


(Mantrack <- (ggplot(Mtrack, aes(x=Date_Type, y=log2FoldChange))+
                theme_bw()+
                geom_boxplot(outlier.shape = NA)+
                theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1, size=5), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
                ylab("Log2 Fold Change from Fall Pre-Manure Soil Abudance")+
                xlab("Sample")))
ggsave(plot = last_plot(), filename = "OverallManureONLYDispersionPlot.eps", units = "in", width = 6.8, height = 6)
###Liz figure
ggplot(mtrackrelabu, aes(x=Date_Type, y=Abundance, fill=Rank2))+
  geom_bar()
plot_bar(mtrackrelabu, x="Date_Type", fill = "Rank2")
