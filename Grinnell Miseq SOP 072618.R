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

library(RColorBrewer)
library(phyloseq)
library(reshape2)
library(DESeq2)


darte_ed_16s <- import_mothur(mothur_shared_file = "stability.opti_mcc.shared", mothur_constaxonomy_file = "stability.cons.taxonomy")

metadata <- read.csv("GrinnellMiseqv2R.csv", row.names=1)
sample <- sample_data(metadata)
sample_data(darte_ed_16s)<- sample

#Rarefaction cruve is from jin
rf <- read.table("stability.opti_mcc.groups.rarefaction.txt", header=T)
rf[1:10,1:10]

melted <- melt(rf, id.vars=c("numsampled"))
head(melted)
pdf("rarefaction_curve.pdf", width=6, height=6)
ggplot(melted, aes(x=numsampled, y=value, group=variable))+
  geom_line()+
  theme_bw()
dev.off()
# set rarefying depth
after_remove_low_depth <- prune_samples(sample_sums(darte_ed_16s) >= 19000, darte_ed_16s)
sample_sums(after_remove_low_depth)
set.seed(1)
rare <- rarefy_even_depth(after_remove_low_depth, sample.size = 19000,rngseed=TRUE)
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
plot_bar(darte_ed_16s_filter_re_g, fill="Rank2", facet_grid= ".~Type") + scale_fill_hue() +scale_x_discrete(limits=c("AH110MT2"))

#Ed attempt at boxplot by phylum
tax <- tax_glom(darte_ed_16s_filter_re_g, taxrank="Rank2")
phylum10 = names(sort(taxa_sums(tax), TRUE)[1:10])
top10=prune_taxa(phylum10, tax)
dat <- psmelt(top10)

(boxplotbyphylum <- ggplot(dat, aes(Date_Type, Abundance))+
  facet_wrap(~Rank2, ncol=5, scales= "free") +
  geom_boxplot(aes(x = Date_Type, y = Abundance, fill= Type), outlier.shape = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),axis.title.x = element_blank()))
ggsave(boxplotbyphylum, filename = "BoxPlotbyPhylum.pdf", units = "in", width = 10, height = 7)

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
fin_table <- group_bar_chart(rare)
fin_table$group = factor(fin_table$group, levels = c("11/12/16 Manure", "11/08/16 Soil" ,"11/15/16 Manure Line","11/15/16 Soil", "2/17/17 Soil", "3/22/17 Soil" ))
(p <- ggplot(fin_table, aes(x=group,y=abundance, fill=phylum))+
    geom_bar(stat="identity",color="black")+labs(y="relative abundance")+ 
    scale_fill_hue()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank())+
    ylab("Relative Abundance"))
ggsave(p, filename = "RelativeStackedBar.pdf", units = "in", width = 6, height = 6)

tax_table(darte_ed_16s_filter_re)[1:10]
otu_table(darte_ed_16s_filter_re)[1:10]
otu_table(darte_ed_16s_filter)



#NMDS with manure
data.selected <- t(otu_table(rare))
mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$Date_Type <- as.factor(sample_data(rare) $Date_Type)
pairwise.adonis(data.selected, as.factor(sample_data(rare) $Date_Type))


(NMDS <- ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=Date_Type)) + geom_point(size=2) +
  geom_hline(yintercept=0.0, colour="grey", lty=2)+
  geom_vline(xintercept=0.0, colour="grey",lty=2) +
  scale_color_brewer(palette="Set1") +
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))
  +labs(color= "Sample"))
  #+stat_ellipse())
ggsave(NMDS, filename = "NMDSeverything.pdf", units = "in", width = 6, height = 6)

# show the NMDS coordination num 
data.sm

#NMDS without manure
remove_manure <-rare
sample_data(remove_manure) <- subset(sample_data(rare), sample_data(rare)$Date_Type != "11/12/16 Manure")

data.selected <- t(otu_table(remove_manure))
mds.all=metaMDS(data.selected,distance="bray", k=2)
data.sm=as.data.frame(scores(mds.all))
data.sm$Date_Type <- as.factor(sample_data(remove_manure) $Date_Type)



(NMDSNoMan<- ggplot(data=data.sm, aes(x=NMDS1, y=NMDS2, color=Date_Type)) + geom_point(size=2) +
  geom_hline(yintercept=0.0, colour="grey", lty=2)+
  geom_vline(xintercept=0.0, colour="grey",lty=2) +
  scale_color_brewer(palette="Set1")+
  theme_bw()+
  theme(legend.position = "bottom") +
  theme(legend.background = element_rect(fill="white", size=0.3, linetype="solid", colour="black"))+
  stat_ellipse()+
  labs(color="Sample"))
ggsave(NMDSNoMan, filename = "NMDSNoManure.pdf", units = "in", width = 8, height = 6)

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
p <- plot_richness(rare, x="Date_Type", measures=c("Shannon"))
(diversity <- ggplot(p$data, aes(x = Date_Type, y = value, fill = Type))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width=0.2)+
  ylab("Shannon Diversity Index")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title.x = element_blank())+
  labs(fill="Sample Type"))
ggsave(diversity, filename = "Diversity.pdf", units = "in", width = 8, height = 6)

divtest <- aov(value ~ Date_Type,p$data)
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
filtered <- filter_taxa(rare, function(x) sum(x) > 1, TRUE) 
#Get rid of small taxa and NTC
NoSmallTax <- filter_taxa(filtered, function(x) sum(x) > .001, TRUE) 
NoNTC <- prune_samples(!(rownames(sample_data(NoSmallTax)) %in% to_remove), NoSmallTax)


diagdds = phyloseq_to_deseq2(NoNTC, ~ Date_Type)

#manually calculate gm mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, test="Wald", fitType="parametric")

#results table
res0 = results(diagdds, contrast=c("Date_Type", "11/12/16 Manure", "11/08/16 Soil"), altHypothesis = "greater")
alpha = 0.05
sigtab0 = res0[which(res0$padj < alpha), ]
#sigtab0 = res0[which(res0$log2FoldChange > 0), ]
sigtab0 = cbind(as(sigtab0, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab0), ], "matrix"))
sigtab0 = sigtab0[c("log2FoldChange","pvalue","Rank2")]
sigtab0$Date_Type<-"11/12/16 Manure"

sigtabL = results(diagdds, contrast=c("Date_Type", "11/15/16 Manure Line", "11/08/16 Soil"), altHypothesis = "greater")
#sigtabL = sigtabL[which(sigtabL$padj < alpha), ]
sigtabL = cbind(as(sigtabL, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtabL), ], "matrix"))
#sigtabL= subset(sigtabL, rownames(sigtabL) %in% c(rownames(sigtab0)))
sigtabL = sigtabL[c("log2FoldChange","pvalue","Rank2")]
sigtabL$Date_Type<-"11/15/16 Manure Line"



sigtab1 = results(diagdds, contrast=c("Date_Type", "11/15/16 Soil", "11/08/16 Soil"), altHypothesis = "greater")
#sigtab1 = sigtab1[which(sigtab1$padj < alpha), ]
sigtab1 = cbind(as(sigtab1, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab1), ], "matrix"))
#sigtab1 = subset(sigtab1, rownames(sigtab1) %in% c(rownames(sigtab0)))
sigtab1 = sigtab1[c("log2FoldChange","pvalue","Rank2")]
sigtab1$Date_Type<-"11/15/16 Soil"

sigtab2 = results(diagdds, contrast=c("Date_Type", "2/17/17 Soil", "11/08/16 Soil"), altHypothesis = "greater")
sigtab2 = cbind(as(sigtab2, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab2), ], "matrix"))
sigtab2 = subset(sigtab2, rownames(sigtab2) %in% c(rownames(sigtab0)))
sigtab2 = sigtab2[c("log2FoldChange","pvalue","Rank2")]
sigtab2$Date_Type<-"2/7/17 Soil"

sigtab3 = results(diagdds, contrast=c("Date_Type", "3/22/17 Soil", "11/08/16 Soil"), altHypothesis = "greater")
sigtab3 = cbind(as(sigtab3, "data.frame"), as(tax_table(darte_ed_16s)[rownames(sigtab3), ], "matrix"))
sigtab3 = subset(sigtab3, rownames(sigtab3) %in% c(rownames(sigtab0)))
sigtab3 = sigtab3[c("log2FoldChange","pvalue","Rank2")]
sigtab3$Date_Type<-"3/22/17 Soil"

sigtab0 <- setNames(cbind(rownames(sigtab0), sigtab0, row.names = NULL), 
                    c("OTU","log2FoldChange","pval","Phyla","Date_Type"))
sigtabL <- setNames(cbind(rownames(sigtabL), sigtabL, row.names = NULL), 
                 c("OTU","log2FoldChange","pval","Phyla","Date_Type"))
sigtab1 <- setNames(cbind(rownames(sigtab1),  sigtab1, row.names = NULL), 
              c("OTU","log2FoldChange", "pval", "Phyla","Date_Type"))
sigtab2 <- setNames(cbind(rownames(sigtab2),  sigtab2, row.names = NULL), 
              c("OTU","log2FoldChange", "pval", "Phyla","Date_Type"))
sigtab3 <- setNames(cbind(rownames(sigtab3),  sigtab3, row.names = NULL), 
              c("OTU","log2FoldChange", "pval", "Phyla","Date_Type"))

rbind <- rbind(sigtabL, sigtab1, sigtab2, sigtab3)
noOTUbind <- rbind[c("log2FoldChange","Phyla","Date_Type", "pval")]

padj <- p.adjust(rbind$pval, "bonferroni", length(rbind$pval))


ggplot(rbind, aes(x=Date_Type, y=log2FoldChange))+
  geom_violin()+
  theme_bw()+
  facet_wrap(~Phyla, scales="free")
  
#Of sig. in manure and soil

sigL_M <- subset(sigtabL, sigtabL$OTU %in% sigtab0$OTU)

sig123<-rbind(sigtab1,sigtab2,sigtab3)

#sigL_Mtrack <- subset(sig123, sig123$OTU %in% sigL_M$OTU)

Mtrack<- subset(rbind, rbind$OTU %in% sigtab0$OTU)
#Mtrackfiltered <- subset(Mtrack, Mtrack$OTU %in% sigL_M$OTU)
Mtrackfiltered1<- Mtrackfiltered[which(Mtrack$pval < 0.05 & Mtrack$Date_Type=="3/22/17 Soil"),] 
Mtrackfiltered2 <- subset(Mtrackfiltered, Mtrack$OTU %in% Mtrackfiltered1$OTU)
                          
nosig <- Mtrackfiltered[which(Mtrackfiltered$pval > 0.05 & Mtrackfiltered$Date_Type=="3/22/17 Soil"), ]
elev <- Mtrackfiltered[which(Mtrackfiltered$pval < 0.05 & Mtrackfiltered$Date_Type=="3/22/17 Soil" & Mtrackfiltered$log2FoldChange>0), ]
dec<-Mtrackfiltered[which(Mtrackfiltered$pval < 0.05 & Mtrackfiltered$Date_Type=="3/22/17 Soil" & Mtrackfiltered$log2FoldChange<0), ]

MtrackPhyla<- aggregate(Mtrack, by=Mtrack[c("Phyla", "Date_Type")], FUN=mean, na.rm=TRUE)
names(MtrackPhyla) <- c("Phyla","Date_Type","x", "log2FoldChange","pval","xx","xxx")

#plot by OTU
ggplot(Mtrack, aes(x=Date_Type, y=log2FoldChange, group=OTU ,fill=Phyla))+
  theme_bw()+
  geom_line()+
  facet_wrap(~Phyla, scales = "free", nrow=2)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1, size=5))+
  ylab("Log2 Fold Change from 11/08/16 Abudance")+
  xlab("Sample")

#condense by phyla
(Mantrack <- (ggplot(Mtrack, aes(x=Date_Type, y=log2FoldChange, fill=Phyla))+
  theme_bw()+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~Phyla, scales = "free", nrow=2)+
  theme(axis.text.x = element_text(angle = 15, hjust = 1, vjust = 1, size=5), legend.position = "none", axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  ylab("Log2 Fold Change from 11/08/16 Soil Abundance")+
  xlab("Sample")))
ggsave(plot = last_plot(), filename = "ManureDispersionPlot.pdf", units = "in", width = 6, height = 6)

#stat test
manres<- aov(log2FoldChange~Phyla+Date_Type,Mtrack)
summary(manres)
TukeyHSD(manres)



#Bar graph summary of sig changes
# Phylum order
x = tapply(sigtab$log2FoldChange, sigtab$Rank2, function(x) max(x))
x = sort(x, TRUE)
sigtab$Rank2 = factor(as.character(sigtab$Rank2), levels=names(x))

testM = results(diagdds, contrast=c("Date_Type", "11/08/16 Soil", "11/12/16 Manure"), altHypothesis = "greater")
testM[1:51,]
y <- sigtab[c("log2FoldChange","Rank2")]
#Plot Phylum diff

ggplot(y, aes(x=Rank2, y=log2FoldChange))+geom_point()

#Lets see if about OTUs unique to manure
manureonly <- subset_samples(rare, Date_Type=="11/12/16 Manure")
premansoil <- subset_samples(rare, Date_Type=="11/08/16 Soil")
soil <- subset_samples(rare, Date_Type=="11/15/16 Manure Line")
 
no0man <- filter_taxa(manureonly, function(x) sum(x) > 1, TRUE)
no0soil <- filter_taxa(premansoil, function(x) sum(x) > 1, TRUE)
nol <- filter_taxa(soil, function(x) sum(x) > 1, TRUE)
uniqueman<- prune_taxa(!(row.names(otu_table(no0man)) %in% row.names(otu_table(no0soil))), no0man)

uniMtrack <- subset(rbind, rbind$OTU %in% row.names(otu_table(uniqueman)))
 

ggplot(uniMtrack, aes(x=Date_Type, y=log2FoldChange, group=OTU, color=Phyla))+
  theme_bw()+
  geom_line()


#export CSV
resOrdered <- res[order(res$pvalue),]
resSig <- subset(resOrdered, padj < 0.1)
write.csv(as.data.frame(sigtab), file="sig_results.csv")
resOrdered

#RDA

#import gene abundance data
RDAdata.gene <- read.csv("ARG_RDA.csv", row.names=1)
spe <- t(otu_table(darte_ed_16s_filter_re_g))
#env <- RDAdata.gene[, c('ermb','ermc')]  # select only two explanatory variables

trimmedRDA<- subset(RDAdata.gene, row.names(RDAdata.gene) %in% row.names(spe))
speRDA<- subset(spe, row.names(spe) %in% row.names(trimmedRDA))



spe.log <- log1p (speRDA)  # species data are in percentage scale which is strongly rightskewed, better to transform them
spe.hell <- decostand (spe.log, 'hell')  # we are planning to do tb-RDA, this is Hellinger pre-transformation
tbRDA <- rda (speRDA ~ ermB+ermC+strB+sul1+incW.repA+int1+Temperature, data = trimmedRDA)  # calculate tb-RDA with two explanatory variables
tbRDA


plot(tbRDA,display=c("sites","bp"),type="text",scaling=-1)

anova.cca(tbRDA, by="term")
anova.cca(tbRDA, by="margin")
anova.cca(tbRDA, by="axis")
permutest(tbRDA)

#RDA for indv species
otu <- spe.hell[, c('Otu00001')]
spRDA <- rda (otu ~ ermB+ermC+strB+sul1+incW.repA+int1, data = trimmedRDA)
spRDA
anova.cca(spRDA, by="term")
permutest(spRDA)

###SHADE LAB METHOD RDA
# Read OTU and map table from txt file
otu=t(otu_table(rare))
RDAdata.gene <- read.csv("ARG_RDA.csv", row.names=1)

#trim to same samples
trimmedRDA<- subset(RDAdata.gene, row.names(RDAdata.gene) %in% row.names(otu))
speRDA<- subset(spe, row.names(otu) %in% row.names(trimmedRDA))

#set RDP identifying OTUs later
rdp=speRDA[,ncol(speRDA)]

otu=speRDA[,-ncol(speRDA)]

#standardize data
map.s=decostand(trimmedRDA, method="standardize")

#build a CCA ordination, using Environmental factors
ord.ca=rda(otu ~ ermB+ermC+strB+sul1+incW.repA+int1+Temperature, data = map.s)
plot(ord.ca, display=(c("sites", "bp")))


# Score only display
scores(ord.ca, choices=c(1,2))
ord.scores=scores(ord.ca, choices=c(1,2))

# calculation of internal structure
str(ord.ca)

#extract axis scores (1,2) from CA plot
sites.sc=ord.scores$sites
sites.sc

# correlation between axis 1,2 & environmental factors
e=envfit(ord.ca, map.s)
e


#anova
anova.cca(ord.ca, by="term")

