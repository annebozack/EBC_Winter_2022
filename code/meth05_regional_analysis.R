#'# Regional DNA methylation analysis using DMRcate
#' Using data preprocessed in our script:  
#' meth01_process_data.R

options(warn=-1)
suppressMessages(library(data.table))
library(stringi)
suppressMessages(library(minfi))
options(warn=0)

#' load the data
load("data/processed.rda")

betas.clean = beta[manifest[manifest$probe_type=="cg" & !chr %in% c("X","Y")]$index,]

#'# Introduction to limma 
#' see [Smyth GK. Stat Appl Genet Mol Biol 2004](https://www.ncbi.nlm.nih.gov/pubmed/16646809).  
suppressMessages(library(limma,minfi))

#' First we need to define a model
model = model.matrix( ~smoker+sex+CD4+CD8+NK+B+MO+GR,data=pheno)
EWAS.limma <- eBayes(lmFit(betas.clean, design=model))
Top<-topTable(EWAS.limma, coef=2, number=Inf, sort.by="p")[1:10,]
Top

#' Bind results with annotation
require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
Annot<-as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
Annot.Tops<- Annot[match(rownames(Top),Annot$Name),]
Annot.Tops<-Annot.Tops[,c("UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island","chr","pos")]
Top<-cbind(Top[,1:5], Annot.Tops)

#' Order by chr and chromosomal position
Top$chr = as.numeric(gsub("chr", "", Top$chr))
Top[order(Top$chr,Top$pos),c(1,6,7,8,9,10)] 

#' Use of M-values reduces heteroscedasticity to meet linear model assumptions, see [Du P, et al. BMC Bioinformatics. 2010](https://pubmed.ncbi.nlm.nih.gov/21118553/). 
rm(Annot.Tops,beta,EWAS.limma,manifest,Top);gc()

#'# Introduction to differential variability analysis
#' see [Phipson and Oshlack. Genome Biol 2014](https://pubmed.ncbi.nlm.nih.gov/25245051/). 
suppressMessages(library(missMethyl))
suppressMessages(library(ChAMP))

#' Impute missing Beta-values (varFit will produce an error with missingness)
sum(is.na(betas.clean))
betas.impute = champ.impute(beta=betas.clean, pd=pheno, k=5, ProbeCutoff=0.2, SampleCutoff=0.1)
betas.impute = betas.impute$beta
gc()

#' coef parameter in varFit states which columns of design matrix correspond to the intercept and variable of interest
#' If Beta-values are used, a lofit transformation is performed within the varFit function
head(model)
EWAS.diffVar <- varFit(betas.impute, design=model, coef=c(1,2))
Top.diffVar <- topVar(EWAS.diffVar, coef=2, number=10, sort =TRUE)
Top.diffVar

#' Bind results with annotation
Annot.Top.diffVar<- Annot[match(rownames(Top.diffVar),Annot$Name),]
Annot.Top.diffVar<-Annot.Top.diffVar[,c("UCSC_RefGene_Name","UCSC_RefGene_Group","Relation_to_Island","chr","pos")]
Top.diffVar<-cbind(Top.diffVar, Annot.Top.diffVar)

#' Order by chr and chromosomal position
Top.diffVar$chr = as.numeric(gsub("chr", "", Top.diffVar$chr))
Top.diffVar[,c(2,6)] = round(Top.diffVar[,c(2,6)],2)
Top.diffVar[order(Top.diffVar$chr, Top.diffVar$pos),c(2,5,6,7,10,11)]

#' Although no CpGs meet FDR significance, we can see differences in variability among the top sites
cpg1 = data.frame(id = colnames(betas.impute), beta = betas.impute[rownames(betas.impute) == 'cg19754622',])
cpg1 = merge(cpg1, pheno[,c('gsm', 'smoker')], by.x = 'id', by.y = 'gsm') 
cpg1$smoker = as.numeric(factor(cpg1$smoker))

cpg2 = data.frame(id = colnames(betas.impute), beta = betas.impute[rownames(betas.impute) == 'cg21173402',])
cpg2 = merge(cpg2, pheno[,c('gsm', 'smoker')], by.x = 'id', by.y = 'gsm') 
cpg2$smoker = as.numeric(factor(cpg2$smoker))

boxplot(cpg1$beta ~ cpg1$smoker, col = c("blue", "red"), outline = F, xlab = 'smoking', ylab = 'Beta-value', names = c("non-smoker", "smoker"));points(jitter(cpg1$smoker, amount = 0.1), cpg1$beta, pch = 16)

boxplot(cpg2$beta ~ cpg2$smoker, col = c("blue", "red"), outline = F, xlab = 'smoking', ylab = 'Beta-value', names = c("non-smoker", "smoker"));points(jitter(cpg1$smoker, amount = 0.1), cpg2$beta, pch = 16)

#' diffVar results may be influenced by outliers
rm(Annot,Annot.Top.diffVar,Top.diffVar,EWAS.diffVar,betas.impute,cpg1,cpg2);gc()

#' Load package for regional analysis "DMRcate"
#' see [Peters et al. Bioinformatics 2015](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/1756-8935-8-6).  
#' Other popular options for conducting Regional DNA methylation analysis in R are Aclust and bumphunter 
suppressMessages(library(DMRcate)) # Popular package for regional DNA methylation analysis
#'Regions are now agglomerated from groups of significant probes 
#'Let's run the regional analysis using the Beta-values from our preprocessed data
myannotation <- cpg.annotate("array", na.omit(betas.clean), analysis.type="differential",arraytype="450K",
                             what="Beta",design=model, coef=2)

#'Regions are now agglomerated from groups of significant probes 
#'where the distance to the next consecutive probe is less than lambda nucleotides away
dmrcoutput.smoking <- dmrcate(myannotation, lambda=1000, C=2)

#'Let's look at the results
results.ranges <- extractRanges(dmrcoutput.smoking, genome = "hg19")
results.ranges

#' Plot the DMR using the Gviz

#' if you are interested in plotting genomic data the Gviz is extremely useful
#'Let's look at the second region
results.ranges[2]

# set up the grouping variables and colours
pheno$smoker<-as.factor(pheno$smoker)
cols = c("magenta","red")[pheno$smoker]
names(cols) = levels(pheno$smoker)[pheno$smoker]

#'Draw the plot a DMR\
#+ fig.width=8, fig.height=6, dpi=300
DMR.plot(ranges=results.ranges, dmr=1, CpGs=betas.clean, phen.col=cols, what = "Beta",
         arraytype = "450K", pch=16, toscale=TRUE, plotmedians=TRUE, 
         genome="hg19", samps=1:nrow(pheno))

#'Clean data
rm(dmrcoutput.smoking,myannotation,results.ranges);gc()


#'# Predicting smoking with EpiSmokEr
#' see [Bollepalli, Sailalitha, et al. Epigenomics 2019](https://pubmed.ncbi.nlm.nih.gov/31466478/).  
suppressMessages(require(EpiSmokEr))
# Make sure rows of pheno match betas column names
rownames(pheno)<-pheno$gsm
identical(colnames(betas.clean),rownames(pheno))

# pheno needs a column for sex,in the format of 1 and 2 representing men and women respectively
pheno$sex<-ifelse(pheno$sex=="m",1,2)
# 121 CpGs are used selected by LASSO along with Sex to get 3 categories (current, former and never smokers)
result <- epismoker(dataset=betas.clean, samplesheet = pheno, method = "SSt")
# Let's look how well the prediction performed
table(pheno$smoker,result$PredictedSmokingStatus)
rm(list = ls());gc()

#' End of script 05